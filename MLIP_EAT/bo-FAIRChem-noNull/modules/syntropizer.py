import math
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
import torch
from ase.data import atomic_numbers

from .objfun import objfun  # <-- only import the class
from .projections import project_simplex_rows, project_flat_all
from .minimizer import minimizer
from .lattice_helpers import ucparams_to_lattice, cell_to_ucparams
from .formation_energy_helpers import FormationEnergyCalculator
from .energy_above_hull import calc_hull, HullEnergyCalculator
from .function_helpers import entropy, softplus


@dataclass
class OptimizationResult:
    history: Dict[str, List[float]]
    anneal_boundaries: List[int]
    final_energy_total_ev: float
    final_formation_energy_ev_per_atom: float
    final_energy_above_hull_ev_per_atom: float
    selected_elements: List[str]
    X_final: np.ndarray              # (n_atoms, n_elements)
    frac_coords: np.ndarray          # (n_atoms, 3)
    lattice_vectors: np.ndarray      # (3, 3)
    symbols: List[str]               # length n_atoms
    db_label: str


class syntropizer:
    def __init__(
        self,
        selected_elements,
        atoms_base,
        calc,
        log,
        mu,
        fix,
        DB_JSON_PATH,
        REFERENCE_JSON_PATH,
        n_anneals=1000,
        length_rel_delta=1.0,
        angle_deg_margin=40,
        angle_min_deg=40,
        angle_max_deg=140,
        tol_cost=1e-3,
        nIterations_tol=5,
        tau_mult=0.8,
        gamma=10,
        Tsallis_q=2.0,
        pbfgs_maxit_anneal=100,
        pbfgs_tol_anneal=1e-3,
        pbfgs_maxit_relax=500,
        pbfgs_tol_relax=5e-3,
        device=None,
    ):
        self.device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")

        self.selected_elements = list(selected_elements)
        self.atoms_base = atoms_base
        self.calc = calc
        self.log = log
        self.mu = mu.to(self.device) if torch.is_tensor(mu) else torch.tensor(mu, device=self.device, dtype=torch.float32)
        self.fix = fix

        self.DB_JSON_PATH = DB_JSON_PATH
        self.REFERENCE_JSON_PATH = REFERENCE_JSON_PATH

        self.site_is_mix = torch.tensor(atoms_base.arrays["site_is_mix"], device=self.device).bool()
        self.site_fixed_Z = torch.tensor(atoms_base.arrays["site_fixed_Z"], device=self.device).long()

        self.mix_indices = torch.nonzero(self.site_is_mix, as_tuple=False).view(-1)
        self.fixed_indices = torch.nonzero(~self.site_is_mix, as_tuple=False).view(-1)

        self.n_atoms = len(atoms_base)
        self.n_mix = int(self.mix_indices.numel())
        self.n_elements = len(self.selected_elements)

        # IMPORTANT: x only exists for mix sites
        self.n_x = self.n_mix * self.n_elements
        self.n_uc = 6

        self.n_anneals = int(n_anneals)
        self.length_rel_delta = float(length_rel_delta)
        self.tau_mult = float(tau_mult)
        self.gamma = float(gamma)
        self.Tsallis_q = float(Tsallis_q)

        self.tol_cost = float(tol_cost)
        self.nIterations_tol = int(nIterations_tol)

        self.pbfgs_maxit_anneal = int(pbfgs_maxit_anneal)
        self.pbfgs_tol_anneal = float(pbfgs_tol_anneal)
        self.pbfgs_maxit_relax = int(pbfgs_maxit_relax)
        self.pbfgs_tol_relax = float(pbfgs_tol_relax)

        # initial frac + ucparams
        frac0_np = atoms_base.get_scaled_positions()
        self.f0 = torch.tensor(frac0_np, dtype=torch.float32, device=self.device)

        cell0_np = atoms_base.cell.array.copy()
        uc0_np = cell_to_ucparams(cell0_np)
        self.uc0 = torch.tensor(uc0_np, dtype=torch.float32, device=self.device)

        # uc bounds
        angle_margin_rad = math.radians(angle_deg_margin)
        angle_min_rad = math.radians(angle_min_deg)
        angle_max_rad = math.radians(angle_max_deg)

        self.uc_lower = self.uc0.clone()
        self.uc_upper = self.uc0.clone()

        for i in range(3):
            self.uc_lower[i] = self.uc0[i] * (1.0 - self.length_rel_delta)
            self.uc_upper[i] = self.uc0[i] * (1.0 + self.length_rel_delta)

        for i in range(3, 6):
            self.uc_lower[i] = max(self.uc0[i].item() - angle_margin_rad, angle_min_rad)
            self.uc_upper[i] = min(self.uc0[i].item() + angle_margin_rad, angle_max_rad)

        self.uc_lower = self.uc_lower.to(self.device)
        self.uc_upper = self.uc_upper.to(self.device)

        # calculators
        self.formation_calc = FormationEnergyCalculator(reference_energies_path=self.REFERENCE_JSON_PATH)
        self.hullcalc = HullEnergyCalculator(self.DB_JSON_PATH, self.selected_elements)

        # tracking
        self.anneal_boundaries: List[int] = []
        self.history: Dict[str, List[float]] = {}  # will be filled from objfun
        self.global_iter = 0

        # optimizer state (persist between phases)
        self.z = None

    def _init_z(self):
        # x0 = 1/self.n_elements * torch.ones(self.n_mix, self.n_elements, device=self.device) + 1e-2 * torch.randn(self.n_mix, self.n_elements, device=self.device)
        dist = torch.distributions.Dirichlet(torch.ones(self.n_elements, device=self.device))
        x0 = dist.sample((self.n_mix,))
        x0 = project_simplex_rows(x0).view(-1)
        if self.fix:
            self.z = x0
        else:
            self.z = torch.cat([x0, self.uc0, self.f0.reshape(-1)], dim=0)

    def run(self) -> OptimizationResult:
        if self.z is None:
            self._init_z()

        self._syntropize_phase()
        return self._relax_and_finalize()

    def _syntropize_phase(self):
        tau_init = 1.0 - (1.0 / self.n_atoms)

        for anneal in range(self.n_anneals):
            self.anneal_boundaries.append(self.global_iter)

            tau = tau_init * (self.tau_mult ** anneal)

            self.log.info("")
            self.log.info("=" * 80)
            self.log.info(f"ANNEALING STEP {anneal + 1}/{self.n_anneals} | tau = {tau:.6e}")
            self.log.info("=" * 80)

            objfun_instance = objfun(
                tau=tau,
                gamma=self.gamma,
                mu=self.mu,
                relax=False,
                fix=self.fix,
                n_atoms=self.n_atoms,
                n_elements=self.n_elements,
                Tsallis_q=self.Tsallis_q,
                selected_elements=self.selected_elements,
                atoms_base=self.atoms_base,
                calc=self.calc,
                tol_cost=self.tol_cost,
                nIterations_tol=self.nIterations_tol,
                log=self.log,
                DB_JSON_PATH=self.DB_JSON_PATH,
                REFERENCE_JSON_PATH=self.REFERENCE_JSON_PATH,
            )

            opt = minimizer(
                fun_jac=objfun_instance.fun_jac,
                x0=self.z,
                args=(),
                project=project_flat_all,
                proj_args=(self.n_mix, self.n_elements, self.uc_lower, self.uc_upper, self.fix),
                callback=objfun_instance.pbfgs_callback,
                BFGS_hist=5,
                maxit=self.pbfgs_maxit_anneal,
                tol=self.pbfgs_tol_anneal,
                max_line_search=5,
            )

            try:
                _, self.z = opt.PBFGS()
            except Exception as e:
                # plateau exceptions should be handled inside objfun via z_last;
                # if anything else happens, prefer last known iterate
                self.log.warning(f"Anneal PBFGS interrupted: {e}")
                if getattr(objfun_instance, "z_last", None) is not None:
                    self.z = objfun_instance.z_last.detach().clone()

            # merge/append history (keep one source of truth)
            self.history = objfun_instance.history
            self.global_iter = int(self.history["iter"][-1]) if len(self.history.get("iter", [])) else self.global_iter

            # entropy-based rounding check
            Xmix = project_simplex_rows(self.z[:self.n_x].view(self.n_mix, self.n_elements))
            ent_total, _ = entropy(Xmix, self.Tsallis_q)
            ent = float((ent_total / self.n_atoms).detach().item())
            if ent < 1e-4:
                X_round = torch.zeros_like(Xmix)
                max_idx = torch.argmax(Xmix, dim=1)
                X_round.scatter_(1, max_idx.unsqueeze(1), 1.0)
                self.z = torch.cat([X_round.view(-1), self.z[self.n_x:]], dim=0)
                self.log.info("✅ Entropy small → rounding X and moving to final relaxation.")
                break

    def _relax_and_finalize(self) -> OptimizationResult:
        objfun_instance = objfun(
            tau=0.0,
            gamma=self.gamma,
            mu=self.mu,
            relax=True,
            fix=self.fix,
            n_atoms=self.n_atoms,
            n_elements=self.n_elements,
            Tsallis_q=self.Tsallis_q,
            selected_elements=self.selected_elements,
            atoms_base=self.atoms_base,
            calc=self.calc,
            tol_cost=self.tol_cost,
            nIterations_tol=self.nIterations_tol,
            log=self.log,
            DB_JSON_PATH=self.DB_JSON_PATH,
            REFERENCE_JSON_PATH=self.REFERENCE_JSON_PATH,
        )
        
        # Only relax if not fixed
        if self.fix == False:
            self.log.info("")
            self.log.info("=" * 80)
            self.log.info("🔧 FINAL PBFGS RELAXATION")
            self.log.info("=" * 80)

            opt = minimizer(
                fun_jac=objfun_instance.fun_jac,
                x0=self.z,
                args=(),
                project=project_flat_all,
                proj_args=(self.n_mix, self.n_elements, self.uc_lower, self.uc_upper, self.fix),
                callback=objfun_instance.pbfgs_callback,
                BFGS_hist=5,
                maxit=self.pbfgs_maxit_relax,
                tol=self.pbfgs_tol_relax,
                max_line_search=5,
            )

            try:
                _, self.z = opt.PBFGS()
            except Exception as e:
                self.log.warning(f"Final relaxation interrupted: {e}")
                if getattr(objfun_instance, "z_last", None) is not None:
                    self.z = objfun_instance.z_last.detach().clone()

            # unpack final state
            x_final = self.z[:self.n_x]
            uc_final = self.z[self.n_x:self.n_x + self.n_uc]
            frac_final = self.z[self.n_x + self.n_uc:].view(self.n_atoms, 3)
        else:
            x_final = self.z
            atoms = self.atoms_base.copy()
            frac_final = atoms.get_scaled_positions()
            frac_final = torch.tensor(frac_final, device=self.device)
            cell = atoms.get_cell()                 # ase.cell.Cell
            lattice = np.array(cell.array, dtype=float)  # shape (3,3)
            uc_final = cell_to_ucparams(lattice)
            uc_final = torch.tensor(uc_final, device=self.device).detach().clone()

        lattice = ucparams_to_lattice(uc_final).detach().cpu().numpy()
        print('final: ', lattice)
        frac_np = frac_final.detach().cpu().numpy()

        # Xf_mix is (n_mix,n_elements)
        Xf_mix = project_simplex_rows(x_final.view(self.n_mix, self.n_elements))
        Xf_mix_np = Xf_mix.detach().cpu().numpy()

        symbols = [None] * self.n_atoms

        # fixed symbols from atoms_base (already correct element labels like "O")
        base_symbols = list(self.atoms_base.get_chemical_symbols())
        for idx in range(self.n_atoms):
            if not bool(self.atoms_base.arrays["site_is_mix"][idx]):
                symbols[idx] = base_symbols[idx]

        # mix symbols from argmax
        elem_idx = np.argmax(Xf_mix_np, axis=1)
        mix_syms = [self.selected_elements[i] for i in elem_idx]
        for k, atom_idx in enumerate(self.mix_indices.detach().cpu().numpy().tolist()):
            symbols[atom_idx] = mix_syms[k]

        n_fairchem_elements = 100
        eat_fairchem = torch.zeros((self.n_atoms, n_fairchem_elements), device=self.device, dtype=Xf_mix.dtype)

        # fixed sites -> one-hot at atomic number
        # site_fixed_Z = 0 for mix sites, element Z for fixed sites
        fixed_idx = torch.nonzero(self.site_fixed_Z > 0, as_tuple=False).view(-1)
        if fixed_idx.numel() > 0:
            Z_fixed = self.site_fixed_Z[fixed_idx].clamp(min=0, max=n_fairchem_elements-1)
            eat_fairchem[fixed_idx, Z_fixed] = 1.0

        # mix sites -> scatter Xmix into selected element Z columns
        fairchem_element_indices = torch.tensor(
            [atomic_numbers[e] for e in self.selected_elements],
            dtype=torch.long,
            device=self.device,
        )
        cols = fairchem_element_indices.unsqueeze(0).expand(self.n_mix, -1)     # (n_mix,n_elements)
        rows = self.mix_indices.unsqueeze(1).expand(-1, self.n_elements)        # (n_mix,n_elements)

        eat_fairchem.scatter_add_(1, cols, torch.zeros_like(cols, dtype=Xf_mix.dtype))  # (no-op; keeps shape known)
        # easiest: use advanced indexing add
        eat_fairchem[rows, cols] += Xf_mix

        atoms_fin = self.atoms_base.copy()
        atoms_fin.set_scaled_positions(frac_np)
        atoms_fin.set_cell(lattice, scale_atoms=True)

        E_total, _, _ = objfun_instance.get_E_F_S(
            atoms_fin,
            eat_weights_np=eat_fairchem.detach().cpu().numpy(),
            enable_eat_grad=False,
        )

        with torch.no_grad():
            # -----------------------------
            # Composition for DB (mix sites only)
            # -----------------------------
            comp_frac_mix = Xf_mix.mean(dim=0)  # (n_elements,) averages over *mix sites only*
            composition = {
                self.selected_elements[i]: float(comp_frac_mix[i].item())
                for i in range(self.n_elements)
                if comp_frac_mix[i].item() > 1e-8
            }

            # -----------------------------
            # Formation energy per atom (FULL structure, includes fixed atoms)
            # -----------------------------
            # N = total atoms in the structure (mix + fixed)
            N = float(self.n_atoms)
            invN = 1.0 / N

            # Overall composition fractions contributed by mix sites, normalized by TOTAL atoms:
            # c_mix_i = (sum over mix sites of Xf_mix[:, i]) / N_total
            c_mix = Xf_mix.sum(dim=0) * invN  # (n_elements,)

            # E_form = E_total/N_total - sum_i c_mix_i * Eref_i - fixed_ref_offset
            # where fixed_ref_offset = sum_{fixed atoms} (Eref_fixed)/N_total  (constant)
            E_total_t = torch.tensor(E_total, device=self.device, dtype=Xf_mix.dtype)
            fixed_off = objfun_instance.fixed_ref_offset
            if not torch.is_tensor(fixed_off):
                fixed_off = torch.tensor(fixed_off, device=self.device, dtype=objfun_instance.ref_energies_tensor.dtype)
            else:
                fixed_off = fixed_off.to(dtype=objfun_instance.ref_energies_tensor.dtype)

            E_form_t = (
                E_total_t * invN
                - torch.dot(c_mix.to(objfun_instance.ref_energies_tensor.dtype), objfun_instance.ref_energies_tensor)
                - fixed_off
            )

            E_form = float(E_form_t.item())

            # -----------------------------
            # Hull (depends only on variable sublattice / mix composition)
            # -----------------------------
            E_hull, _ = calc_hull(Xf_mix, self.hullcalc)
            E_above = float(E_form - E_hull)


        # label for DB (you can replace with a stable hash later)
        label = " ".join([f"{el}{composition.get(el, 0.0):.3f}" for el in self.selected_elements])

        # expose history
        self.history = objfun_instance.history

        return OptimizationResult(
            history=self.history,
            anneal_boundaries=self.anneal_boundaries,
            final_energy_total_ev=float(E_total),
            final_formation_energy_ev_per_atom=float(E_form),
            final_energy_above_hull_ev_per_atom=float(E_above),
            selected_elements=self.selected_elements,
            X_final=Xf_mix_np,
            frac_coords=frac_np,
            lattice_vectors=lattice,
            symbols=symbols,
            db_label=label,
        )
