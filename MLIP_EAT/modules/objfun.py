# objfun.py
from __future__ import annotations

import json
import numpy as np
import torch
from collections import Counter

from ase.data import atomic_numbers

from .projections import project_simplex_rows
from .function_helpers import entropy, softplus
from .lattice_helpers import ucparams_to_lattice, jacobian_ucparams_to_lattice, cell_to_ucparams
from .energy_above_hull import calc_hull, HullEnergyCalculator
from .formation_energy_helpers import FormationEnergyCalculator


class PBFGSPlateauReached(RuntimeError):
    """Raised by callback to early-stop PBFGS when cost plateaus."""
    pass


class objfun:
    """
    Provides:
      - fun_jac(z) -> (cost, grad)
      - pbfgs_callback(...)
      - history + last_metrics bookkeeping

    z layout:
      [ x_raw (n_atoms*n_elements),
        uc    (6),
        frac  (n_atoms*3) ]
    """

    def __init__(
        self,
        tau: float,
        gamma: float,
        mu: torch.Tensor,
        relax: bool,
        fix: bool,
        n_atoms: int,
        n_elements: int,
        Tsallis_q: float,
        selected_elements: list[str],
        atoms_base,
        calc,
        tol_cost: float,
        nIterations_tol: int,
        log,
        DB_JSON_PATH: str,
        REFERENCE_JSON_PATH: str,
        device: torch.device | None = None,
        dtype: torch.dtype = torch.float32,
    ):

        # Pick device consistently (prefer mu.device, else passed device, else cpu)
        if device is None:
            if torch.is_tensor(mu):
                device = mu.device
            else:
                device = torch.device("cpu")
        self.device = device
        self.dtype = dtype
        self.tau = float(tau)
        self.gamma = float(gamma)
        self.mu = mu
        self.relax = bool(relax)
        self.fix = bool(fix)

        self.Tsallis_q = float(Tsallis_q)
        self.selected_elements = list(selected_elements)

        self.atoms_base = atoms_base
        self.calc = calc
        self.log = log

        self.n_atoms = len(atoms_base)
        self.tol_cost = float(tol_cost)
        self.nIterations_tol = int(nIterations_tol)

        self.site_is_mix = torch.tensor(atoms_base.arrays["site_is_mix"], device=self.device).bool()
        self.site_fixed_Z = torch.tensor(atoms_base.arrays["site_fixed_Z"], device=self.device).long()

        self.mix_indices = torch.nonzero(self.site_is_mix, as_tuple=False).view(-1)
        self.fixed_indices = torch.nonzero(~self.site_is_mix, as_tuple=False).view(-1)
        self.n_mix = int(self.mix_indices.numel())
        self.n_elements = len(self.selected_elements)

        # IMPORTANT: x only exists for mix sites
        self.n_x = self.n_mix * self.n_elements
        self.n_uc = 6

        # ---- diagnostics / state
        self.iter = 1
        self.plateau_count = 0
        self.last_cost_for_plateau = None
        self.z_last = None  # last iterate seen by callback

        self.last_metrics = {
            "entropy_per_atom": None,
            "entropy_total": None,
            "cost": None,
            "total_energy": None,
            "formation_energy": None,
            "energy_above_hull": None,
            "cell_volume": None,
            "force_norm": None,
            "fmax": None,
        }

        self.history = {
            "iter": [],
            "entropy": [],
            "cost": [],
            "total_energy": [],
            "formation_energy": [],
            "energy_above_hull": [],
            "cell_volume": [],
            "force_norm": [],
            "fmax": [],
        }

        # ---- hull + formation energy calculators
        self.hullcalc = HullEnergyCalculator(DB_JSON_PATH, self.selected_elements)

        with open(REFERENCE_JSON_PATH, "r") as f:
            ref_energies = json.load(f)

        self.formation_calc = FormationEnergyCalculator(
            reference_energies_path=REFERENCE_JSON_PATH
        )

        # tensor in design-space order
        self.ref_energies_tensor = torch.tensor(
            [ref_energies[e] for e in self.selected_elements],
            dtype=self.dtype,
            device=self.device,
        )

        # after self.device / self.dtype are set, and after site_fixed_Z is available
        fixed_Z = torch.tensor(atoms_base.arrays["site_fixed_Z"], device=self.device, dtype=torch.long)
        fixed_Z = fixed_Z[fixed_Z > 0]  # only fixed sites

        # build a constant offset: sum_{fixed atoms} Eref(Z) / N
        fixed_ref_offset = 0.0
        if fixed_Z.numel() > 0:
            # map Z -> element symbol -> ref energy
            # easiest: use chemical symbols (robust and readable)
            fixed_symbols = [s for s, is_mix in zip(atoms_base.get_chemical_symbols(),
                                                    atoms_base.arrays["site_is_mix"])
                            if int(is_mix) == 0]
            # count fixed atoms by element
            counts = Counter(fixed_symbols)

            # scalar offset = sum_e (count_e/N)*Eref_e
            off = 0.0
            for el, cnt in counts.items():
                if el not in ref_energies:
                    raise ValueError(f"Reference energy missing for fixed element {el}")
                off += (cnt / float(len(atoms_base))) * float(ref_energies[el])

            fixed_ref_offset = torch.tensor(off, device=self.device, dtype=self.dtype)

        self.fixed_ref_offset = fixed_ref_offset


    # ------------------------------------------------------------------
    # Core: objective and gradient
    # ------------------------------------------------------------------
    def fun_jac(self, z_flat: torch.Tensor):
        # Ensure tensor on correct device/dtype
        if not torch.is_tensor(z_flat):
            z_flat = torch.tensor(z_flat, device=self.device, dtype=self.dtype)
        else:
            z_flat = z_flat.to(device=self.device, dtype=self.dtype)

        # fresh leaf tensors
        if self.fix:
            x_raw = z_flat.detach().clone().requires_grad_(True)
            atoms = self.atoms_base.copy()
            frac = atoms.get_scaled_positions()
            frac = torch.tensor(frac, device=self.device, dtype=self.dtype).detach().clone().requires_grad_(True)
            cell = atoms.get_cell()                 # ase.cell.Cell
            lattice = np.array(cell.array, dtype=float)  # shape (3,3)
            uc = cell_to_ucparams(lattice)
            uc = torch.tensor(uc, device=self.device, dtype=self.dtype).detach().clone().requires_grad_(True)
        else:
            x_raw = z_flat[: self.n_x].detach().clone().requires_grad_(True)
            uc = z_flat[self.n_x : self.n_x + self.n_uc].detach().clone().requires_grad_(True)
            frac = (
                z_flat[self.n_x + self.n_uc :]
                .detach()
                .clone()
                .view(self.n_atoms, 3)
                .requires_grad_(True)
            )

        # -------------------------
        # EAT weights (project → X)
        # -------------------------
        Xmix = project_simplex_rows(x_raw.view(self.n_mix, self.n_elements))

        entropy_total, grad_entropy = entropy(Xmix, self.Tsallis_q)  # total entropy
        grad_entropy_flat = grad_entropy.reshape(-1)

        entropy_cost, d_entropy_cost_d_entropy_total = softplus(
            entropy_total - self.n_atoms * self.tau, self.gamma
        )
        grad_entropy_cost = d_entropy_cost_d_entropy_total * grad_entropy_flat

        n_fairchem_elements = 100
        eat_fairchem = torch.zeros((self.n_atoms, n_fairchem_elements), device=self.device, dtype=Xmix.dtype)

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

        eat_fairchem.scatter_add_(1, cols, torch.zeros_like(cols, dtype=Xmix.dtype))  # (no-op; keeps shape known)
        # easiest: use advanced indexing add
        eat_fairchem[rows, cols] += Xmix

        # -------------------------
        # Lattice
        # -------------------------
        cell = ucparams_to_lattice(uc)  # torch (3,3)

        # -------------------------
        # FairChem energy + grads
        # -------------------------
        atoms = self.atoms_base.copy()
        if self.fix == False:
            atoms.set_scaled_positions(frac.detach().cpu().numpy())
            atoms.set_cell(cell.detach().cpu().numpy(), scale_atoms=True)

        E, F_cart, S = self.get_E_F_S(
            atoms,
            eat_weights_np=eat_fairchem.detach().cpu().numpy(),
            enable_eat_grad=True,
        )

        if "eat_grad" not in getattr(self.calc, "results", {}):
            raise RuntimeError(
                "FAIRChem did not return 'eat_grad'. "
                "Your fairchem build must support enable_eat_grad."
            )

        eat_grad_full = torch.tensor(self.calc.results["eat_grad"], device=self.device, dtype=Xmix.dtype)  # (n_atoms,100)

        # restrict to mix rows
        eat_grad_mix = eat_grad_full[self.mix_indices]                 # (n_mix,100)

        # gather only selected-element columns
        cols = fairchem_element_indices.unsqueeze(0).expand(self.n_mix, -1)      # (n_mix,n_elements)
        grad_Xmix = eat_grad_mix.gather(1, cols)                                 # (n_mix,n_elements)
        grad_x_fairchem = grad_Xmix.reshape(-1)

        # fractional gradient from forces: dE/d(frac) = -(F_cart @ cell^T)
        F_t = torch.tensor(F_cart, device=self.device, dtype=Xmix.dtype)  # (n_atoms,3)
        grad_frac_fairchem = -(F_t @ cell.T).reshape(-1)
        # stress -> dE/dH
        V = float(atoms.get_volume())
        S_t = torch.tensor(S, device=self.device, dtype=Xmix.dtype)  # (3,3)
        HinvT = torch.inverse(cell).T
        dE_dH = V * (HinvT @ S_t)

        # chain rule to ucparams
        J = jacobian_ucparams_to_lattice(uc)  # (3,3,6)
        grad_uc_fairchem = torch.zeros(6, device=self.device, dtype=Xmix.dtype)
        for k in range(6):
            grad_uc_fairchem[k] = torch.sum(dE_dH * J[:, :, k])

        E_t = torch.tensor(E, device=self.device, dtype=Xmix.dtype)

        # -------------------------
        # Formation energy (value + grads)
        # -------------------------
        comp_frac = Xmix.mean(dim=0)  # metals only, average over mix sublattice
        composition = {
            self.selected_elements[i]: float(comp_frac[i].item())
            for i in range(self.n_elements)
            if comp_frac[i].item() > 1e-8
        }

        N = float(self.n_atoms)
        invN = 1.0 / N

        # overall composition fractions contributed by mix sites:
        # c_mix_i = (sum over mix sites of Xmix[:,i]) / N
        c_mix = Xmix.sum(dim=0) * invN  # (n_elements,)

        # formation energy per atom including fixed atoms:
        # E_form = E_total/N - sum_i c_mix_i * Eref_i - fixed_ref_offset
        E_form = (E_t * invN) - torch.dot(c_mix.to(self.dtype), self.ref_energies_tensor.to(self.dtype)) - self.fixed_ref_offset


        E_hull, grad_x_hull = calc_hull(Xmix, self.hullcalc)  # E_hull: float/None, grad_x_hull: torch (n_x,)
        if E_hull is None or grad_x_hull is None:
            # if hull fails, degrade gracefully but keep optimization moving
            E_hull = float("nan")
            grad_x_hull = torch.zeros(self.n_x, device=self.device, dtype=Xmix.dtype)

        E_hull_t = torch.tensor(E_hull, device=self.device, dtype=Xmix.dtype)
        E_above_hull = E_form - E_hull_t

        invN = 1.0 / float(self.n_atoms)

        # d( E_total / N ) / d*
        grad_x_energy_part = invN * grad_x_fairchem
        grad_uc_energy_part = invN * grad_uc_fairchem
        grad_frac_energy_part = invN * grad_frac_fairchem

        # d/dXmix[a,i] of (- sum_i ( (sum_a Xmix[a,i])/N ) * Eref_i ) = -(1/N)*Eref_i
        grad_X_comp = (-invN) * self.ref_energies_tensor.to(dtype=Xmix.dtype).view(1, -1).expand(self.n_mix, -1)
        grad_x_comp = grad_X_comp.reshape(-1)


        # formation-energy gradients
        grad_x_form = grad_x_energy_part + grad_x_comp
        grad_uc_form = grad_uc_energy_part
        grad_frac_form = grad_frac_energy_part

        # energy-above-hull gradients
        grad_x_above_hull = grad_x_form - grad_x_hull
        _ = grad_x_above_hull  # available if you want to use it in cost later

        # -------------------------
        # Composition penalty
        # -------------------------
        penalty = 10 * 1/2 * torch.sum((comp_frac - self.mu)**2)
        grad_penalty = torch.zeros((self.n_mix, self.n_elements), device=self.device, dtype=Xmix.dtype)
        for i in range(self.n_elements):
            grad_penalty[:, i] = 10 * (comp_frac[i] - self.mu[i].to(device=self.device, dtype=Xmix.dtype)) / self.n_atoms
        grad_penalty = grad_penalty.reshape(-1)

        # -------------------------
        # Total cost (entropy + grand)
        # -------------------------
        cost = entropy_cost + E_form + penalty
        grad_x = grad_entropy_cost + grad_x_form + grad_penalty
        grad_uc = grad_uc_form
        grad_f = grad_frac_form

        # stage gating
        if self.relax:
            grad_x = grad_x * 0.0  # no EAT change during relax

        if self.fix == False:
            grad = torch.cat([grad_x, grad_uc, grad_f], dim=0)
        else:
            grad = grad_x

        if torch.isnan(grad).any() or torch.isinf(grad).any():
            raise RuntimeError("NaN/Inf detected in gradient.")

        # -------------------------
        # diagnostics
        # -------------------------
        with torch.no_grad():
            cell_volume = float(atoms.get_volume())

            gf = grad_f.view(self.n_atoms, 3)
            force_norm = float(torch.norm(gf, dim=1).mean().item())

            fmax = float(torch.norm(F_t, dim=1).max().item())

            self.last_metrics["cell_volume"] = cell_volume
            self.last_metrics["force_norm"] = force_norm
            self.last_metrics["fmax"] = fmax
            self.last_metrics["entropy_total"] = float(entropy_total.detach().item())
            self.last_metrics["entropy_per_atom"] = float((entropy_total / self.n_atoms).detach().item())
            self.last_metrics["cost"] = float(cost.detach().item())
            self.last_metrics["total_energy"] = float(E)
            self.last_metrics["formation_energy"] = float(E_form.detach().item())
            self.last_metrics["energy_above_hull"] = float(E_above_hull.detach().item())

        return cost.detach(), grad.detach()

    # ------------------------------------------------------------------
    # FairChem call
    # ------------------------------------------------------------------
    def get_E_F_S(self, atoms, eat_weights_np=None, enable_eat_grad: bool = False):
        if eat_weights_np is not None:
            atoms.new_array("eat_weights", np.asarray(eat_weights_np))
        if enable_eat_grad:
            atoms.info["enable_eat_grad"] = True
        else:
            atoms.info.pop("enable_eat_grad", None)

        atoms.calc = self.calc
        E = atoms.get_potential_energy()
        F = atoms.get_forces()
        S = atoms.get_stress(voigt=False)
        return E, F, S

    # ------------------------------------------------------------------
    # PBFGS callback (logging + plateau-stop)
    # ------------------------------------------------------------------
    def pbfgs_callback(self, z, f, err, k, converged, *args):
        self.z_last = z.detach().clone() if torch.is_tensor(z) else z

        f_val = f.item() if torch.is_tensor(f) else float(f)

        ent = self.last_metrics.get("entropy_per_atom", None)
        cost = self.last_metrics.get("cost", None)
        total_energy = self.last_metrics.get("total_energy", None)
        formation_energy = self.last_metrics.get("formation_energy", None)
        energy_above_hull = self.last_metrics.get("energy_above_hull", None)
        cell_volume = self.last_metrics.get("cell_volume", None)
        force_norm = self.last_metrics.get("force_norm", None)
        fmax = self.last_metrics.get("fmax", None)

        self.history["iter"].append(self.iter)
        self.history["entropy"].append(ent if ent is not None else np.nan)
        self.history["cost"].append(cost if cost is not None else np.nan)
        self.history["total_energy"].append(total_energy if total_energy is not None else np.nan)
        self.history["formation_energy"].append(formation_energy if formation_energy is not None else np.nan)
        self.history["energy_above_hull"].append(energy_above_hull if energy_above_hull is not None else np.nan)
        self.history["cell_volume"].append(cell_volume if cell_volume is not None else np.nan)
        self.history["force_norm"].append(force_norm if force_norm is not None else np.nan)
        self.history["fmax"].append(fmax if fmax is not None else np.nan)

        # Plateau detection on cost
        if cost is not None:
            if self.last_cost_for_plateau is not None:
                if abs(cost - self.last_cost_for_plateau) <= self.tol_cost:
                    self.plateau_count += 1
                else:
                    self.plateau_count = 0
            self.last_cost_for_plateau = cost

            if self.plateau_count >= self.nIterations_tol:
                self.log.info(
                    f"[PBFGS iter {k:03d}] Plateau convergence reached: "
                    f"|Δcost| <= {self.tol_cost} for {self.plateau_count} steps. "
                    f"Stopping PBFGS early for this anneal."
                )
                raise PBFGSPlateauReached()

        self.iter += 1

        if ent is None or cost is None:
            self.log.info(f"[PBFGS iter {k:03d}] f={f_val:.6f}, PG-error={err:.3e}")
            return

        self.log.info(
            f"[PBFGS iter {k:03d}] "
            f"cost={cost:.6f}, E_above_hull={energy_above_hull:.6f} eV/atom, "
            f"V={cell_volume:.2f} Å^3, fmax={fmax:.3e}, "
            f"ent/atom={ent:.3e}, PG-error={err:.3e}"
        )
