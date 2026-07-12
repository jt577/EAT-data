# objfun.py
from __future__ import annotations

import json
from collections import Counter

import numpy as np
import torch
from ase.data import atomic_numbers

from .energy_above_hull import HullEnergyCalculator, calc_hull
from .function_helpers import null_total_entropy
from .lattice_helpers import cell_to_ucparams, jacobian_ucparams_to_lattice, ucparams_to_lattice
from .projections import projector


class PBFGSPlateauReached(RuntimeError):
    """Raised by callback to early-stop PBFGS when cost plateaus."""


class objfun:
    """
    Provides:
      - fun_jac(z) -> (cost, grad)
      - pbfgs_callback(...)
      - history + last_metrics bookkeeping

    z layout:
      [ x_raw (n_mix*n_elements),
        uc    (6),
        frac  (n_atoms*3) ]
    """

    def __init__(
        self,
        tau: float,
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
        null_normalization_eps: float = 1e-8,
        composition_penalty_multiplier: float = 1.0,
        device: torch.device | None = None,
        dtype: torch.dtype = torch.float32,
    ):
        if device is None:
            if torch.is_tensor(mu):
                device = mu.device
            else:
                device = torch.device("cpu")

        self.device = device
        self.dtype = dtype
        self.tau = float(tau)
        self.mu = mu
        self.relax = bool(relax)
        self.fix = bool(fix)
        self.Tsallis_q = float(Tsallis_q)
        self.selected_elements = list(selected_elements)
        self.atoms_base = atoms_base
        self.calc = calc
        self.log = log
        self.null_normalization_eps = float(null_normalization_eps)
        self.composition_penalty_multiplier = float(composition_penalty_multiplier)

        self.n_atoms = len(atoms_base)
        self.tol_cost = float(tol_cost)
        self.nIterations_tol = int(nIterations_tol)

        self.site_is_mix = torch.tensor(atoms_base.arrays["site_is_mix"], device=self.device).bool()
        self.site_is_null = torch.tensor(
            atoms_base.arrays.get("site_is_null", np.zeros(len(atoms_base), dtype=np.int32)),
            device=self.device,
        ).bool()
        self.site_fixed_Z = torch.tensor(atoms_base.arrays["site_fixed_Z"], device=self.device).long()

        self.mix_indices = torch.nonzero(self.site_is_mix, as_tuple=False).view(-1)
        self.fixed_indices = torch.nonzero(~self.site_is_mix, as_tuple=False).view(-1)
        self.mix_is_null = self.site_is_null[self.mix_indices]

        self.n_mix = int(self.mix_indices.numel())
        self.n_null_mix = int(self.mix_is_null.sum().item())
        self.n_elements = len(self.selected_elements)
        self.n_x = self.n_mix * self.n_elements
        self.n_uc = 6

        self.iter = 1
        self.energy_evaluations = 0
        self.plateau_count = 0
        self.last_cost_for_plateau = None
        self.z_last = None
        self.last_composition_vector = None

        self.last_metrics = {
            "energy_evaluations": None,
            "entropy_per_atom": None,
            "entropy_total": None,
            "species_entropy_total": None,
            "null_entropy_total": None,
            "cost": None,
            "total_energy": None,
            "formation_energy": None,
            "energy_above_hull": None,
            "cell_volume": None,
            "force_norm": None,
            "fmax": None,
            "discreteness_dist": None,
            "mix_total_occupancy": None,
            "total_atom_count": None,
            "composition_penalty": None,
            "objective_no_entropy": None,
        }

        self.history = {
            "iter": [],
            "energy_evaluations": [],
            "entropy": [],
            "cost": [],
            "total_energy": [],
            "formation_energy": [],
            "energy_above_hull": [],
            "cell_volume": [],
            "force_norm": [],
            "fmax": [],
            "discreteness_dist": [],
            "mix_total_occupancy": [],
            "total_atom_count": [],
            "composition_penalty": [],
            "objective_no_entropy": [],
        }

        self.all_history = {k: [] for k in self.last_metrics}
        self.all_history["iter"] = []
        self.all_history["frac_coords"] = []
        self.all_history["cart_positions_angstrom"] = []
        self.all_history["lattice_vectors_angstrom"] = []
        self.all_history["species_weights"] = []
        self.all_history["site_occupancies"] = []
        self.all_history["vacuum_probabilities"] = []
        self.all_history["dominant_species"] = []
        self.all_history["trajectory"] = []

        self.base_symbols = list(self.atoms_base.get_chemical_symbols())
        self.species_order = list(self.selected_elements)
        for atom_idx in self.fixed_indices.detach().cpu().numpy().tolist():
            symbol = self.base_symbols[atom_idx]
            if symbol not in self.species_order:
                self.species_order.append(symbol)
        self.species_to_index = {symbol: i for i, symbol in enumerate(self.species_order)}
        self.mix_species_cols = np.array([self.species_to_index[symbol] for symbol in self.selected_elements], dtype=int)

        self.hullcalc = HullEnergyCalculator(DB_JSON_PATH, self.selected_elements)

        with open(REFERENCE_JSON_PATH, "r", encoding="utf-8") as f:
            ref_energies = json.load(f)

        self.ref_energies_tensor = torch.tensor(
            [ref_energies[e] for e in self.selected_elements],
            dtype=self.dtype,
            device=self.device,
        )

        fixed_symbols = [self.base_symbols[i] for i in self.fixed_indices.detach().cpu().numpy().tolist()]
        self.fixed_atom_count = float(len(fixed_symbols))
        self.fixed_ref_energy_sum = torch.tensor(
            sum(float(ref_energies[sym]) for sym in fixed_symbols),
            dtype=self.dtype,
            device=self.device,
        )
        fixed_selected_counts = Counter(sym for sym in fixed_symbols if sym in self.selected_elements)
        self.fixed_selected_counts = torch.tensor(
            [float(fixed_selected_counts.get(el, 0.0)) for el in self.selected_elements],
            dtype=self.dtype,
            device=self.device,
        )

        self.projector = projector(
            self.n_mix,
            self.n_elements,
            True,
            self.mu,
            None,
            None,
            mix_is_null=self.mix_is_null,
        )

    def _selected_counts_and_composition(self, Xmix: torch.Tensor):
        counts = Xmix.sum(dim=0) + self.fixed_selected_counts.to(dtype=Xmix.dtype)
        total = torch.sum(counts)
        total_safe = torch.clamp(total, min=self.null_normalization_eps)
        composition = counts / total_safe
        return counts, total, total_safe, composition

    def _lift_composition_gradient(self, grad_c: torch.Tensor, composition: torch.Tensor, total_safe: torch.Tensor):
        grad_rows = grad_c.view(1, -1).expand(self.n_mix, -1) / total_safe
        if self.n_null_mix > 0:
            offset = torch.dot(grad_c, composition) / total_safe
            grad_rows = grad_rows.clone()
            grad_rows[self.mix_is_null] -= offset
        return grad_rows.reshape(-1)

    def _snapshot_from_z(self, z_like) -> dict:
        if not torch.is_tensor(z_like):
            z_tensor = torch.tensor(z_like, device=self.device, dtype=self.dtype)
        else:
            z_tensor = z_like.to(device=self.device, dtype=self.dtype)

        x_proj = self.projector.project(z_tensor[: self.n_x]).view(self.n_mix, self.n_elements)

        if self.fix:
            atoms = self.atoms_base.copy()
            frac_np = np.array(atoms.get_scaled_positions(), dtype=float)
            cell_np = np.array(atoms.get_cell().array, dtype=float)
        else:
            uc = z_tensor[self.n_x : self.n_x + self.n_uc]
            frac = z_tensor[self.n_x + self.n_uc :].view(self.n_atoms, 3)
            frac_np = frac.detach().cpu().numpy()
            cell_np = ucparams_to_lattice(uc).detach().cpu().numpy()

        cart_np = frac_np @ cell_np
        weights_np = np.zeros((self.n_atoms, len(self.species_order)), dtype=float)
        occupancies_np = np.zeros(self.n_atoms, dtype=float)
        vacuum_np = np.zeros(self.n_atoms, dtype=float)

        if self.n_mix > 0:
            mix_rows = self.mix_indices.detach().cpu().numpy()
            x_proj_np = x_proj.detach().cpu().numpy()
            weights_np[np.ix_(mix_rows, self.mix_species_cols)] = x_proj_np
            occupancies_np[mix_rows] = np.sum(x_proj_np, axis=1)
            null_mask_np = self.mix_is_null.detach().cpu().numpy().astype(bool)
            if np.any(null_mask_np):
                vacuum_np[mix_rows[null_mask_np]] = np.clip(1.0 - occupancies_np[mix_rows[null_mask_np]], 0.0, 1.0)

        for atom_idx in self.fixed_indices.detach().cpu().numpy().tolist():
            symbol = self.base_symbols[atom_idx]
            weights_np[atom_idx, self.species_to_index[symbol]] = 1.0
            occupancies_np[atom_idx] = 1.0

        dominant_idx = np.argmax(weights_np, axis=1)
        dominant_species = [self.species_order[int(i)] for i in dominant_idx.tolist()]

        return {
            "frac_coords": frac_np.tolist(),
            "cart_positions_angstrom": cart_np.tolist(),
            "lattice_vectors_angstrom": cell_np.tolist(),
            "species_weights": weights_np.tolist(),
            "site_occupancies": occupancies_np.tolist(),
            "vacuum_probabilities": vacuum_np.tolist(),
            "dominant_species": dominant_species,
        }

    def fun_jac(self, z_flat: torch.Tensor):
        if not torch.is_tensor(z_flat):
            z_flat = torch.tensor(z_flat, device=self.device, dtype=self.dtype)
        else:
            z_flat = z_flat.to(device=self.device, dtype=self.dtype)

        if self.fix:
            x_raw = z_flat.detach().clone().requires_grad_(True)
            atoms = self.atoms_base.copy()
            frac = torch.tensor(atoms.get_scaled_positions(), device=self.device, dtype=self.dtype).detach().clone().requires_grad_(True)
            lattice = np.array(atoms.get_cell().array, dtype=float)
            uc = torch.tensor(cell_to_ucparams(lattice), device=self.device, dtype=self.dtype).detach().clone().requires_grad_(True)
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

        Xmix = self.projector.project(x_raw).view(self.n_mix, self.n_elements)

        if self.n_mix > 0:
            X_entropy = Xmix.detach().clone().requires_grad_(True)
            entropy_total, species_entropy_total, null_entropy_total = null_total_entropy(
                X_entropy,
                self.Tsallis_q,
                self.mix_is_null,
                eps=self.null_normalization_eps,
            )
            grad_entropy_flat = torch.autograd.grad(entropy_total, X_entropy)[0].reshape(-1)
        else:
            entropy_total = torch.zeros((), device=self.device, dtype=self.dtype)
            species_entropy_total = torch.zeros((), device=self.device, dtype=self.dtype)
            null_entropy_total = torch.zeros((), device=self.device, dtype=self.dtype)
            grad_entropy_flat = torch.zeros(self.n_x, device=self.device, dtype=self.dtype)

        entropy_cost = self.tau * entropy_total
        grad_entropy_cost = self.tau * grad_entropy_flat

        n_fairchem_elements = 100
        eat_fairchem = torch.zeros((self.n_atoms, n_fairchem_elements), device=self.device, dtype=Xmix.dtype)

        fixed_idx = torch.nonzero(self.site_fixed_Z > 0, as_tuple=False).view(-1)
        if fixed_idx.numel() > 0:
            Z_fixed = self.site_fixed_Z[fixed_idx].clamp(min=0, max=n_fairchem_elements - 1)
            eat_fairchem[fixed_idx, Z_fixed] = 1.0

        fairchem_element_indices = torch.tensor(
            [atomic_numbers[e] for e in self.selected_elements],
            dtype=torch.long,
            device=self.device,
        )
        if self.n_mix > 0:
            cols = fairchem_element_indices.unsqueeze(0).expand(self.n_mix, -1)
            rows = self.mix_indices.unsqueeze(1).expand(-1, self.n_elements)
            eat_fairchem[rows, cols] += Xmix

        cell = ucparams_to_lattice(uc)
        atoms = self.atoms_base.copy()
        if self.fix is False:
            atoms.set_scaled_positions(frac.detach().cpu().numpy())
            atoms.set_cell(cell.detach().cpu().numpy(), scale_atoms=True)

        E, F_cart, S = self.get_E_F_S(
            atoms,
            eat_weights_np=eat_fairchem.detach().cpu().numpy(),
            enable_eat_grad=True,
        )

        if "eat_grad" not in getattr(self.calc, "results", {}):
            raise RuntimeError(
                "FAIRChem did not return 'eat_grad'. Your fairchem build must support enable_eat_grad."
            )

        eat_grad_full = torch.tensor(self.calc.results["eat_grad"], device=self.device, dtype=Xmix.dtype)
        eat_grad_mix = eat_grad_full[self.mix_indices]
        cols = fairchem_element_indices.unsqueeze(0).expand(self.n_mix, -1)
        grad_Xmix = eat_grad_mix.gather(1, cols)
        grad_x_fairchem = grad_Xmix.reshape(-1)

        F_t = torch.tensor(F_cart, device=self.device, dtype=Xmix.dtype)
        grad_frac_fairchem = -(F_t @ cell.T).reshape(-1)

        volume = float(atoms.get_volume())
        S_t = torch.tensor(S, device=self.device, dtype=Xmix.dtype)
        HinvT = torch.inverse(cell).T
        dE_dH = volume * (HinvT @ S_t)

        J = jacobian_ucparams_to_lattice(uc)
        grad_uc_fairchem = torch.zeros(6, device=self.device, dtype=Xmix.dtype)
        for k in range(6):
            grad_uc_fairchem[k] = torch.sum(dE_dH * J[:, :, k])

        E_t = torch.tensor(E, device=self.device, dtype=Xmix.dtype)
        mix_counts = Xmix.sum(dim=0)
        total_atom_count = self.fixed_atom_count + float(torch.sum(mix_counts).detach().item())
        total_atom_count_t = torch.tensor(total_atom_count, device=self.device, dtype=Xmix.dtype)
        total_atom_count_safe = torch.clamp(total_atom_count_t, min=self.null_normalization_eps)

        ref_total = torch.dot(mix_counts.to(self.dtype), self.ref_energies_tensor.to(self.dtype)) + self.fixed_ref_energy_sum.to(self.dtype)
        numerator = E_t.to(self.dtype) - ref_total
        E_form = numerator / total_atom_count_safe

        ref_row = self.ref_energies_tensor.to(dtype=Xmix.dtype).view(1, -1).expand(self.n_mix, -1)
        grad_x_form = (grad_Xmix - ref_row) / total_atom_count_safe
        if self.n_null_mix > 0:
            grad_x_form = grad_x_form.clone()
            grad_x_form[self.mix_is_null] -= E_form.to(dtype=Xmix.dtype) / total_atom_count_safe
        grad_x_form = grad_x_form.reshape(-1)
        grad_uc_form = grad_uc_fairchem / total_atom_count_safe
        grad_frac_form = grad_frac_fairchem / total_atom_count_safe

        selected_counts, selected_total, selected_total_safe, composition_selected = self._selected_counts_and_composition(Xmix)
        E_hull, grad_x_hull = calc_hull(
            Xmix,
            self.hullcalc,
            null_mask=self.mix_is_null,
            fixed_counts=self.fixed_selected_counts,
            min_total=self.null_normalization_eps,
        )
        if E_hull is None or grad_x_hull is None:
            E_hull = float("nan")
            grad_x_hull = torch.zeros(self.n_x, device=self.device, dtype=Xmix.dtype)

        E_hull_t = torch.tensor(E_hull, device=self.device, dtype=Xmix.dtype)
        E_above_hull = E_form.to(dtype=Xmix.dtype) - E_hull_t

        penalty_coeff = self.composition_penalty_multiplier
        comp_delta = composition_selected - self.mu.to(device=self.device, dtype=Xmix.dtype)
        penalty = penalty_coeff * self.n_mix * 0.5 * torch.sum(comp_delta**2)
        grad_penalty_c = penalty_coeff * self.n_mix * comp_delta
        grad_penalty = self._lift_composition_gradient(grad_penalty_c, composition_selected, selected_total_safe)

        cost = entropy_cost + E_form.to(dtype=Xmix.dtype) + penalty
        grad_x = grad_entropy_cost + grad_x_form + grad_penalty
        grad_uc = grad_uc_form
        grad_f = grad_frac_form

        if self.relax:
            grad_x = grad_x * 0.0

        if self.fix is False:
            grad = torch.cat([grad_x, grad_uc, grad_f], dim=0)
        else:
            grad = grad_x

        if torch.isnan(grad).any() or torch.isinf(grad).any():
            raise RuntimeError("NaN/Inf detected in gradient.")

        with torch.no_grad():
            cell_volume = float(atoms.get_volume())
            gf = grad_f.view(self.n_atoms, 3)
            force_norm = float(torch.norm(gf, dim=1).mean().item())
            fmax = float(torch.norm(F_t, dim=1).max().item())

            discreteness_sq = torch.zeros((), device=self.device, dtype=Xmix.dtype)
            if self.n_mix > 0:
                normal_mask = ~self.mix_is_null
                if torch.any(normal_mask):
                    normal_probs = Xmix[normal_mask]
                    normal_target = torch.zeros_like(normal_probs)
                    normal_target.scatter_(1, torch.argmax(normal_probs, dim=1, keepdim=True), 1.0)
                    discreteness_sq = discreteness_sq + torch.sum((normal_probs - normal_target) ** 2)
                if torch.any(self.mix_is_null):
                    null_species = Xmix[self.mix_is_null]
                    null_empty = torch.clamp(1.0 - null_species.sum(dim=1, keepdim=True), min=0.0)
                    null_probs = torch.cat([null_species, null_empty], dim=1)
                    null_target = torch.zeros_like(null_probs)
                    null_target.scatter_(1, torch.argmax(null_probs, dim=1, keepdim=True), 1.0)
                    discreteness_sq = discreteness_sq + torch.sum((null_probs - null_target) ** 2)
            discreteness_dist = float(torch.sqrt(discreteness_sq).item())

            self.last_composition_vector = composition_selected.detach().cpu().numpy()
            self.last_metrics["energy_evaluations"] = int(self.energy_evaluations)
            self.last_metrics["cell_volume"] = cell_volume
            self.last_metrics["force_norm"] = force_norm
            self.last_metrics["fmax"] = fmax
            self.last_metrics["entropy_total"] = float(entropy_total.detach().item())
            self.last_metrics["species_entropy_total"] = float(species_entropy_total.detach().item())
            self.last_metrics["null_entropy_total"] = float(null_entropy_total.detach().item())
            self.last_metrics["entropy_per_atom"] = float((entropy_total / total_atom_count_safe).detach().item())
            self.last_metrics["cost"] = float(cost.detach().item())
            self.last_metrics["total_energy"] = float(E)
            self.last_metrics["formation_energy"] = float(E_form.detach().item())
            self.last_metrics["energy_above_hull"] = float(E_above_hull.detach().item())
            self.last_metrics["discreteness_dist"] = discreteness_dist
            self.last_metrics["mix_total_occupancy"] = float(torch.sum(mix_counts).detach().item())
            self.last_metrics["total_atom_count"] = float(total_atom_count_safe.detach().item())
            self.last_metrics["composition_penalty"] = float(penalty.detach().item())
            self.last_metrics["objective_no_entropy"] = float((E_form.to(dtype=Xmix.dtype) + penalty).detach().item())

        _ = grad_x_hull
        _ = selected_counts
        _ = selected_total

        return cost.detach(), grad.detach()

    def get_E_F_S(self, atoms, eat_weights_np=None, enable_eat_grad: bool = False):
        if eat_weights_np is not None:
            if "eat_weights" in atoms.arrays:
                atoms.arrays["eat_weights"][:] = np.asarray(eat_weights_np)
            else:
                atoms.new_array("eat_weights", np.asarray(eat_weights_np))
        if enable_eat_grad:
            atoms.info["enable_eat_grad"] = True
        else:
            atoms.info.pop("enable_eat_grad", None)

        atoms.calc = self.calc
        self.energy_evaluations += 1
        E = atoms.get_potential_energy()
        F = atoms.get_forces()
        S = atoms.get_stress(voigt=False)
        return E, F, S

    def pbfgs_callback(self, z, f, err, k, converged, *args):
        self.z_last = z.detach().clone() if torch.is_tensor(z) else z
        snapshot = self._snapshot_from_z(self.z_last)

        self.all_history["iter"].append(self.iter)
        for key, value in self.last_metrics.items():
            self.all_history[key].append(value if value is not None else np.nan)
        self.all_history["frac_coords"].append(snapshot["frac_coords"])
        self.all_history["cart_positions_angstrom"].append(snapshot["cart_positions_angstrom"])
        self.all_history["lattice_vectors_angstrom"].append(snapshot["lattice_vectors_angstrom"])
        self.all_history["species_weights"].append(snapshot["species_weights"])
        self.all_history["site_occupancies"].append(snapshot["site_occupancies"])
        self.all_history["vacuum_probabilities"].append(snapshot["vacuum_probabilities"])
        self.all_history["dominant_species"].append(snapshot["dominant_species"])
        self.all_history["trajectory"].append(
            {
                "iter": int(self.iter),
                "energy_evaluations": int(self.energy_evaluations),
                "frac_coords": snapshot["frac_coords"],
                "cart_positions_angstrom": snapshot["cart_positions_angstrom"],
                "lattice_vectors_angstrom": snapshot["lattice_vectors_angstrom"],
                "species_weights": snapshot["species_weights"],
                "site_occupancies": snapshot["site_occupancies"],
                "vacuum_probabilities": snapshot["vacuum_probabilities"],
                "dominant_species": snapshot["dominant_species"],
                "cost": self.last_metrics.get("cost", np.nan),
                "total_energy": self.last_metrics.get("total_energy", np.nan),
                "formation_energy": self.last_metrics.get("formation_energy", np.nan),
                "energy_above_hull": self.last_metrics.get("energy_above_hull", np.nan),
                "entropy_per_atom": self.last_metrics.get("entropy_per_atom", np.nan),
                "cell_volume": self.last_metrics.get("cell_volume", np.nan),
                "mix_total_occupancy": self.last_metrics.get("mix_total_occupancy", np.nan),
                "total_atom_count": self.last_metrics.get("total_atom_count", np.nan),
            }
        )

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
        self.history["energy_evaluations"].append(int(self.energy_evaluations))
        self.history["entropy"].append(ent if ent is not None else np.nan)
        self.history["cost"].append(cost if cost is not None else np.nan)
        self.history["total_energy"].append(total_energy if total_energy is not None else np.nan)
        self.history["formation_energy"].append(formation_energy if formation_energy is not None else np.nan)
        self.history["energy_above_hull"].append(energy_above_hull if energy_above_hull is not None else np.nan)
        self.history["cell_volume"].append(cell_volume if cell_volume is not None else np.nan)
        self.history["force_norm"].append(force_norm if force_norm is not None else np.nan)
        self.history["fmax"].append(fmax if fmax is not None else np.nan)
        self.history["discreteness_dist"].append(
            self.last_metrics["discreteness_dist"] if self.last_metrics["discreteness_dist"] is not None else np.nan
        )
        self.history["mix_total_occupancy"].append(
            self.last_metrics["mix_total_occupancy"] if self.last_metrics["mix_total_occupancy"] is not None else np.nan
        )
        self.history["total_atom_count"].append(
            self.last_metrics["total_atom_count"] if self.last_metrics["total_atom_count"] is not None else np.nan
        )

        if cost is not None:
            if self.last_cost_for_plateau is not None:
                if abs(cost - self.last_cost_for_plateau) <= self.tol_cost:
                    self.plateau_count += 1
                else:
                    self.plateau_count = 0
            self.last_cost_for_plateau = cost

            if self.plateau_count >= self.nIterations_tol:
                self.log.info(
                    f"[PBFGS iter {int(k):03d}] Plateau convergence reached: "
                    f"|Δcost| <= {self.tol_cost} for {self.plateau_count} steps. "
                    f"Stopping PBFGS early for this anneal."
                )
                raise PBFGSPlateauReached()

        self.iter += 1

        if ent is None or cost is None:
            self.log.info(f"[PBFGS iter {k:03d}] f={f_val:.6f}, PG-error={err:.3e}")
            return

        comp1 = float(self.last_composition_vector[0]) if self.last_composition_vector is not None else float("nan")
        self.log.info(
            f"[PBFGS iter {k:03d}] "
            f"cost={cost:.6f}, E_above_hull={energy_above_hull:.6f} eV/atom, "
            f"V={cell_volume:.2f} Å^3, fmax={fmax:.3e}, "
            f"ent/atom={ent:.3e}, PG-error={err:.3e}, comp1={comp1:.3f}"
        )
