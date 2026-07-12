import math
from dataclasses import dataclass
from typing import Dict, List

import numpy as np
import torch
from ase.data import atomic_numbers

from .energy_above_hull import calc_hull
from .function_helpers import null_total_entropy
from .null_geometry import perturb_fractional_positions, repair_null_overlaps
from .lattice_helpers import cell_to_ucparams, ucparams_to_lattice
from .minimizer import minimizer
from .objfun import PBFGSPlateauReached, objfun
from .projections import projector


@dataclass
class OptimizationResult:
    history: Dict[str, List[float]]
    all_history: Dict[str, object]
    anneal_boundaries: List[int]
    final_energy_total_ev: float
    final_formation_energy_ev_per_atom: float
    final_energy_above_hull_ev_per_atom: float
    selected_elements: List[str]
    X_final: np.ndarray
    site_occupancies: np.ndarray
    composition_selected: np.ndarray
    frac_coords: np.ndarray
    lattice_vectors: np.ndarray
    symbols: List[str]
    db_label: str
    energy_evaluations: int
    cum_global_iter: int
    aborted: bool = False
    abort_reason: str = ""


def _project_capped_simplex(v: torch.Tensor, target: float, iters: int = 80) -> torch.Tensor:
    x = v.reshape(-1)
    n = x.numel()
    target_eff = float(target)
    if target_eff < 0.0 or target_eff > float(n):
        raise ValueError(f"Infeasible target occupancy sum: {target_eff} for {n} sites.")

    if target_eff == 0.0:
        return torch.zeros_like(x)
    if target_eff == float(n):
        return torch.ones_like(x)

    tau_lo = (x - 1.0).min().item()
    tau_hi = x.max().item()
    for _ in range(iters):
        tau = 0.5 * (tau_lo + tau_hi)
        y = torch.clamp(x - tau, 0.0, 1.0)
        if y.sum().item() > target_eff:
            tau_lo = tau
        else:
            tau_hi = tau
    return torch.clamp(x - 0.5 * (tau_lo + tau_hi), 0.0, 1.0)


def _project_ucparams_tensor(
    ucparams: torch.Tensor,
    cell_length_min: float,
    cell_length_max: float,
    cell_angle_min_deg: float,
    cell_angle_max_deg: float,
) -> torch.Tensor:
    out = ucparams.clone()
    out[:3] = torch.clamp(out[:3], float(cell_length_min), float(cell_length_max))
    angle_min = math.radians(float(cell_angle_min_deg))
    angle_max = math.radians(float(cell_angle_max_deg))
    out[3:] = torch.clamp(out[3:], angle_min, angle_max)
    return out


class syntropizer:
    def __init__(
        self,
        selected_elements,
        atoms_base,
        calc,
        log,
        mu,
        fix,
        alpha,
        DB_JSON_PATH,
        REFERENCE_JSON_PATH,
        n_anneals=1000,
        length_rel_delta=1.0,
        angle_deg_margin=40,
        angle_min_deg=40,
        angle_max_deg=140,
        tol_cost=1e-3,
        nIterations_tol=5,
        tau_init=1.0,
        tau_mult=0.9,
        first_stage_maxit=None,
        stage_maxit=100,
        pbfgs_tol=1e-3,
        Tsallis_q=2.0,
        pbfgs_maxit_relax=500,
        pbfgs_tol_relax=5e-3,
        null_min_distance=1.5,
        null_normalization_eps=1e-8,
        null_threshold=1e-3,
        initial_occupancy_sum=None,
        cell_init_mode="input_cell",
        random_cell_length_min=2.0,
        random_cell_length_max=5.0,
        random_cell_angle_min_deg=60.0,
        random_cell_angle_max_deg=120.0,
        cell_bound_length_min=0.5,
        cell_bound_length_max=100.0,
        cell_bound_angle_min_deg=20.0,
        cell_bound_angle_max_deg=160.0,
        composition_penalty_multiplier=1.0,
        escape_enabled=True,
        escape_trials=3,
        escape_pos_sigma=0.3,
        escape_empty_prob_sigma=0.5,
        escape_species_prob_sigma=0.5,
        escape_cell_length_sigma=0.0,
        escape_cell_angle_sigma_deg=0.0,
        escape_trial_maxit=50,
        escape_stagnation_window=3,
        escape_min_shock_interval=1,
        escape_max_attempts_per_stage=3,
        escape_accept_tol=1.0e-3,
        restart_parent_result=None,
        restart_memory_strength=0.75,
        initial_lattice_vectors=None,
        abort_energy_above_hull=None,
        device=None,
    ):
        self.device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")

        self.selected_elements = list(selected_elements)
        self.atoms_base = atoms_base
        self.calc = calc
        self.log = log
        self.mu = mu.to(self.device) if torch.is_tensor(mu) else torch.tensor(mu, device=self.device, dtype=torch.float32)
        self.fix = fix
        self.alpha = alpha

        self.DB_JSON_PATH = DB_JSON_PATH
        self.REFERENCE_JSON_PATH = REFERENCE_JSON_PATH

        self.site_is_mix = torch.tensor(atoms_base.arrays["site_is_mix"], device=self.device).bool()
        self.site_is_null = torch.tensor(
            atoms_base.arrays.get("site_is_null", np.zeros(len(atoms_base), dtype=np.int32)),
            device=self.device,
        ).bool()
        self.site_fixed_Z = torch.tensor(atoms_base.arrays["site_fixed_Z"], device=self.device).long()

        self.mix_indices = torch.nonzero(self.site_is_mix, as_tuple=False).view(-1)
        self.fixed_indices = torch.nonzero(~self.site_is_mix, as_tuple=False).view(-1)
        self.mix_is_null = self.site_is_null[self.mix_indices]
        self.base_symbols = list(atoms_base.get_chemical_symbols())
        self.species_order = list(self.selected_elements)
        for atom_idx in self.fixed_indices.detach().cpu().numpy().tolist():
            symbol = self.base_symbols[atom_idx]
            if symbol not in self.species_order:
                self.species_order.append(symbol)
        self.species_to_index = {symbol: i for i, symbol in enumerate(self.species_order)}
        self.mix_species_cols = np.array([self.species_to_index[symbol] for symbol in self.selected_elements], dtype=int)

        self.n_atoms = len(atoms_base)
        self.n_mix = int(self.mix_indices.numel())
        self.n_elements = len(self.selected_elements)
        self.n_x = self.n_mix * self.n_elements
        self.n_uc = 6

        self.n_anneals = int(n_anneals)
        self.length_rel_delta = float(length_rel_delta)
        self.tau_init = float(tau_init)
        self.tau_mult = float(tau_mult)
        self.first_stage_maxit = int(stage_maxit if first_stage_maxit is None else first_stage_maxit)
        self.stage_maxit = int(stage_maxit)
        self.pbfgs_tol = float(pbfgs_tol)
        self.Tsallis_q = float(Tsallis_q)
        self.null_min_distance = float(null_min_distance)
        self.null_normalization_eps = float(null_normalization_eps)
        self.null_threshold = float(null_threshold)
        self.initial_occupancy_sum = initial_occupancy_sum
        self.cell_init_mode = str(cell_init_mode)
        self.random_cell_length_min = float(random_cell_length_min)
        self.random_cell_length_max = float(random_cell_length_max)
        self.random_cell_angle_min_deg = float(random_cell_angle_min_deg)
        self.random_cell_angle_max_deg = float(random_cell_angle_max_deg)
        self.cell_bound_length_min = float(cell_bound_length_min)
        self.cell_bound_length_max = float(cell_bound_length_max)
        self.cell_bound_angle_min_deg = float(cell_bound_angle_min_deg)
        self.cell_bound_angle_max_deg = float(cell_bound_angle_max_deg)
        self.composition_penalty_multiplier = float(composition_penalty_multiplier)

        self.escape_enabled = bool(escape_enabled)
        self.escape_trials = int(escape_trials)
        self.escape_pos_sigma = float(escape_pos_sigma)
        self.escape_empty_prob_sigma = float(escape_empty_prob_sigma)
        self.escape_species_prob_sigma = float(escape_species_prob_sigma)
        self.escape_cell_length_sigma = float(escape_cell_length_sigma)
        self.escape_cell_angle_sigma_deg = float(escape_cell_angle_sigma_deg)
        self.escape_trial_maxit = int(escape_trial_maxit)
        self.escape_stagnation_window = int(escape_stagnation_window)
        if self.escape_stagnation_window < 2:
            raise ValueError("escape_stagnation_window must be >= 2")
        self.escape_min_shock_interval = int(escape_min_shock_interval)
        if self.escape_min_shock_interval < 1:
            raise ValueError("escape_min_shock_interval must be >= 1")
        self.escape_max_attempts_per_stage = int(escape_max_attempts_per_stage)
        if self.escape_max_attempts_per_stage < 0:
            raise ValueError("escape_max_attempts_per_stage must be >= 0")
        self.escape_accept_tol = float(escape_accept_tol)
        if self.escape_accept_tol < 0.0:
            raise ValueError("escape_accept_tol must be >= 0")
        self.restart_parent_result = restart_parent_result
        self.restart_memory_strength = float(restart_memory_strength)
        self.initial_lattice_vectors = initial_lattice_vectors
        self.abort_energy_above_hull = None if abort_energy_above_hull is None else float(abort_energy_above_hull)

        self.tol_cost = float(tol_cost)
        self.nIterations_tol = int(nIterations_tol)

        self.pbfgs_maxit_relax = int(pbfgs_maxit_relax)
        self.pbfgs_tol_relax = float(pbfgs_tol_relax)
        if self.first_stage_maxit <= 0:
            raise ValueError("first_stage_maxit must be > 0")
        if self.stage_maxit <= 0:
            raise ValueError("stage_maxit must be > 0")
        if self.pbfgs_tol <= 0.0:
            raise ValueError("pbfgs_tol must be > 0")
        if self.null_threshold < 0.0:
            raise ValueError("null_threshold must be >= 0")

        frac0_np = atoms_base.get_scaled_positions()
        self.f0 = torch.tensor(frac0_np, dtype=torch.float32, device=self.device)

        cell0_np = atoms_base.cell.array.copy()
        uc0_np = cell_to_ucparams(cell0_np)
        self.uc0 = torch.tensor(uc0_np, dtype=torch.float32, device=self.device)
        self.initial_uc_override = None
        if self.initial_lattice_vectors is not None:
            initial_cell_np = np.asarray(self.initial_lattice_vectors, dtype=float)
            if initial_cell_np.shape != (3, 3):
                raise ValueError("initial_lattice_vectors must be a 3x3 lattice matrix")
            initial_uc_np = cell_to_ucparams(initial_cell_np)
            self.initial_uc_override = torch.tensor(initial_uc_np, dtype=torch.float32, device=self.device)

        self.uc_lower = torch.tensor(
            [
                self.cell_bound_length_min,
                self.cell_bound_length_min,
                self.cell_bound_length_min,
                math.radians(self.cell_bound_angle_min_deg),
                math.radians(self.cell_bound_angle_min_deg),
                math.radians(self.cell_bound_angle_min_deg),
            ],
            dtype=torch.float32,
            device=self.device,
        )
        self.uc_upper = torch.tensor(
            [
                self.cell_bound_length_max,
                self.cell_bound_length_max,
                self.cell_bound_length_max,
                math.radians(self.cell_bound_angle_max_deg),
                math.radians(self.cell_bound_angle_max_deg),
                math.radians(self.cell_bound_angle_max_deg),
            ],
            dtype=torch.float32,
            device=self.device,
        )

        def _assert_uc_within_bounds(uc, label):
            below = bool(torch.any(uc < self.uc_lower - 1e-6).item())
            above = bool(torch.any(uc > self.uc_upper + 1e-6).item())
            if below or above:
                raise ValueError(
                    f"{label} cell is outside absolute cell_bounds. "
                    "Increase cell_bounds.length_min/max or cell_bounds.angle_min_deg/max_deg."
                )

        if self.cell_init_mode != "random":
            _assert_uc_within_bounds(self.uc0, "input")
        if self.initial_uc_override is not None:
            _assert_uc_within_bounds(self.initial_uc_override, "population.initial_lattice_vectors")

        self.anneal_boundaries: List[int] = []
        self.history: Dict[str, List[float]] = {}
        self.global_iter = 0
        self.energy_evaluations = 0
        self.z = None

        self.projector = projector(
            self.n_mix,
            self.n_elements,
            self.fix,
            self.mu,
            self.uc_lower,
            self.uc_upper,
            mix_is_null=self.mix_is_null,
        )

        self.all_history: Dict[str, List] = {}
        self.all_history["schema_version"] = 3
        self.all_history["escape_rounds"] = []
        self.all_history["stages"] = []
        self.all_history["initial_state"] = None
        self.all_history["final_relaxation"] = None

    def _escape_reseed_candidates(self, min_count: int = 0):
        n_candidates = max(int(min_count), int(min_count) * 16, 64)
        return np.random.random((n_candidates, 3))

    def _periodic_pairwise_dist2_frac(self, frac_a: np.ndarray, frac_b: np.ndarray, cell_np: np.ndarray) -> np.ndarray:
        frac_a = np.asarray(frac_a, dtype=float)
        frac_b = np.asarray(frac_b, dtype=float)
        if frac_a.size == 0 or frac_b.size == 0:
            return np.empty((len(frac_a), len(frac_b)), dtype=float)
        ds = frac_a[:, None, :] - frac_b[None, :, :]
        ds = ds - np.round(ds)
        dr = ds @ np.asarray(cell_np, dtype=float)
        return np.sum(dr**2, axis=-1)

    def _make_random_initial_uc(self):
        uc_np = np.empty(6, dtype=float)
        uc_np[:3] = np.random.uniform(self.random_cell_length_min, self.random_cell_length_max, size=3)
        uc_np[3:] = np.radians(
            np.random.uniform(self.random_cell_angle_min_deg, self.random_cell_angle_max_deg, size=3)
        )
        uc_t = torch.tensor(uc_np, device=self.device, dtype=torch.float32)
        return _project_ucparams_tensor(
            uc_t,
            self.random_cell_length_min,
            self.random_cell_length_max,
            self.random_cell_angle_min_deg,
            self.random_cell_angle_max_deg,
        )

    def _round_projected_weights(self, Xmix):
        Xout = torch.zeros_like(Xmix)
        if self.n_mix == 0:
            return Xout

        max_idx = torch.argmax(Xmix, dim=1)
        null_mask = self.mix_is_null

        if torch.any(~null_mask):
            rows = torch.nonzero(~null_mask, as_tuple=False).view(-1)
            Xout[rows, max_idx[rows]] = 1.0

        if torch.any(null_mask):
            rows = torch.nonzero(null_mask, as_tuple=False).view(-1)
            empty = torch.clamp(1.0 - Xmix[rows].sum(dim=1), min=0.0)
            max_species_prob = Xmix[rows, max_idx[rows]]
            keep = max_species_prob >= empty
            if torch.any(keep):
                keep_rows = rows[keep]
                Xout[keep_rows, max_idx[keep_rows]] = 1.0

        return Xout

    def _initialize_species_weights(self):
        dist = torch.distributions.Dirichlet(
            self.alpha * torch.ones(self.n_elements, dtype=torch.float32, device=self.device)
        )
        X0 = dist.sample((self.n_mix,))

        if self.n_mix > 0 and torch.any(self.mix_is_null):
            null_rows = torch.nonzero(self.mix_is_null, as_tuple=False).view(-1)
            n_null = int(null_rows.numel())
            target_sum = float(n_null if self.initial_occupancy_sum is None else self.initial_occupancy_sum)
            target_sum = max(0.0, min(float(n_null), target_sum))
            occ0 = _project_capped_simplex(
                torch.full((n_null,), target_sum / float(max(1, n_null)), device=self.device, dtype=torch.float32),
                target_sum,
            )
            X0 = X0.clone()
            X0[null_rows] = X0[null_rows] * occ0.unsqueeze(1)

        return X0

    def _zero_close_null_occupancies(
        self,
        Xmix: torch.Tensor,
        frac: torch.Tensor,
        uc: torch.Tensor,
        *,
        label: str,
    ) -> tuple[torch.Tensor, list[int]]:
        if self.n_mix == 0 or (not torch.any(self.mix_is_null)) or self.null_min_distance <= 0.0:
            return Xmix, []

        x_rows = Xmix.detach().clone()
        frac_np = frac.detach().cpu().numpy()
        cell_np = ucparams_to_lattice(uc).detach().cpu().numpy()
        cutoff2 = float(self.null_min_distance) ** 2

        mix_atom_indices = self.mix_indices.detach().cpu().numpy()
        fixed_atom_indices = self.fixed_indices.detach().cpu().numpy()
        null_local_indices = np.where(self.mix_is_null.detach().cpu().numpy().astype(bool))[0]
        null_atom_indices = mix_atom_indices[null_local_indices]
        null_local_by_atom = {int(atom_idx): int(local_idx) for local_idx, atom_idx in zip(null_local_indices, null_atom_indices)}
        normal_mix_local = np.where((~self.mix_is_null).detach().cpu().numpy().astype(bool))[0]
        normal_atom_indices = mix_atom_indices[normal_mix_local]
        non_null_atom_indices = np.concatenate([fixed_atom_indices, normal_atom_indices]).astype(int)

        zeroed: set[int] = set()
        while True:
            occ_np = x_rows.sum(dim=1).detach().cpu().numpy()
            active_null_local = null_local_indices[occ_np[null_local_indices] > self.null_threshold]
            active_null_atoms = mix_atom_indices[active_null_local]
            considered_atoms = np.concatenate([non_null_atom_indices, active_null_atoms]).astype(int)
            if considered_atoms.size < 2:
                break

            considered_frac = frac_np[considered_atoms]
            dist2 = self._periodic_pairwise_dist2_frac(considered_frac, considered_frac, cell_np)
            close_pairs = np.column_stack(np.triu_indices(len(considered_atoms), k=1))[dist2[np.triu_indices(len(considered_atoms), k=1)] < cutoff2]
            if close_pairs.size == 0:
                break

            to_zero: set[int] = set()
            for pair in close_pairs:
                atom_i = int(considered_atoms[int(pair[0])])
                atom_j = int(considered_atoms[int(pair[1])])
                local_i = null_local_by_atom.get(atom_i)
                local_j = null_local_by_atom.get(atom_j)
                if local_i is None and local_j is None:
                    continue
                if local_i is not None and local_j is None:
                    to_zero.add(local_i)
                    continue
                if local_i is None and local_j is not None:
                    to_zero.add(local_j)
                    continue

                occ_i = float(occ_np[local_i])
                occ_j = float(occ_np[local_j])
                if occ_i < occ_j:
                    to_zero.add(local_i)
                elif occ_j < occ_i:
                    to_zero.add(local_j)
                else:
                    to_zero.add(max(local_i, local_j))

            to_zero -= zeroed
            if not to_zero:
                break
            zero_idx = torch.as_tensor(sorted(to_zero), device=self.device, dtype=torch.long)
            x_rows[zero_idx] = 0.0
            zeroed.update(to_zero)

        zeroed_list = sorted(int(i) for i in zeroed)
        if zeroed_list:
            self.log.info(
                f"{label}: kept all {int(self.mix_is_null.sum().item())} null sites; "
                f"set {len(zeroed_list)} close-overlap null occupancies to zero."
            )
        return x_rows, zeroed_list

    def _initialize_from_parent_result(self):
        parent = self.restart_parent_result
        if parent is None:
            return None

        X_parent = torch.tensor(parent.X_final, dtype=torch.float32, device=self.device)
        if X_parent.shape != (self.n_mix, self.n_elements):
            raise ValueError(
                "restart_parent_result.X_final has shape "
                f"{tuple(X_parent.shape)}, expected {(self.n_mix, self.n_elements)}"
            )

        mem = max(0.0, min(1.0, float(self.restart_memory_strength)))
        dist = torch.distributions.Dirichlet(
            self.alpha * torch.ones(self.n_elements, dtype=torch.float32, device=self.device)
        )
        species_random = dist.sample((self.n_mix,))

        parent_species = torch.argmax(X_parent, dim=1)
        species_onehot = torch.zeros_like(X_parent)
        if self.n_mix > 0:
            species_onehot.scatter_(1, parent_species.view(-1, 1), 1.0)

        parent_occ = torch.clamp(X_parent.sum(dim=1), min=0.0, max=1.0)
        row_memory = torch.full((self.n_mix, 1), mem, dtype=torch.float32, device=self.device)
        if self.n_mix > 0 and torch.any(self.mix_is_null):
            null_rows = torch.nonzero(self.mix_is_null, as_tuple=False).view(-1)
            null_parent_alive = parent_occ[null_rows] > self.null_threshold
            row_memory[null_rows] = torch.where(
                null_parent_alive.view(-1, 1),
                row_memory[null_rows],
                torch.zeros_like(row_memory[null_rows]),
            )

        species_pref = (1.0 - row_memory) * species_random + row_memory * species_onehot
        if self.escape_species_prob_sigma > 0.0 and self.n_mix > 0:
            species_pref = species_pref + self.escape_species_prob_sigma * torch.randn_like(species_pref)
        species_dist = self.projector._project_simplex_rows(species_pref) if self.n_mix > 0 else species_pref
        X0 = species_dist.clone()

        if self.n_mix > 0 and torch.any(self.mix_is_null):
            null_rows = torch.nonzero(self.mix_is_null, as_tuple=False).view(-1)
            n_null = int(null_rows.numel())
            target_sum = float(n_null if self.initial_occupancy_sum is None else self.initial_occupancy_sum)
            target_sum = max(0.0, min(float(n_null), target_sum))
            base_occ = torch.full(
                (n_null,),
                target_sum / float(max(1, n_null)),
                dtype=torch.float32,
                device=self.device,
            )
            parent_alive = (parent_occ[null_rows] > self.null_threshold).to(dtype=torch.float32)
            occ_scores = (1.0 - mem) * base_occ + mem * parent_alive
            if self.escape_empty_prob_sigma > 0.0:
                occ_scores = occ_scores + self.escape_empty_prob_sigma * torch.randn_like(occ_scores)
            occ0 = _project_capped_simplex(occ_scores, target_sum)
            X0[null_rows] = species_dist[null_rows] * occ0.view(-1, 1)

        if self.fix:
            return X0.view(-1)

        parent_frac = torch.tensor(parent.frac_coords, dtype=torch.float32, device=self.device)
        if parent_frac.shape != (self.n_atoms, 3):
            raise ValueError(
                "restart_parent_result.frac_coords has shape "
                f"{tuple(parent_frac.shape)}, expected {(self.n_atoms, 3)}"
            )

        parent_uc = torch.tensor(
            cell_to_ucparams(np.asarray(parent.lattice_vectors, dtype=float)),
            dtype=torch.float32,
            device=self.device,
        )
        uc_init = parent_uc.clone()
        if self.escape_cell_length_sigma > 0.0:
            uc_init[:3] = uc_init[:3] + self.escape_cell_length_sigma * torch.randn(3, device=self.device)
        if self.escape_cell_angle_sigma_deg > 0.0:
            uc_init[3:] = uc_init[3:] + math.radians(self.escape_cell_angle_sigma_deg) * torch.randn(3, device=self.device)
        uc_init = torch.max(torch.min(uc_init, self.uc_upper), self.uc_lower)

        frac_init = self.f0.clone()
        keep_atom_mask = torch.zeros(self.n_atoms, dtype=torch.bool, device=self.device)
        if self.fixed_indices.numel() > 0:
            keep_atom_mask[self.fixed_indices] = True
        if self.n_mix > 0:
            normal_mix_rows = torch.nonzero(~self.mix_is_null, as_tuple=False).view(-1)
            if normal_mix_rows.numel() > 0:
                keep_atom_mask[self.mix_indices[normal_mix_rows]] = True
            alive_mix_rows = torch.nonzero(
                self.mix_is_null & (parent_occ > self.null_threshold),
                as_tuple=False,
            ).view(-1)
            if alive_mix_rows.numel() > 0:
                keep_atom_mask[self.mix_indices[alive_mix_rows]] = True
        frac_init[keep_atom_mask] = parent_frac[keep_atom_mask]

        cell_np = ucparams_to_lattice(uc_init).detach().cpu().numpy()
        frac_np = frac_init.detach().cpu().numpy()
        if self.escape_pos_sigma > 0.0:
            frac_np = perturb_fractional_positions(frac_np, cell_np, self.escape_pos_sigma)

        if self.n_mix > 0 and torch.any(self.mix_is_null) and self.null_min_distance > 0.0:
            x_rows = X0.clone()
            occ_np = x_rows.sum(dim=1).detach().cpu().numpy()
            mix_atom_indices = self.mix_indices.detach().cpu().numpy()
            null_local_indices = np.where(self.mix_is_null.detach().cpu().numpy().astype(bool))[0]
            null_atom_indices = mix_atom_indices[null_local_indices]
            priority = occ_np[null_local_indices] if occ_np.size > 0 else np.zeros(len(null_local_indices), dtype=float)
            repair = repair_null_overlaps(
                frac_np,
                cell_np,
                null_indices=null_atom_indices,
                candidate_frac=self._escape_reseed_candidates(min_count=int(len(null_atom_indices))),
                min_distance=self.null_min_distance,
                priority_scores=priority,
            )
            frac_np = repair.frac_all
            if repair.deactivated_null_local_indices:
                deactivate_mix_local = null_local_indices[np.asarray(repair.deactivated_null_local_indices, dtype=int)]
                x_rows = x_rows.clone()
                x_rows[torch.as_tensor(deactivate_mix_local, device=self.device, dtype=torch.long)] = 0.0
                X0 = x_rows

        frac_t = torch.tensor(frac_np, dtype=torch.float32, device=self.device)
        return torch.cat([X0.view(-1), uc_init, frac_t.reshape(-1)], dim=0)

    def _init_z(self):
        z_restart = self._initialize_from_parent_result()
        if z_restart is not None:
            self.z = self.projector.project(z_restart)
            self.all_history["restart_parent"] = {
                "memory_strength": float(self.restart_memory_strength),
                "formation_energy_ev_per_atom": float(self.restart_parent_result.final_formation_energy_ev_per_atom),
                "energy_above_hull_ev_per_atom": float(self.restart_parent_result.final_energy_above_hull_ev_per_atom),
            }
            return

        X0 = self._initialize_species_weights()

        if self.cell_init_mode == "random":
            uc_init = self._make_random_initial_uc()
        else:
            uc_init = self.uc0.clone()
        if self.initial_uc_override is not None:
            uc_init = self.initial_uc_override.clone()

        X0, zeroed_close_nulls = self._zero_close_null_occupancies(
            X0,
            self.f0,
            uc_init,
            label="initial null placement",
        )
        self.all_history["initial_zeroed_close_nulls"] = [int(i) for i in zeroed_close_nulls]
        x0 = X0.view(-1)

        if self.fix:
            self.z = x0
            return

        self.z = torch.cat([x0, uc_init, self.f0.reshape(-1)], dim=0)

    def _state_snapshot_from_z(
        self,
        z_state: torch.Tensor,
        *,
        stage_index: int,
        tau: float,
        run_kind: str,
        attempt_index: int,
        iter_idx: int,
        note: str | None = None,
        reseeded_nulls: int | None = None,
        metrics: dict | None = None,
    ) -> dict:
        Xmix, uc, frac = self._unpack_state(self.projector.project(z_state))
        frac_np = frac.detach().cpu().numpy()
        cell_np = ucparams_to_lattice(uc).detach().cpu().numpy()
        cart_np = frac_np @ cell_np

        weights_np = np.zeros((self.n_atoms, len(self.species_order)), dtype=float)
        occupancies_np = np.zeros(self.n_atoms, dtype=float)
        vacuum_np = np.zeros(self.n_atoms, dtype=float)

        if self.n_mix > 0:
            mix_rows = self.mix_indices.detach().cpu().numpy()
            x_proj_np = Xmix.detach().cpu().numpy()
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

        snapshot = {
            "stage": int(stage_index + 1),
            "tau": float(tau),
            "run_kind": str(run_kind),
            "attempt_index": int(attempt_index),
            "iter": int(iter_idx),
            "frac_coords": frac_np.tolist(),
            "cart_positions_angstrom": cart_np.tolist(),
            "lattice_vectors_angstrom": cell_np.tolist(),
            "species_weights": weights_np.tolist(),
            "site_occupancies": occupancies_np.tolist(),
            "vacuum_probabilities": vacuum_np.tolist(),
            "dominant_species": dominant_species,
        }
        if metrics is not None:
            for key, value in metrics.items():
                snapshot[key] = value
        if note is not None:
            snapshot["note"] = str(note)
        if reseeded_nulls is not None:
            snapshot["reseeded_nulls"] = int(reseeded_nulls)
        return snapshot

    def _make_parallel_run_record(
        self,
        *,
        z_start: torch.Tensor,
        z_out: torch.Tensor,
        objfun_instance,
        tau: float,
        stage_index: int,
        run_kind: str,
        attempt_index: int,
        f_true_history,
        formation_energy_history,
        cost_history,
        failed: bool = False,
        failure_reason: str | None = None,
    ) -> dict:
        def _safe_metric(name: str) -> float:
            value = objfun_instance.last_metrics.get(name, np.nan)
            if value is None:
                return float("nan")
            try:
                return float(value)
            except (TypeError, ValueError):
                return float("nan")

        trajectory = []
        raw_traj = objfun_instance.all_history.get("trajectory", [])
        for idx, snapshot in enumerate(raw_traj):
            snap = dict(snapshot)
            snap["stage"] = int(stage_index + 1)
            snap["tau"] = float(tau)
            snap["run_kind"] = str(run_kind)
            snap["attempt_index"] = int(attempt_index)
            snap["iter"] = int(snap.get("iter", idx + 1))
            trajectory.append(snap)

        final_metrics = {
            "cost": _safe_metric("cost"),
            "total_energy": _safe_metric("total_energy"),
            "formation_energy": _safe_metric("formation_energy"),
            "energy_above_hull": _safe_metric("energy_above_hull"),
            "entropy_per_atom": _safe_metric("entropy_per_atom"),
            "cell_volume": _safe_metric("cell_volume"),
            "mix_total_occupancy": _safe_metric("mix_total_occupancy"),
            "total_atom_count": _safe_metric("total_atom_count"),
            "composition_penalty": _safe_metric("composition_penalty"),
            "objective_no_entropy": _safe_metric("objective_no_entropy"),
        }

        initial_state = self._state_snapshot_from_z(
            z_start,
            stage_index=stage_index,
            tau=tau,
            run_kind=run_kind,
            attempt_index=attempt_index,
            iter_idx=0,
            note="pre_pbfgs_state",
        )
        if trajectory:
            final_state = dict(trajectory[-1])
            final_state["note"] = "post_pbfgs_final_state"
        else:
            final_state = self._state_snapshot_from_z(
                z_out,
                stage_index=stage_index,
                tau=tau,
                run_kind=run_kind,
                attempt_index=attempt_index,
                iter_idx=0,
                note="post_pbfgs_final_state",
                metrics=final_metrics,
            )

        stagnated, stagnation_change = self._objective_no_entropy_stagnated(f_true_history)
        energy_evaluations = int(getattr(objfun_instance, "energy_evaluations", 0))

        return {
            "run_kind": str(run_kind),
            "attempt_index": int(attempt_index),
            "tau": float(tau),
            "energy_evaluations": energy_evaluations,
            "failed": bool(failed),
            "failure_reason": None if failure_reason is None else str(failure_reason),
            "initial_state": initial_state,
            "trajectory": trajectory,
            "cost_history": list(cost_history) if cost_history is not None else [],
            "cost_final": (
                float(cost_history[-1])
                if cost_history
                else float(final_metrics["cost"])
            ),
            "f_true_history": list(f_true_history),
            "f_true_final": (
                float(f_true_history[-1])
                if f_true_history
                else float(final_metrics["objective_no_entropy"])
            ),
            "formation_energy_history": list(formation_energy_history),
            "formation_energy_final": (
                float(formation_energy_history[-1])
                if formation_energy_history
                else float(final_metrics["formation_energy"])
            ),
            "entropy_per_atom_final": float(final_metrics["entropy_per_atom"]),
            "composition_penalty_final": float(final_metrics["composition_penalty"]),
            "objective_no_entropy_final": float(final_metrics["objective_no_entropy"]),
            "stagnated": bool(stagnated),
            "stagnation_change": float(stagnation_change),
            "stagnation_window": int(self.escape_stagnation_window),
            "stagnation_tol": float(self.pbfgs_tol),
            "final_state": final_state,
        }

    def run(self) -> OptimizationResult:
        if self.z is None:
            self._init_z()
        self.all_history["initial_state"] = self._state_snapshot_from_z(
            self.z,
            stage_index=0,
            tau=float(self.tau_init),
            run_kind="syntrop_initial",
            attempt_index=0,
            iter_idx=0,
            note="syntropization_initial_state",
        )

        abort_reason = self._syntropize_phase()
        if abort_reason is not None:
            return self._build_aborted_result(abort_reason)
        return self._relax_and_finalize()

    def _make_objfun(self, tau: float, relax: bool):
        return objfun(
            tau=tau,
            mu=self.mu,
            relax=relax,
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
            null_normalization_eps=self.null_normalization_eps,
            composition_penalty_multiplier=self.composition_penalty_multiplier,
        )

    def _evaluate_state_metrics(self, z_state: torch.Tensor, tau: float) -> dict:
        z_eval = self.projector.project(z_state)
        objfun_instance = self._make_objfun(tau=tau, relax=False)
        try:
            objfun_instance.fun_jac(z_eval)
        finally:
            self.energy_evaluations += int(getattr(objfun_instance, "energy_evaluations", 0))

        metrics = {}
        for key, value in objfun_instance.last_metrics.items():
            if value is None:
                metrics[key] = float("nan")
                continue
            try:
                metrics[key] = float(value)
            except (TypeError, ValueError):
                metrics[key] = float("nan")

        Xmix = self.projector._project_rows(z_eval[: self.n_x].view(self.n_mix, self.n_elements))
        mix_counts = Xmix.sum(dim=0)
        selected_counts = mix_counts + objfun_instance.fixed_selected_counts.to(dtype=Xmix.dtype)
        selected_total = torch.clamp(torch.sum(selected_counts), min=self.null_normalization_eps)
        metrics["composition_selected"] = (selected_counts / selected_total).detach().cpu().numpy()
        return metrics

    def _unpack_state(self, z_state: torch.Tensor):
        x = self.projector._project_rows(z_state[: self.n_x].view(self.n_mix, self.n_elements))
        if self.fix:
            frac = torch.tensor(self.atoms_base.get_scaled_positions(), device=self.device, dtype=torch.float32)
            uc = torch.tensor(cell_to_ucparams(np.asarray(self.atoms_base.cell.array, dtype=float)), device=self.device, dtype=torch.float32)
        else:
            uc = z_state[self.n_x:self.n_x + self.n_uc]
            frac = z_state[self.n_x + self.n_uc:].view(self.n_atoms, 3)
        return x, uc, frac

    def _objective_no_entropy_stagnated(self, f_true_history) -> tuple[bool, float]:
        window = int(self.escape_stagnation_window)
        if len(f_true_history) < window:
            return False, float("inf")
        recent = np.asarray(f_true_history[-window:], dtype=float)
        if not np.all(np.isfinite(recent)):
            return False, float("nan")
        max_step = float(np.max(np.abs(np.diff(recent))))
        span = float(abs(recent[-1] - recent[0]))
        change_metric = max(max_step, span)
        return bool(change_metric <= float(self.pbfgs_tol)), change_metric

    def _maybe_round_and_finish_stage(self) -> bool:
        Xmix = self.projector.project(self.z)[:self.n_x].view(self.n_mix, self.n_elements)
        ent_total, _, _ = null_total_entropy(
            Xmix,
            self.Tsallis_q,
            self.mix_is_null,
            eps=self.null_normalization_eps,
        )
        mix_atom_count = torch.clamp(torch.sum(Xmix), min=self.null_normalization_eps)
        ent = float((ent_total / mix_atom_count).detach().item())
        if ent < 1e-4:
            X_round = self._round_projected_weights(Xmix)
            self.z = torch.cat([X_round.view(-1), self.z[self.n_x:]], dim=0)
            self.log.info("Entropy small -> rounding null/species weights and moving to final relaxation.")
            return True
        return False

    def _reseed_dead_nulls_for_escape(self, z_state: torch.Tensor):
        if self.n_mix == 0 or (not torch.any(self.mix_is_null)):
            return z_state, 0, np.zeros(0, dtype=int)

        Xmix, uc, frac = self._unpack_state(z_state)
        occ_np = Xmix.sum(dim=1).detach().cpu().numpy()
        null_mask_np = self.mix_is_null.detach().cpu().numpy().astype(bool)
        dead_local = np.where(null_mask_np & (occ_np <= self.null_threshold))[0]
        if dead_local.size == 0:
            return z_state, 0, np.zeros(0, dtype=int)

        frac_np = frac.detach().cpu().numpy().copy()
        cell_np = ucparams_to_lattice(uc).detach().cpu().numpy()
        mix_atom_indices = self.mix_indices.detach().cpu().numpy()
        dead_atom_indices = mix_atom_indices[dead_local]
        keep_mask = np.ones(len(frac_np), dtype=bool)
        keep_mask[dead_atom_indices] = False
        blocked_frac = frac_np[keep_mask]
        reseeded_frac = self._escape_reseed_candidates(min_count=int(dead_local.size))
        repair = repair_null_overlaps(
            np.vstack([blocked_frac, frac_np[dead_atom_indices]]),
            cell_np,
            null_indices=np.arange(len(blocked_frac), len(blocked_frac) + len(dead_atom_indices), dtype=int),
            candidate_frac=reseeded_frac,
            min_distance=self.null_min_distance,
        )
        placed_frac = repair.frac_all[len(blocked_frac):]
        frac_np[dead_atom_indices] = placed_frac
        deactivated_mix_local = (
            dead_local[np.asarray(repair.deactivated_null_local_indices, dtype=int)]
            if repair.deactivated_null_local_indices
            else np.zeros(0, dtype=int)
        )

        frac_t = torch.tensor(frac_np, device=self.device, dtype=torch.float32)
        return (
            torch.cat([z_state[:self.n_x], z_state[self.n_x:self.n_x + self.n_uc], frac_t.reshape(-1)], dim=0),
            int(dead_local.size),
            np.asarray(deactivated_mix_local, dtype=int),
        )

    def _perturb_escape_state(self, z_state: torch.Tensor):
        z_trial, reseeded_nulls, deactivated_mix_local = self._reseed_dead_nulls_for_escape(z_state)
        Xmix, uc, frac = self._unpack_state(z_trial)

        x_rows = Xmix.detach().clone()
        if self.n_mix > 0:
            null_mask = self.mix_is_null
            if torch.any(~null_mask) and self.escape_species_prob_sigma > 0.0:
                normal_rows = torch.nonzero(~null_mask, as_tuple=False).view(-1)
                shocked = x_rows[normal_rows] + self.escape_species_prob_sigma * torch.randn_like(x_rows[normal_rows])
                x_rows[normal_rows] = self.projector._project_simplex_rows(shocked)

            if torch.any(null_mask):
                null_rows = torch.nonzero(null_mask, as_tuple=False).view(-1)
                species_probs = x_rows[null_rows].clone()
                empty_prob = torch.clamp(1.0 - species_probs.sum(dim=1), min=0.0, max=1.0)
                if self.escape_empty_prob_sigma > 0.0:
                    empty_prob = empty_prob + self.escape_empty_prob_sigma * torch.randn_like(empty_prob)
                    empty_prob = torch.clamp(empty_prob, min=0.0, max=1.0)
                if self.escape_species_prob_sigma > 0.0:
                    species_probs = species_probs + self.escape_species_prob_sigma * torch.randn_like(species_probs)
                    species_probs = torch.clamp(species_probs, min=0.0)

                target_occ = torch.clamp(1.0 - empty_prob, min=0.0, max=1.0)
                species_sum = species_probs.sum(dim=1, keepdim=True)
                scaled_species = torch.where(
                    species_sum > self.null_normalization_eps,
                    species_probs / torch.clamp(species_sum, min=self.null_normalization_eps) * target_occ.unsqueeze(1),
                    torch.full_like(species_probs, 1.0 / float(max(1, self.n_elements))) * target_occ.unsqueeze(1),
                )
                x_rows[null_rows] = scaled_species

        x_raw = x_rows.reshape(-1)
        if deactivated_mix_local.size > 0:
            x_rows = x_raw.view(self.n_mix, self.n_elements).clone()
            x_rows[torch.as_tensor(deactivated_mix_local, device=self.device, dtype=torch.long)] = 0.0
            x_raw = x_rows.reshape(-1)

        if self.fix:
            return x_raw, int(reseeded_nulls)

        uc_trial = uc.detach().clone()
        if self.escape_cell_length_sigma > 0.0:
            uc_trial[:3] = uc_trial[:3] + self.escape_cell_length_sigma * torch.randn(3, device=self.device, dtype=uc_trial.dtype)
        if self.escape_cell_angle_sigma_deg > 0.0:
            uc_trial[3:] = uc_trial[3:] + math.radians(self.escape_cell_angle_sigma_deg) * torch.randn(3, device=self.device, dtype=uc_trial.dtype)

        uc_trial = torch.max(torch.min(uc_trial, self.uc_upper), self.uc_lower)

        frac_np = frac.detach().cpu().numpy()
        cell_np = ucparams_to_lattice(uc_trial).detach().cpu().numpy()
        frac_trial = perturb_fractional_positions(frac_np, cell_np, self.escape_pos_sigma)

        x_rows = x_raw.view(self.n_mix, self.n_elements)
        occ_np = x_rows.sum(dim=1).detach().cpu().numpy() if self.n_mix > 0 else np.zeros(0, dtype=float)
        mix_atom_indices = self.mix_indices.detach().cpu().numpy()
        null_local_indices = np.where(self.mix_is_null.detach().cpu().numpy().astype(bool))[0]
        null_atom_indices = mix_atom_indices[null_local_indices]
        if null_atom_indices.size > 0 and self.null_min_distance > 0.0:
            priority = occ_np[null_local_indices] if occ_np.size > 0 else np.zeros(len(null_local_indices), dtype=float)
            repair = repair_null_overlaps(
                frac_trial,
                cell_np,
                null_indices=null_atom_indices,
                candidate_frac=self._escape_reseed_candidates(min_count=int(len(null_atom_indices))),
                min_distance=self.null_min_distance,
                priority_scores=priority,
            )
            frac_trial = repair.frac_all
            if repair.deactivated_null_local_indices:
                deactivate_mix_local = null_local_indices[np.asarray(repair.deactivated_null_local_indices, dtype=int)]
                x_rows = x_rows.clone()
                x_rows[torch.as_tensor(deactivate_mix_local, device=self.device, dtype=torch.long)] = 0.0
                x_raw = x_rows.reshape(-1)

        frac_t = torch.tensor(frac_trial, device=self.device, dtype=torch.float32)

        return torch.cat([x_raw, uc_trial, frac_t.reshape(-1)], dim=0), int(reseeded_nulls)

    def _run_pbfgs_stage(
        self,
        z0: torch.Tensor,
        tau: float,
        *,
        stage_index: int,
        run_kind: str,
        attempt_index: int,
        maxit: int,
        merge_history: bool,
    ):
        z_start = self.projector.project(z0)
        objfun_instance = self._make_objfun(tau=tau, relax=False)
        cost_history = []
        f_true_history = []
        formation_energy_history = []
        local_escape_events = []
        run_failed = False
        failure_reason = None

        def callback(z, f, err, k, converged):
            raw_cost = objfun_instance.last_metrics.get("cost", np.nan)
            raw_f_true = objfun_instance.last_metrics.get("objective_no_entropy", np.nan)
            raw_formation = objfun_instance.last_metrics.get("formation_energy", np.nan)
            cost = float("nan") if raw_cost is None else float(raw_cost)
            f_true = float("nan") if raw_f_true is None else float(raw_f_true)
            formation_energy = float("nan") if raw_formation is None else float(raw_formation)
            cost_history.append(cost)
            f_true_history.append(f_true)
            formation_energy_history.append(formation_energy)
            objfun_instance.pbfgs_callback(z, f, err, k, converged)

        opt = minimizer(
            fun_jac=objfun_instance.fun_jac,
            x0=z_start,
            args=(),
            project=self.projector.project,
            proj_args=(),
            callback=callback,
            BFGS_hist=5,
            maxit=maxit,
            tol=self.pbfgs_tol,
            max_line_search=5,
        )

        try:
            _, z_out = opt.PBFGS()
        except PBFGSPlateauReached:
            run_failed = False
            failure_reason = "plateau_reached"
            self.log.info(f"{run_kind} PBFGS reached the plateau stop criterion; keeping the latest iterate.")
            if getattr(objfun_instance, "z_last", None) is not None:
                z_out = objfun_instance.z_last.detach().clone()
            else:
                z_out = z_start
        except Exception as e:
            run_failed = True
            failure_reason = str(e)
            self.log.warning(f"{run_kind} PBFGS interrupted: {e}")
            if getattr(objfun_instance, "z_last", None) is not None:
                z_out = objfun_instance.z_last.detach().clone()
            else:
                z_out = z_start

        run_info = self._make_parallel_run_record(
            z_start=z_start,
            z_out=z_out,
            objfun_instance=objfun_instance,
            tau=tau,
            stage_index=stage_index,
            run_kind=run_kind,
            attempt_index=attempt_index,
            f_true_history=f_true_history,
            formation_energy_history=formation_energy_history,
            cost_history=cost_history,
            failed=run_failed,
            failure_reason=failure_reason,
        )
        self.energy_evaluations += int(run_info.get("energy_evaluations", 0))

        if merge_history:
            self.history = objfun_instance.history
            self.global_iter = int(self.history["iter"][-1]) if len(self.history.get("iter", [])) else self.global_iter
            self.cum_global_iter += self.global_iter
            for key, values in objfun_instance.all_history.items():
                if key not in self.all_history:
                    self.all_history[key] = []
                self.all_history[key].extend(values)
            if local_escape_events:
                self.all_history["escape_rounds"].extend(local_escape_events)

        return z_out, run_info

    def _try_escape_trials(self, z_start: torch.Tensor, tau: float, base_run_info: dict, anneal_index: int, escape_round: int):
        base_cost = float(base_run_info["cost_final"])
        base_formation = float(base_run_info["formation_energy_final"])
        base_entropy = float(base_run_info["entropy_per_atom_final"])
        base_penalty = float(base_run_info["composition_penalty_final"])
        base_objective = float(base_run_info["objective_no_entropy_final"])

        best_z = z_start
        best_cost = float(base_cost)
        best_objective = float(base_objective)
        best_trial = None
        best_accept_mode = None
        trial_records = []
        trial_states: dict[int, torch.Tensor] = {}
        objective_candidates = []

        for trial_idx in range(self.escape_trials):
            z_trial0, reseeded_nulls = self._perturb_escape_state(z_start)
            z_trial, run_info = self._run_pbfgs_stage(
                z_trial0,
                tau,
                stage_index=anneal_index,
                run_kind="escape_trial",
                attempt_index=trial_idx + 1,
                maxit=self.escape_trial_maxit,
                merge_history=False,
            )
            trial_index = int(trial_idx + 1)
            trial_states[trial_index] = z_trial
            trial_cost = float(run_info["cost_final"])
            trial_objective = float(run_info["objective_no_entropy_final"])
            trial_formation = float(run_info["formation_energy_final"])
            trial_entropy = float(run_info["entropy_per_atom_final"])
            trial_penalty = float(run_info["composition_penalty_final"])

            objective_accept_pass = (
                np.isfinite(best_objective)
                and np.isfinite(trial_objective)
                and trial_objective < best_objective - float(self.escape_accept_tol)
            )

            self.log.info(
                f"[escape anneal={anneal_index + 1:02d} round={escape_round}] "
                f"trial {trial_idx + 1}/{self.escape_trials} | "
                f"cost={trial_cost:+.6f} | "
                f"formation_energy={trial_formation:+.6f} eV/atom | "
                f"objective_no_entropy={trial_objective:+.6f} | "
                f"penalty={trial_penalty:+.6f} | "
                f"ent/atom={trial_entropy:+.6f} | "
                f"objective_no_entropy_lowered={objective_accept_pass} | "
                f"reseeded_nulls={reseeded_nulls}"
            )
            record = {
                "trial_index": trial_index,
                "anneal": int(anneal_index + 1),
                "escape_round": int(escape_round),
                "reseeded_nulls": int(reseeded_nulls),
                "shock_state": self._state_snapshot_from_z(
                    z_trial0,
                    stage_index=anneal_index,
                    tau=tau,
                    run_kind="escape_trial_shock",
                    attempt_index=trial_index,
                    iter_idx=0,
                    note="post_shock_pre_relaxation",
                    reseeded_nulls=int(reseeded_nulls),
                ),
                "relaxation": run_info,
                "cost_final": float(trial_cost),
                "formation_energy_final": float(trial_formation),
                "objective_no_entropy_final": float(trial_objective),
                "composition_penalty_final": float(trial_penalty),
                "entropy_per_atom_final": float(trial_entropy),
                "objective_accept_pass": bool(objective_accept_pass),
                "objective_accept_tol": float(self.escape_accept_tol),
                "accepted": False,
            }
            trial_records.append(record)
            if objective_accept_pass:
                objective_candidates.append(record)

        selected_record = None
        if objective_candidates:
            selected_record = min(
                objective_candidates,
                key=lambda rec: (
                    float(rec["objective_no_entropy_final"]),
                    int(rec["trial_index"]),
                ),
            )
            best_accept_mode = "objective_no_entropy_lowering"

        if selected_record is not None:
            best_trial = int(selected_record["trial_index"])
            best_z = trial_states[best_trial]
            best_cost = float(selected_record["cost_final"])
            best_objective = float(selected_record["objective_no_entropy_final"])
            for record in trial_records:
                if record["trial_index"] == int(best_trial):
                    record["accepted"] = True

        return {
            "accepted": selected_record is not None,
            "best_z": best_z,
            "best_cost": best_cost,
            "best_formation_energy": (
                float(selected_record["formation_energy_final"]) if selected_record is not None else float(base_formation)
            ),
            "best_objective_no_entropy": (
                float(selected_record["objective_no_entropy_final"]) if selected_record is not None else float(base_objective)
            ),
            "best_entropy_per_atom": (
                float(selected_record["entropy_per_atom_final"]) if selected_record is not None else float(base_entropy)
            ),
            "best_composition_penalty": (
                float(selected_record["composition_penalty_final"]) if selected_record is not None else float(base_penalty)
            ),
            "best_trial": best_trial,
            "accept_mode": best_accept_mode,
            "trial_records": trial_records,
        }

    def _syntropize_phase(self):
        self.cum_global_iter = 0
        last_escape_anneal = None
        for anneal in range(self.n_anneals):
            self.anneal_boundaries.append(self.global_iter)
            tau = self.tau_init * (self.tau_mult ** anneal)
            stage_record = {
                "stage": int(anneal + 1),
                "tau": float(tau),
                "optimization_runs": [],
                "escape_rounds": [],
            }

            self.log.info("")
            self.log.info("=" * 80)
            self.log.info(f"ANNEALING STEP {anneal + 1}/{self.n_anneals} | tau = {tau:.6e}")
            self.log.info("=" * 80)

            escape_attempts_used = 0
            stage_failed = False
            while True:
                self.z, run_info = self._run_pbfgs_stage(
                    self.z,
                    tau,
                    stage_index=anneal,
                    run_kind="main" if escape_attempts_used == 0 else "post_escape",
                    attempt_index=escape_attempts_used,
                    maxit=self.first_stage_maxit if anneal == 0 and escape_attempts_used == 0 else self.stage_maxit,
                    merge_history=True,
                )
                stage_record["accepted"] = True
                stage_record["optimization_runs"].append(run_info)
                if run_info.get("failed", False):
                    stage_record["accepted"] = False
                    stage_failed = True
                    self.all_history["stages"].append(stage_record)
                    self.log.warning(
                        f"[stage anneal={anneal + 1:02d}] PBFGS failed "
                        f"({run_info.get('failure_reason', 'unknown error')}). Moving to the next tau stage."
                    )
                    break

                if not self.escape_enabled:
                    break

                if not run_info.get("stagnated", False):
                    self.log.info(
                        f"[escape anneal={anneal + 1:02d}] objective_no_entropy still changing "
                        f"(change={run_info.get('stagnation_change', float('nan')):.3e}, "
                        f"tol={run_info.get('stagnation_tol', self.pbfgs_tol):.3e}); skipping shock."
                    )
                    break

                shock_blocked_by_interval = (
                    escape_attempts_used == 0
                    and last_escape_anneal is not None
                    and (anneal - last_escape_anneal) < self.escape_min_shock_interval
                )
                if shock_blocked_by_interval:
                    stages_since = int(anneal - last_escape_anneal)
                    stages_remaining = int(self.escape_min_shock_interval - stages_since)
                    self.log.info(
                        f"[escape anneal={anneal + 1:02d}] objective_no_entropy stagnated, "
                        f"but shock cooldown is active "
                        f"({stages_since}/{self.escape_min_shock_interval} stages since last shock episode; "
                        f"next allowed in {stages_remaining} stage(s)). Skipping shock."
                    )
                    break

                if escape_attempts_used >= self.escape_max_attempts_per_stage:
                    self.log.info(
                        f"[escape anneal={anneal + 1:02d}] reached maximum accepted escape rounds "
                        f"for this stage ({self.escape_max_attempts_per_stage})."
                    )
                    break

                escape_round = int(escape_attempts_used + 1)
                self.log.info(
                    f"[escape anneal={anneal + 1:02d}] objective_no_entropy stagnated over the last "
                    f"{run_info.get('stagnation_window', 3)} iterations "
                    f"(change={run_info.get('stagnation_change', float('nan')):.3e}, "
                    f"tol={run_info.get('stagnation_tol', self.pbfgs_tol):.3e}). "
                    f"Running escape round {escape_round}/{self.escape_max_attempts_per_stage}."
                )
                escape_result = self._try_escape_trials(
                    self.z,
                    tau,
                    base_run_info=run_info,
                    anneal_index=anneal,
                    escape_round=escape_round,
                )
                if escape_attempts_used == 0:
                    last_escape_anneal = anneal

                escape_record = {
                    "anneal": int(anneal + 1),
                    "escape_round": escape_round,
                    "base_cost": float(run_info["cost_final"]),
                    "base_formation_energy": float(run_info["formation_energy_final"]),
                    "base_objective_no_entropy": float(run_info["objective_no_entropy_final"]),
                    "base_entropy_per_atom": float(run_info["entropy_per_atom_final"]),
                    "base_composition_penalty": float(run_info["composition_penalty_final"]),
                    "accepted": bool(escape_result["accepted"]),
                    "best_trial": None if escape_result["best_trial"] is None else int(escape_result["best_trial"]),
                    "best_cost": float(escape_result["best_cost"]),
                    "best_formation_energy": float(escape_result["best_formation_energy"]),
                    "best_objective_no_entropy": float(escape_result["best_objective_no_entropy"]),
                    "best_entropy_per_atom": float(escape_result["best_entropy_per_atom"]),
                    "best_composition_penalty": float(escape_result["best_composition_penalty"]),
                    "accept_mode": None if escape_result["accept_mode"] is None else str(escape_result["accept_mode"]),
                    "trial_records": escape_result["trial_records"],
                }
                self.all_history["escape_rounds"].append(dict(escape_record))
                stage_record["escape_rounds"].append(escape_record)

                if not escape_result["accepted"]:
                    self.log.info(
                        f"[escape anneal={anneal + 1:02d}] no perturbation trial lowered "
                        f"objective_no_entropy by at least {self.escape_accept_tol:.3e}."
                    )
                    break

                self.z = escape_result["best_z"]
                escape_attempts_used += 1
                self.log.info(
                    f"[escape anneal={anneal + 1:02d}] accepted trial {escape_result['best_trial']} "
                    f"via {escape_result['accept_mode']} | "
                    f"objective_no_entropy={escape_result['best_objective_no_entropy']:+.6f} | "
                    f"formation_energy={escape_result['best_formation_energy']:+.6f} eV/atom | "
                    f"penalty={escape_result['best_composition_penalty']:+.6f} | "
                    f"ent/atom={escape_result['best_entropy_per_atom']:+.6f} | "
                    f"cost={escape_result['best_cost']:+.6f}. Re-running PBFGS in the same anneal stage."
                )

            if stage_failed:
                if self._maybe_round_and_finish_stage():
                    break
                continue

            stage_metrics = self._evaluate_state_metrics(self.z, tau=tau)
            stage_record["stage_end_metrics"] = {
                "cost": float(stage_metrics["cost"]),
                "formation_energy": float(stage_metrics["formation_energy"]),
                "energy_above_hull": float(stage_metrics["energy_above_hull"]),
                "objective_no_entropy": float(stage_metrics["objective_no_entropy"]),
                "entropy_per_atom": float(stage_metrics["entropy_per_atom"]),
            }
            self.log.info(
                f"[stage anneal={anneal + 1:02d}] end metrics | "
                f"cost={stage_metrics['cost']:+.6f}, "
                f"E_above_hull={stage_metrics['energy_above_hull']:+.6f} eV/atom, "
                f"E_form={stage_metrics['formation_energy']:+.6f} eV/atom"
            )

            if (
                self.abort_energy_above_hull is not None
                and np.isfinite(stage_metrics["energy_above_hull"])
                and stage_metrics["energy_above_hull"] > self.abort_energy_above_hull
            ):
                reason = (
                    f"stage {anneal + 1} energy_above_hull={stage_metrics['energy_above_hull']:.6f} "
                    f"exceeded abort threshold {self.abort_energy_above_hull:.6f}"
                )
                stage_record["aborted"] = True
                stage_record["abort_reason"] = reason
                self.all_history["stages"].append(stage_record)
                self.all_history["aborted"] = True
                self.all_history["abort_reason"] = reason
                self.log.warning(f"[stage anneal={anneal + 1:02d}] aborting trial: {reason}")
                return reason

            self.all_history["stages"].append(stage_record)
            if self._maybe_round_and_finish_stage():
                break
        return None

    def _build_symbols(self, Xmix_np: np.ndarray, site_occupancies_mix: np.ndarray) -> list[str]:
        symbols = [None] * self.n_atoms
        base_symbols = list(self.atoms_base.get_chemical_symbols())
        for idx in range(self.n_atoms):
            if not bool(self.atoms_base.arrays["site_is_mix"][idx]):
                symbols[idx] = base_symbols[idx]

        elem_idx = np.argmax(Xmix_np, axis=1) if self.n_mix > 0 else np.zeros(0, dtype=int)
        mix_syms = [self.selected_elements[i] for i in elem_idx]
        mix_rows = self.mix_indices.detach().cpu().numpy().tolist()
        mix_is_null_np = self.mix_is_null.detach().cpu().numpy().astype(bool)
        for k, atom_idx in enumerate(mix_rows):
            if mix_is_null_np[k] and site_occupancies_mix[k] <= self.null_threshold:
                symbols[atom_idx] = "X"
            else:
                symbols[atom_idx] = mix_syms[k]
        return symbols

    def _build_aborted_result(self, reason: str) -> OptimizationResult:
        metrics = self._evaluate_state_metrics(self.z, tau=0.0)
        z_proj = self.projector.project(self.z)
        Xmix = z_proj[:self.n_x].view(self.n_mix, self.n_elements)
        Xmix_np = Xmix.detach().cpu().numpy()
        site_occupancies_mix = np.sum(Xmix_np, axis=1)
        _, uc_state, frac_state = self._unpack_state(z_proj)
        lattice = ucparams_to_lattice(uc_state).detach().cpu().numpy()
        frac_np = frac_state.detach().cpu().numpy()
        composition_selected = np.asarray(metrics["composition_selected"], dtype=float)
        symbols = self._build_symbols(Xmix_np, site_occupancies_mix)
        label = " ".join([f"{el}{composition_selected[i]:.3f}" for i, el in enumerate(self.selected_elements)])

        self.all_history["aborted"] = True
        self.all_history["abort_reason"] = str(reason)
        self.all_history["final_relaxation"] = None
        self.all_history["selected_elements"] = list(self.selected_elements)
        self.all_history["species_order"] = list(self.species_order)
        self.all_history["site_is_mix"] = self.atoms_base.arrays["site_is_mix"].astype(int).tolist()
        self.all_history["site_is_null"] = self.atoms_base.arrays["site_is_null"].astype(int).tolist()
        self.all_history["site_fixed_Z"] = self.atoms_base.arrays["site_fixed_Z"].astype(int).tolist()
        self.all_history["atom_symbols_placeholder"] = list(self.atoms_base.get_chemical_symbols())
        self.all_history["final_symbols"] = list(symbols)
        self.all_history["final_site_occupancies_mix"] = site_occupancies_mix.tolist()
        self.all_history["final_selected_composition"] = composition_selected.tolist()
        self.all_history["final_energy_evaluations"] = int(self.energy_evaluations)

        return OptimizationResult(
            history=self.history,
            all_history=self.all_history,
            anneal_boundaries=self.anneal_boundaries,
            final_energy_total_ev=float(metrics["total_energy"]),
            final_formation_energy_ev_per_atom=float(metrics["formation_energy"]),
            final_energy_above_hull_ev_per_atom=float(metrics["energy_above_hull"]),
            selected_elements=self.selected_elements,
            X_final=Xmix_np,
            site_occupancies=site_occupancies_mix,
            composition_selected=composition_selected,
            frac_coords=frac_np,
            lattice_vectors=lattice,
            symbols=symbols,
            db_label=label,
            energy_evaluations=int(self.energy_evaluations),
            cum_global_iter=self.cum_global_iter,
            aborted=True,
            abort_reason=str(reason),
        )

    def _relax_and_finalize(self) -> OptimizationResult:
        objfun_instance = self._make_objfun(tau=0.0, relax=True)
        z_relax_start = self.z.detach().clone()
        final_cost_history = []
        final_f_true_history = []
        final_formation_history = []

        if self.fix is False:
            self.log.info("")
            self.log.info("=" * 80)
            self.log.info("FINAL PBFGS RELAXATION")
            self.log.info("=" * 80)

            def final_callback(z, f, err, k, converged):
                final_cost_history.append(float(objfun_instance.last_metrics.get("cost", np.nan)))
                final_f_true_history.append(float(objfun_instance.last_metrics.get("objective_no_entropy", np.nan)))
                final_formation_history.append(float(objfun_instance.last_metrics.get("formation_energy", np.nan)))
                objfun_instance.pbfgs_callback(z, f, err, k, converged)

            opt = minimizer(
                fun_jac=objfun_instance.fun_jac,
                x0=self.z,
                args=(),
                project=self.projector.project,
                proj_args=(),
                callback=final_callback,
                BFGS_hist=5,
                maxit=self.pbfgs_maxit_relax,
                tol=self.pbfgs_tol_relax,
                max_line_search=5,
            )

            try:
                _, self.z = opt.PBFGS()
            except PBFGSPlateauReached:
                self.log.info("Final relaxation reached the plateau stop criterion; keeping the latest iterate.")
                if getattr(objfun_instance, "z_last", None) is not None:
                    self.z = objfun_instance.z_last.detach().clone()
            except Exception as e:
                self.log.warning(f"Final relaxation interrupted: {e}")
                if getattr(objfun_instance, "z_last", None) is not None:
                    self.z = objfun_instance.z_last.detach().clone()

            uc_final = self.z[self.n_x:self.n_x + self.n_uc]
            frac_final = self.z[self.n_x + self.n_uc:].view(self.n_atoms, 3)
        else:
            atoms = self.atoms_base.copy()
            frac_final = torch.tensor(atoms.get_scaled_positions(), device=self.device)
            lattice = np.array(atoms.get_cell().array, dtype=float)
            uc_final = torch.tensor(cell_to_ucparams(lattice), device=self.device).detach().clone()

        lattice = ucparams_to_lattice(uc_final).detach().cpu().numpy()
        frac_np = frac_final.detach().cpu().numpy()

        Xf_mix = self.projector.project(self.z)[:self.n_x].view(self.n_mix, self.n_elements)
        Xf_mix_np = Xf_mix.detach().cpu().numpy()
        site_occupancies_mix = np.sum(Xf_mix_np, axis=1)

        symbols = [None] * self.n_atoms
        base_symbols = list(self.atoms_base.get_chemical_symbols())
        for idx in range(self.n_atoms):
            if not bool(self.atoms_base.arrays["site_is_mix"][idx]):
                symbols[idx] = base_symbols[idx]

        elem_idx = np.argmax(Xf_mix_np, axis=1) if self.n_mix > 0 else np.zeros(0, dtype=int)
        mix_syms = [self.selected_elements[i] for i in elem_idx]
        mix_rows = self.mix_indices.detach().cpu().numpy().tolist()
        mix_is_null_np = self.mix_is_null.detach().cpu().numpy().astype(bool)
        for k, atom_idx in enumerate(mix_rows):
            if mix_is_null_np[k] and site_occupancies_mix[k] <= self.null_threshold:
                symbols[atom_idx] = "X"
            else:
                symbols[atom_idx] = mix_syms[k]

        n_fairchem_elements = 100
        eat_fairchem = torch.zeros((self.n_atoms, n_fairchem_elements), device=self.device, dtype=Xf_mix.dtype)

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
            mix_counts = Xf_mix.sum(dim=0)
            total_atom_count = torch.clamp(
                torch.tensor(
                    float(len(self.fixed_indices)) + float(torch.sum(mix_counts).item()),
                    device=self.device,
                    dtype=Xf_mix.dtype,
                ),
                min=self.null_normalization_eps,
            )
            ref_total = (
                torch.dot(mix_counts.to(objfun_instance.ref_energies_tensor.dtype), objfun_instance.ref_energies_tensor)
                + objfun_instance.fixed_ref_energy_sum.to(objfun_instance.ref_energies_tensor.dtype)
            )
            E_total_t = torch.tensor(E_total, device=self.device, dtype=objfun_instance.ref_energies_tensor.dtype)
            E_form_t = (E_total_t - ref_total) / total_atom_count.to(dtype=objfun_instance.ref_energies_tensor.dtype)
            E_form = float(E_form_t.item())

            selected_counts = mix_counts + objfun_instance.fixed_selected_counts.to(dtype=Xf_mix.dtype)
            selected_total = torch.clamp(torch.sum(selected_counts), min=self.null_normalization_eps)
            composition_selected = (selected_counts / selected_total).detach().cpu().numpy()

            E_hull, _ = calc_hull(
                Xf_mix,
                objfun_instance.hullcalc,
                null_mask=self.mix_is_null,
                fixed_counts=objfun_instance.fixed_selected_counts,
                min_total=self.null_normalization_eps,
            )
            E_above = float(E_form - E_hull)

        label = " ".join([f"{el}{composition_selected[i]:.3f}" for i, el in enumerate(self.selected_elements)])

        self.history = objfun_instance.history
        for k, v in objfun_instance.all_history.items():
            if k not in self.all_history:
                self.all_history[k] = []
            self.all_history[k].extend(v)

        self.all_history["final_relaxation"] = self._make_parallel_run_record(
            z_start=z_relax_start,
            z_out=self.z,
            objfun_instance=objfun_instance,
            tau=0.0,
            stage_index=len(self.anneal_boundaries),
            run_kind="final_relaxation",
            attempt_index=0,
            f_true_history=final_f_true_history,
            formation_energy_history=final_formation_history,
            cost_history=final_cost_history,
        )
        self.energy_evaluations += int(getattr(objfun_instance, "energy_evaluations", 0))

        self.all_history["selected_elements"] = list(self.selected_elements)
        self.all_history["species_order"] = list(self.species_order)
        self.all_history["site_is_mix"] = self.atoms_base.arrays["site_is_mix"].astype(int).tolist()
        self.all_history["site_is_null"] = self.atoms_base.arrays["site_is_null"].astype(int).tolist()
        self.all_history["site_fixed_Z"] = self.atoms_base.arrays["site_fixed_Z"].astype(int).tolist()
        self.all_history["atom_symbols_placeholder"] = list(self.atoms_base.get_chemical_symbols())
        self.all_history["final_symbols"] = list(symbols)
        self.all_history["final_site_occupancies_mix"] = site_occupancies_mix.tolist()
        self.all_history["final_selected_composition"] = composition_selected.tolist()
        self.all_history["final_energy_evaluations"] = int(self.energy_evaluations)

        return OptimizationResult(
            history=self.history,
            all_history=self.all_history,
            anneal_boundaries=self.anneal_boundaries,
            final_energy_total_ev=float(E_total),
            final_formation_energy_ev_per_atom=float(E_form),
            final_energy_above_hull_ev_per_atom=float(E_above),
            selected_elements=self.selected_elements,
            X_final=Xf_mix_np,
            site_occupancies=site_occupancies_mix,
            composition_selected=composition_selected,
            frac_coords=frac_np,
            lattice_vectors=lattice,
            symbols=symbols,
            db_label=label,
            energy_evaluations=int(self.energy_evaluations),
            cum_global_iter=self.cum_global_iter,
        )
