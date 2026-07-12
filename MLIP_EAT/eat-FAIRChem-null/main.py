#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

# Avoid PyTorch scanning entry-point backends at import time.
os.environ.setdefault("TORCH_DEVICE_BACKEND_AUTOLOAD", "0")

import torch
import yaml

from fairchem.core import FAIRChemCalculator, pretrained_mlip
from modules.atoms_helpers import build_base_atoms
from modules.energy_above_hull import HullEnergyCalculator, PhaseRecord
from modules.logging_utils import setup_logging
from modules.plotting import plot_phase_diagram_from_db
from modules.seeding import load_rng_state, save_rng_state
from modules.syntropizer import OptimizationResult, syntropizer


def _db_label_from_result(result: OptimizationResult) -> str:
    if getattr(result, "db_label", None):
        return result.db_label
    comp = np.asarray(result.composition_selected, dtype=float)
    return " ".join(f"{el}{comp[i]:.3f}" for i, el in enumerate(result.selected_elements))


def _composition_tag(elements: list[str], composition: list[float]) -> str:
    comp = np.asarray(composition, dtype=float)
    return "_".join(f"{el}{comp[i]:.3f}" for i, el in enumerate(elements))


def _fmt_energy(value: float) -> str:
    try:
        value = float(value)
    except (TypeError, ValueError):
        return "nan"
    if not np.isfinite(value):
        return "nan"
    return f"{value:+.6f}"


def _append_progress_status(path: Path, message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with path.open("a", encoding="utf-8") as f:
        f.write(f"{timestamp} | {message}\n")


def ensure_endmembers_in_db_from_yaml(
    db_path: Path,
    design_space: list[str],
    endmembers: list[dict[str, Any]],
    *,
    overwrite: bool = False,
) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if not db_path.exists():
        db_path.write_text("{}", encoding="utf-8")

    db = json.loads(db_path.read_text(encoding="utf-8") or "{}")

    for em in endmembers:
        fractions = [float(x) for x in em["fractions"]]
        e_form = float(em["formation_energy_per_atom"])
        label = " ".join([f"{design_space[j]}{fractions[j]:.3f}" for j in range(len(design_space))])

        if (label in db) and (not overwrite):
            continue

        db[label] = {
            "elements": design_space,
            "fractions": fractions,
            "formation_energy_per_atom": e_form,
        }

    db_path.write_text(json.dumps(db, indent=2), encoding="utf-8")


def _automatic_elemental_endmembers(selected_elements: list[str]) -> list[dict[str, Any]]:
    endmembers = []
    n_elements = len(selected_elements)
    for i in range(n_elements):
        fractions = [0.0] * n_elements
        fractions[i] = 1.0
        endmembers.append(
            {
                "fractions": fractions,
                "formation_energy_per_atom": 0.0,
                "metadata": {"source": "automatic_elemental_endmember"},
            }
        )
    return endmembers


def _validate_endmember(em: Any, selected_elements: list[str], *, label: str) -> dict[str, Any]:
    if not isinstance(em, dict):
        raise ValueError(f"config: {label} must be a dict, got {type(em)}")
    if "fractions" not in em or "formation_energy_per_atom" not in em:
        raise ValueError(f"config: {label} needs 'fractions' and 'formation_energy_per_atom'")

    fr = em["fractions"]
    if not isinstance(fr, list) or len(fr) != len(selected_elements):
        raise ValueError(
            f"config: {label}.fractions must have length {len(selected_elements)} "
            f"(selected_elements order). Got: {fr}"
        )
    fractions = [float(v) for v in fr]
    s = float(sum(fractions))
    if abs(s - 1.0) > 1e-6:
        raise ValueError(f"config: {label}.fractions must sum to 1. Got sum={s} fr={fr}")

    out = dict(em)
    out["fractions"] = fractions
    out["formation_energy_per_atom"] = float(em["formation_energy_per_atom"])
    return out


def _merge_endmembers(selected_elements: list[str], yaml_endmembers: Any) -> list[dict[str, Any]]:
    endmembers = _automatic_elemental_endmembers(selected_elements)

    if yaml_endmembers is not None:
        if not isinstance(yaml_endmembers, list):
            raise ValueError("config: endmembers must be a list when provided")
        for i, em in enumerate(yaml_endmembers):
            endmembers.append(_validate_endmember(em, selected_elements, label=f"endmembers[{i}]"))

    deduped: list[dict[str, Any]] = []
    seen: set[tuple[float, ...]] = set()
    for em in endmembers:
        key = tuple(round(float(v), 12) for v in em["fractions"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(em)
    return deduped


@dataclass(frozen=True)
class RunInputs:
    out_dir: Path
    n_trials: int
    compositions: list[list[float]]
    fix: bool
    alpha: float
    Tsallis_q: float
    tau_init: float
    tau_mult: float
    first_stage_maxit: int
    stage_maxit: int
    pbfgs_tol: float
    composition_penalty_multiplier: float
    selected_elements: list[str]
    lattice_vectors: list[list[float]]
    frac_coords_init: Any
    supercell_size: tuple[int, int, int]
    db_path: Path
    reference_json_path: Path
    endmembers: list[dict[str, Any]]
    null_default_for_mix: bool
    null_num_sites: int | None
    null_min_distance: float
    null_initial_occupancy_sum: float | None
    null_normalization_eps: float
    null_threshold: float
    cell_init_mode: str
    random_cell_length_min: float
    random_cell_length_max: float
    random_cell_angle_min_deg: float
    random_cell_angle_max_deg: float
    cell_bound_length_min: float
    cell_bound_length_max: float
    cell_bound_angle_min_deg: float
    cell_bound_angle_max_deg: float
    escape_enabled: bool
    escape_trials: int
    escape_pos_sigma: float
    escape_empty_prob_sigma: float
    escape_species_prob_sigma: float
    escape_cell_length_sigma: float
    escape_cell_angle_sigma_deg: float
    escape_trial_maxit: int
    escape_stagnation_window: int
    escape_min_shock_interval: int
    escape_max_attempts_per_stage: int
    escape_accept_tol: float
    population_enabled: bool
    population_max_generations: int
    population_initial_runs: int
    population_parents_per_generation: int
    population_children_per_parent: int
    population_keep_archive: int
    population_memory_strength: float
    population_selection_metric: str
    population_initial_lattice_vectors: list[list[float]] | None
    population_seed_parent_xsf: Path | None
    population_seed_parent_symbols: list[str] | None
    abort_energy_above_hull: float | None
    seed: int | None


def _load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping/object. Got: {type(data)}")
    return data


def _as_mapping(section: Any, *, name: str) -> dict[str, Any]:
    if section is None:
        return {}
    if isinstance(section, dict):
        return dict(section)
    if isinstance(section, list):
        out: dict[str, Any] = {}
        for i, item in enumerate(section):
            if not isinstance(item, dict):
                raise ValueError(f"config: {name}[{i}] must be a dict, got {type(item)}")
            out.update(item)
        return out
    raise ValueError(f"config: {name} must be a dict or list[dict], got {type(section)}")


def _normalize_composition(comp: Any, *, n_elements: int, label: str) -> list[float]:
    if not isinstance(comp, list):
        raise ValueError(f"config: {label} must be a list, got {type(comp)}")
    if len(comp) != n_elements:
        raise ValueError(f"config: {label} must have length {n_elements}, got {len(comp)}: {comp}")

    arr = np.asarray([float(v) for v in comp], dtype=float)
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"config: {label} contains non-finite values: {comp}")
    if np.any(arr < -1e-10):
        raise ValueError(f"config: {label} contains negative fractions: {comp}")

    arr = np.clip(arr, 0.0, None)
    s = float(arr.sum())
    if s <= 0.0:
        raise ValueError(f"config: {label} must have a positive sum: {comp}")
    if abs(s - 1.0) > 1e-6:
        raise ValueError(f"config: {label} must sum to 1.0, got sum={s}: {comp}")
    arr = arr / s
    return arr.tolist()


def _dedupe_compositions(compositions: list[list[float]]) -> list[list[float]]:
    seen: set[tuple[float, ...]] = set()
    out: list[list[float]] = []
    for comp in compositions:
        key = tuple(round(float(v), 12) for v in comp)
        if key in seen:
            continue
        seen.add(key)
        out.append([float(v) for v in comp])
    return out


def _optional_lattice_vectors(value: Any, *, label: str) -> list[list[float]] | None:
    if value is None:
        return None
    if not isinstance(value, list) or len(value) != 3:
        raise ValueError(f"config: {label} must be a 3x3 list of lattice vectors")
    out = []
    for i, row in enumerate(value):
        if not isinstance(row, list) or len(row) != 3:
            raise ValueError(f"config: {label}[{i}] must be a length-3 list")
        out.append([float(x) for x in row])
    return out


def _resolve_optional_path(value: Any, *, infile: Path, label: str) -> Path | None:
    if value is None:
        return None
    p = Path(str(value)).expanduser()
    p = p if p.is_absolute() else (infile.parent / p)
    p = p.resolve()
    if not p.exists():
        raise FileNotFoundError(f"config: {label} not found: {p}")
    return p


def _read_xsf_structure(path: Path) -> tuple[list[list[float]], list[str], list[list[float]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        primvec_idx = lines.index("PRIMVEC")
        primcoord_idx = lines.index("PRIMCOORD")
    except ValueError as exc:
        raise ValueError(f"XSF file must contain PRIMVEC and PRIMCOORD sections: {path}") from exc

    lattice = []
    for j in range(3):
        parts = lines[primvec_idx + 1 + j].split()
        if len(parts) < 3:
            raise ValueError(f"Bad PRIMVEC row {j + 1} in {path}")
        lattice.append([float(parts[0]), float(parts[1]), float(parts[2])])

    header = lines[primcoord_idx + 1].split()
    if not header:
        raise ValueError(f"Bad PRIMCOORD header in {path}")
    n_atoms = int(header[0])

    symbols = []
    cart = []
    for line in lines[primcoord_idx + 2 : primcoord_idx + 2 + n_atoms]:
        parts = line.split()
        if len(parts) < 4:
            raise ValueError(f"Bad PRIMCOORD atom row in {path}: {line!r}")
        symbols.append(parts[0])
        cart.append([float(parts[1]), float(parts[2]), float(parts[3])])

    if len(symbols) != n_atoms:
        raise ValueError(f"PRIMCOORD expected {n_atoms} atoms but found {len(symbols)} in {path}")

    lattice_np = np.asarray(lattice, dtype=float)
    cart_np = np.asarray(cart, dtype=float)
    frac_np = np.linalg.solve(lattice_np.T, cart_np.T).T
    frac_np = np.mod(frac_np, 1.0)
    return lattice, symbols, frac_np.tolist()


def _integer_simplex_tuples(total: int, dims: int):
    if dims == 1:
        yield (int(total),)
        return
    for i in range(int(total) + 1):
        for rest in _integer_simplex_tuples(int(total) - i, dims - 1):
            yield (int(i),) + rest


def _simplex_grid_compositions(
    *,
    n_elements: int,
    subdivisions: int,
    include_boundary: bool,
    include_vertices: bool,
) -> list[list[float]]:
    subdivisions = int(subdivisions)
    if subdivisions <= 0:
        raise ValueError("config: composition_grid.subdivisions must be > 0")

    comps: list[list[float]] = []
    for counts in _integer_simplex_tuples(subdivisions, n_elements):
        counts_arr = np.asarray(counts, dtype=float)
        nonzero = int(np.count_nonzero(counts_arr))
        if not include_boundary and nonzero < n_elements:
            continue
        if not include_vertices and nonzero == 1:
            continue
        comps.append((counts_arr / float(subdivisions)).tolist())

    if not comps:
        raise ValueError("config: composition_grid produced no compositions; relax boundary/vertex filters")
    return comps


def _random_simplex_compositions(
    *,
    n_elements: int,
    n_points: int,
    seed: int | None,
    alpha: Any = 1.0,
    include_vertices: bool = False,
) -> list[list[float]]:
    n_points = int(n_points)
    if n_points < 0:
        raise ValueError("config: composition_random.n_points must be >= 0")

    if isinstance(alpha, list):
        if len(alpha) != n_elements:
            raise ValueError(f"config: composition_random.alpha must have length {n_elements}")
        alpha_arr = np.asarray([float(v) for v in alpha], dtype=float)
    else:
        alpha_arr = np.full(n_elements, float(alpha), dtype=float)
    if np.any(alpha_arr <= 0.0):
        raise ValueError("config: composition_random.alpha entries must be > 0")

    rng = np.random.default_rng(seed)
    comps: list[list[float]] = []
    if include_vertices:
        comps.extend(np.eye(n_elements, dtype=float).tolist())
    if n_points > 0:
        comps.extend(rng.dirichlet(alpha_arr, size=n_points).tolist())

    comps = _dedupe_compositions(comps)
    if not comps:
        raise ValueError("config: composition_random produced no compositions")
    return comps


def _load_compositions(cfg: dict[str, Any], selected_elements: list[str]) -> list[list[float]]:
    n_elements = len(selected_elements)
    if n_elements < 2:
        raise ValueError("config: selected_elements must contain at least two elements")

    has_manual = "compositions" in cfg
    has_grid = "composition_grid" in cfg
    has_random = "composition_random" in cfg
    mode_raw = cfg.get("composition_mode", None)

    if mode_raw is None:
        if has_grid and not has_manual and not has_random:
            mode = "grid"
        elif has_random and not has_manual and not has_grid:
            mode = "random"
        else:
            mode = "manual"
    else:
        mode = str(mode_raw).strip().lower()

    if mode not in {"manual", "grid", "random"}:
        raise ValueError("config: composition_mode must be one of {'manual', 'grid', 'random'}")

    if mode == "manual":
        compositions = cfg.get("compositions", None)
        if not isinstance(compositions, list) or not compositions:
            raise ValueError("config: compositions must be a non-empty list of lists")
        return _dedupe_compositions(
            [
                _normalize_composition(comp, n_elements=n_elements, label=f"compositions[{i}]")
                for i, comp in enumerate(compositions)
            ]
        )

    if mode == "grid":
        grid_cfg = _as_mapping(cfg.get("composition_grid", None), name="composition_grid")
        if "subdivisions" not in grid_cfg:
            raise ValueError("config: composition_grid.subdivisions is required")
        include_boundary = bool(grid_cfg.get("include_boundary", grid_cfg.get("include_edges", True)))
        include_vertices = bool(grid_cfg.get("include_vertices", include_boundary))
        return _dedupe_compositions(
            _simplex_grid_compositions(
                n_elements=n_elements,
                subdivisions=int(grid_cfg["subdivisions"]),
                include_boundary=include_boundary,
                include_vertices=include_vertices,
            )
        )

    random_cfg = _as_mapping(cfg.get("composition_random", None), name="composition_random")
    if "n_points" not in random_cfg:
        raise ValueError("config: composition_random.n_points is required")
    seed_raw = random_cfg.get("seed", cfg.get("seed", None))
    seed = None if seed_raw is None else int(seed_raw)
    return _random_simplex_compositions(
        n_elements=n_elements,
        n_points=int(random_cfg["n_points"]),
        seed=seed,
        alpha=random_cfg.get("alpha", 1.0),
        include_vertices=bool(random_cfg.get("include_vertices", False)),
    )


def _subst_a(obj: Any, a: float) -> Any:
    if isinstance(obj, str):
        if obj.strip() == "{a}":
            return float(a)
        return obj
    if isinstance(obj, list):
        return [_subst_a(x, a) for x in obj]
    if isinstance(obj, dict):
        return {k: _subst_a(v, a) for k, v in obj.items()}
    return obj


def _default_random_reference_lattice(length_min: float, length_max: float) -> list[list[float]]:
    length_ref = 0.5 * (float(length_min) + float(length_max))
    return [
        [length_ref, 0.0, 0.0],
        [0.0, length_ref, 0.0],
        [0.0, 0.0, length_ref],
    ]


def _resolve_frac_coords_init(
    frac_coords_init: list[Any] | None,
    null_num_sites: int | None,
) -> list[list[Any]]:
    rows = list(frac_coords_init or [])
    if null_num_sites is None:
        if not rows:
            raise ValueError("config: either frac_coords_init or nulls.num_sites must be provided")
        return rows

    retained_rows = []
    for row in rows:
        if len(row) == 3:
            label = "mix"
            x, y, z = row
        elif len(row) == 4:
            label, x, y, z = row
        else:
            raise ValueError(f"Bad frac_coords_init row: {row}")
        label_lower = str(label).lower()
        if label_lower not in {"mix", "null"}:
            retained_rows.append([label, float(x), float(y), float(z)])

    null_frac = np.random.random((int(null_num_sites), 3))
    null_rows = [["null", float(p[0]), float(p[1]), float(p[2])] for p in null_frac]
    return retained_rows + null_rows


def load_inputs(infile: Path) -> RunInputs:
    cfg = _load_yaml(infile)
    main_dir = Path(__file__).resolve().parent

    outdir = cfg.get("outdir", "out")
    out_dir = Path(outdir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    n_trials = int(cfg["n_trials"])

    fix = bool(cfg.get("fix", True))
    selected_elements = list(cfg["selected_elements"])
    compositions = _load_compositions(cfg, selected_elements)

    init_cfg = _as_mapping(cfg.get("initialization", None), name="initialization")
    syn_cfg = _as_mapping(cfg.get("syntropization", None), name="syntropization")
    null_cfg = _as_mapping(cfg.get("nulls", None), name="nulls")
    escape_cfg = _as_mapping(cfg.get("escape", None), name="escape")
    population_cfg = _as_mapping(cfg.get("population", None), name="population")
    cell_bounds_cfg = _as_mapping(cfg.get("cell_bounds", None), name="cell_bounds")

    alpha = float(init_cfg.get("alpha", 1.0))
    Tsallis_q = float(syn_cfg.get("Tsallis_q", 2.0))
    tau_mult = float(syn_cfg.get("tau_mult", 0.9))
    tau_init = float(syn_cfg.get("tau_init", 1.0))
    stage_maxit = int(syn_cfg.get("stage_maxit", 100))
    first_stage_maxit = int(syn_cfg.get("first_stage_maxit", stage_maxit))
    pbfgs_tol = float(syn_cfg.get("pbfgs_tol", 1e-3))
    composition_penalty_multiplier = float(syn_cfg.get("composition_penalty_multiplier", 1.0))

    if stage_maxit <= 0:
        raise ValueError("config: syntropization.stage_maxit must be > 0")
    if first_stage_maxit <= 0:
        raise ValueError("config: syntropization.first_stage_maxit must be > 0")
    if pbfgs_tol <= 0.0:
        raise ValueError("config: syntropization.pbfgs_tol must be > 0")

    null_enabled = bool(null_cfg.get("enabled", False))
    null_default_for_mix = bool(null_cfg.get("default_for_mix", null_enabled))
    null_num_sites_raw = null_cfg.get("num_sites", None)
    null_num_sites = None if null_num_sites_raw is None else int(null_num_sites_raw)
    null_min_distance = float(null_cfg.get("min_distance", 1.5))
    null_initial_occupancy_sum_raw = null_cfg.get("initial_occupancy_sum", None)
    null_initial_occupancy_sum = None if null_initial_occupancy_sum_raw is None else float(null_initial_occupancy_sum_raw)
    null_normalization_eps = float(null_cfg.get("normalization_eps", 1e-8))
    null_threshold = float(null_cfg.get("threshold", 1e-3))
    if null_min_distance < 0.0:
        raise ValueError("config: nulls.min_distance must be >= 0")
    if null_threshold < 0.0:
        raise ValueError("config: nulls.threshold must be >= 0")

    cell_init_mode = str(init_cfg.get("cell_init_mode", "input_cell"))
    if cell_init_mode not in {"input_cell", "random"}:
        raise ValueError("config: initialization.cell_init_mode must be one of {'input_cell','random'}")
    random_cell_length_min = float(init_cfg.get("random_cell_length_min", 2.0))
    random_cell_length_max = float(init_cfg.get("random_cell_length_max", 5.0))
    random_cell_angle_min_deg = float(init_cfg.get("random_cell_angle_min_deg", 60.0))
    random_cell_angle_max_deg = float(init_cfg.get("random_cell_angle_max_deg", 120.0))
    if random_cell_length_min > random_cell_length_max:
        raise ValueError("config: initialization.random_cell_length_min must be <= random_cell_length_max")
    if random_cell_angle_min_deg > random_cell_angle_max_deg:
        raise ValueError("config: initialization.random_cell_angle_min_deg must be <= random_cell_angle_max_deg")
    cell_bound_length_min = float(cell_bounds_cfg.get("length_min", 0.5))
    cell_bound_length_max = float(cell_bounds_cfg.get("length_max", 100.0))
    cell_bound_angle_min_deg = float(cell_bounds_cfg.get("angle_min_deg", 20.0))
    cell_bound_angle_max_deg = float(cell_bounds_cfg.get("angle_max_deg", 160.0))
    if cell_bound_length_min <= 0.0:
        raise ValueError("config: cell_bounds.length_min must be > 0")
    if cell_bound_length_min > cell_bound_length_max:
        raise ValueError("config: cell_bounds.length_min must be <= cell_bounds.length_max")
    if cell_bound_angle_min_deg <= 0.0 or cell_bound_angle_max_deg >= 180.0:
        raise ValueError("config: cell_bounds angles must be between 0 and 180 degrees")
    if cell_bound_angle_min_deg > cell_bound_angle_max_deg:
        raise ValueError("config: cell_bounds.angle_min_deg must be <= cell_bounds.angle_max_deg")
    if cell_init_mode == "random":
        if random_cell_length_min < cell_bound_length_min or random_cell_length_max > cell_bound_length_max:
            raise ValueError(
                "config: initialization random_cell_length_min/max must lie within "
                "cell_bounds.length_min/max"
            )
        if random_cell_angle_min_deg < cell_bound_angle_min_deg or random_cell_angle_max_deg > cell_bound_angle_max_deg:
            raise ValueError(
                "config: initialization random_cell_angle_min_deg/max_deg must lie within "
                "cell_bounds.angle_min_deg/max_deg"
            )

    escape_enabled = bool(escape_cfg.get("enabled", True))
    escape_trials = int(escape_cfg.get("trials", 3))
    escape_pos_sigma = float(escape_cfg.get("position_sigma", 0.3))
    escape_legacy_occ_sigma = escape_cfg.get("occupancy_sigma", 0.5)
    escape_empty_prob_sigma = float(escape_cfg.get("empty_probability_sigma", escape_legacy_occ_sigma))
    escape_species_prob_sigma = float(escape_cfg.get("species_probability_sigma", escape_legacy_occ_sigma))
    escape_cell_length_sigma = float(escape_cfg.get("cell_length_sigma", 0.0))
    escape_cell_angle_sigma_deg = float(escape_cfg.get("cell_angle_sigma_deg", 0.0))
    escape_trial_maxit = int(escape_cfg.get("trial_maxit", 50))
    escape_stagnation_window = int(escape_cfg.get("stagnation_window", 3))
    escape_min_shock_interval = int(escape_cfg.get("min_shock_interval", 1))
    escape_max_attempts_per_stage = int(escape_cfg.get("max_attempts_per_stage", 3))
    escape_accept_tol = float(escape_cfg.get("accept_tol", 1.0e-3))
    if escape_empty_prob_sigma < 0.0:
        raise ValueError("config: escape.empty_probability_sigma must be >= 0")
    if escape_species_prob_sigma < 0.0:
        raise ValueError("config: escape.species_probability_sigma must be >= 0")
    if escape_stagnation_window < 2:
        raise ValueError("config: escape.stagnation_window must be >= 2")
    if escape_min_shock_interval < 1:
        raise ValueError("config: escape.min_shock_interval must be >= 1")
    if escape_max_attempts_per_stage < 0:
        raise ValueError("config: escape.max_attempts_per_stage must be >= 0")
    if escape_accept_tol < 0.0:
        raise ValueError("config: escape.accept_tol must be >= 0")

    population_enabled = bool(population_cfg.get("enabled", False))
    population_max_generations = int(population_cfg.get("max_generations", 0))
    population_initial_runs = int(population_cfg.get("initial_runs", n_trials))
    population_parents_per_generation = int(population_cfg.get("parents_per_generation", max(1, min(n_trials, 4))))
    population_children_per_parent = int(population_cfg.get("children_per_parent", 1))
    population_keep_archive = int(population_cfg.get("keep_archive", max(n_trials, population_parents_per_generation)))
    population_memory_strength = float(population_cfg.get("memory_strength", 0.75))
    population_selection_metric = str(population_cfg.get("selection_metric", "formation_energy")).strip().lower()
    population_initial_lattice_vectors = _optional_lattice_vectors(
        population_cfg.get("initial_lattice_vectors", None),
        label="population.initial_lattice_vectors",
    )
    population_seed_parent_xsf = _resolve_optional_path(
        population_cfg.get("seed_parent_xsf", None),
        infile=infile,
        label="population.seed_parent_xsf",
    )
    population_seed_parent_symbols = None
    if population_max_generations < 0:
        raise ValueError("config: population.max_generations must be >= 0")
    if population_initial_runs <= 0:
        raise ValueError("config: population.initial_runs must be > 0")
    if population_parents_per_generation <= 0:
        raise ValueError("config: population.parents_per_generation must be > 0")
    if population_children_per_parent <= 0:
        raise ValueError("config: population.children_per_parent must be > 0")
    if population_keep_archive <= 0:
        raise ValueError("config: population.keep_archive must be > 0")
    if not (0.0 <= population_memory_strength <= 1.0):
        raise ValueError("config: population.memory_strength must be in [0, 1]")
    if population_selection_metric not in {"formation_energy", "energy_above_hull", "total_energy"}:
        raise ValueError(
            "config: population.selection_metric must be one of "
            "{'formation_energy','energy_above_hull','total_energy'}"
        )
    if population_seed_parent_xsf is not None and not population_enabled:
        raise ValueError("config: population.seed_parent_xsf requires population.enabled=true")

    abort_energy_above_hull_raw = syn_cfg.get("abort_energy_above_hull", None)
    abort_energy_above_hull = None if abort_energy_above_hull_raw is None else float(abort_energy_above_hull_raw)

    supercell_list = cfg["supercell_size"]
    if len(supercell_list) != 3:
        raise ValueError("config: supercell_size must have length 3")
    supercell_size = (int(supercell_list[0]), int(supercell_list[1]), int(supercell_list[2]))

    frac_coords_init_raw = cfg.get("frac_coords_init", None)
    lattice_vectors_raw = cfg.get("lattice_vectors", None)
    if population_seed_parent_xsf is not None:
        seed_lattice, seed_symbols, seed_frac = _read_xsf_structure(population_seed_parent_xsf)
        missing_seed_elements = sorted(set(seed_symbols) - set(selected_elements))
        if missing_seed_elements:
            raise ValueError(
                "config: population.seed_parent_xsf contains elements not in selected_elements: "
                f"{missing_seed_elements}"
            )
        lattice_vectors_raw = seed_lattice
        frac_coords_init_raw = [["null", float(x), float(y), float(z)] for x, y, z in seed_frac]
        null_num_sites = None
        population_seed_parent_symbols = list(seed_symbols)

    if lattice_vectors_raw is None:
        if population_initial_lattice_vectors is not None:
            lattice_vectors_raw = population_initial_lattice_vectors
        elif cell_init_mode != "random":
            raise ValueError(
                "config: lattice_vectors is required only for initialization.cell_init_mode='input_cell' "
                "when population.initial_lattice_vectors is not provided"
            )
        else:
            lattice_vectors_raw = _default_random_reference_lattice(random_cell_length_min, random_cell_length_max)
    elif "lattice_a" in cfg:
        a = float(cfg["lattice_a"])
        lattice_vectors_raw = _subst_a(lattice_vectors_raw, a)
    lattice_vectors = [[float(x) for x in row] for row in lattice_vectors_raw]
    frac_coords_init = _resolve_frac_coords_init(
        frac_coords_init=frac_coords_init_raw,
        null_num_sites=null_num_sites,
    )

    hull_db_filename = cfg.get("hull_db_filename", "hull_db.json")
    db_path = (out_dir / hull_db_filename).resolve()

    ref_override = cfg.get("reference_json_path", None)
    if ref_override is None:
        reference_json_path = (main_dir / "data" / "reference_energies.json").resolve()
    else:
        p = Path(ref_override).expanduser()
        reference_json_path = (p if p.is_absolute() else (infile.parent / p)).resolve()

    if not reference_json_path.exists():
        raise FileNotFoundError(f"Reference energies JSON not found: {reference_json_path}")
    with reference_json_path.open("r", encoding="utf-8") as f:
        reference_energies = json.load(f)
    missing_reference_elements = [el for el in selected_elements if el not in reference_energies]
    if missing_reference_elements:
        raise ValueError(
            "Reference energies JSON is missing selected elements: "
            f"{missing_reference_elements}. Path: {reference_json_path}"
        )

    endmembers = _merge_endmembers(selected_elements, cfg.get("endmembers", None))

    seed_raw = cfg.get("seed", None)
    seed = None if seed_raw is None else int(seed_raw)

    return RunInputs(
        out_dir=out_dir,
        n_trials=n_trials,
        compositions=compositions,
        fix=fix,
        alpha=alpha,
        Tsallis_q=Tsallis_q,
        tau_init=tau_init,
        tau_mult=tau_mult,
        first_stage_maxit=first_stage_maxit,
        stage_maxit=stage_maxit,
        pbfgs_tol=pbfgs_tol,
        composition_penalty_multiplier=composition_penalty_multiplier,
        selected_elements=selected_elements,
        lattice_vectors=lattice_vectors,
        frac_coords_init=frac_coords_init,
        supercell_size=supercell_size,
        db_path=db_path,
        reference_json_path=reference_json_path,
        endmembers=endmembers,
        null_default_for_mix=null_default_for_mix,
        null_num_sites=null_num_sites,
        null_min_distance=null_min_distance,
        null_initial_occupancy_sum=null_initial_occupancy_sum,
        null_normalization_eps=null_normalization_eps,
        null_threshold=null_threshold,
        cell_init_mode=cell_init_mode,
        random_cell_length_min=random_cell_length_min,
        random_cell_length_max=random_cell_length_max,
        random_cell_angle_min_deg=random_cell_angle_min_deg,
        random_cell_angle_max_deg=random_cell_angle_max_deg,
        cell_bound_length_min=cell_bound_length_min,
        cell_bound_length_max=cell_bound_length_max,
        cell_bound_angle_min_deg=cell_bound_angle_min_deg,
        cell_bound_angle_max_deg=cell_bound_angle_max_deg,
        escape_enabled=escape_enabled,
        escape_trials=escape_trials,
        escape_pos_sigma=escape_pos_sigma,
        escape_empty_prob_sigma=escape_empty_prob_sigma,
        escape_species_prob_sigma=escape_species_prob_sigma,
        escape_cell_length_sigma=escape_cell_length_sigma,
        escape_cell_angle_sigma_deg=escape_cell_angle_sigma_deg,
        escape_trial_maxit=escape_trial_maxit,
        escape_stagnation_window=escape_stagnation_window,
        escape_min_shock_interval=escape_min_shock_interval,
        escape_max_attempts_per_stage=escape_max_attempts_per_stage,
        escape_accept_tol=escape_accept_tol,
        population_enabled=population_enabled,
        population_max_generations=population_max_generations,
        population_initial_runs=population_initial_runs,
        population_parents_per_generation=population_parents_per_generation,
        population_children_per_parent=population_children_per_parent,
        population_keep_archive=population_keep_archive,
        population_memory_strength=population_memory_strength,
        population_selection_metric=population_selection_metric,
        population_initial_lattice_vectors=population_initial_lattice_vectors,
        population_seed_parent_xsf=population_seed_parent_xsf,
        population_seed_parent_symbols=population_seed_parent_symbols,
        abort_energy_above_hull=abort_energy_above_hull,
        seed=seed,
    )


def _score_result(result: OptimizationResult, metric: str) -> float:
    if result.aborted:
        return float("inf")
    if metric == "energy_above_hull":
        value = result.final_energy_above_hull_ev_per_atom
    elif metric == "total_energy":
        value = result.final_energy_total_ev
    else:
        value = result.final_formation_energy_ev_per_atom
    value = float(value)
    return value if np.isfinite(value) else float("inf")


def _population_entry_sort_key(entry: dict[str, Any]) -> tuple[float, int, int]:
    return (float(entry["score"]), int(entry["generation"]), int(entry["run_index"]))


def _seed_parent_entry_from_atoms(
    *,
    atoms_base,
    seed_symbols: list[str],
    selected_elements: list[str],
    seed_xsf_path: Path,
    composition_tag: str,
) -> dict[str, Any]:
    site_is_mix = np.asarray(atoms_base.arrays["site_is_mix"], dtype=bool)
    mix_indices = np.nonzero(site_is_mix)[0]
    if len(mix_indices) != len(seed_symbols):
        raise ValueError(
            "seed_parent_xsf requires one variable site per XSF atom. "
            f"Got {len(mix_indices)} variable sites and {len(seed_symbols)} XSF atoms."
        )

    element_to_col = {el: i for i, el in enumerate(selected_elements)}
    X_final = np.zeros((len(seed_symbols), len(selected_elements)), dtype=float)
    for i, symbol in enumerate(seed_symbols):
        if symbol not in element_to_col:
            raise ValueError(f"seed_parent_xsf symbol {symbol!r} is not in selected_elements={selected_elements}")
        X_final[i, element_to_col[symbol]] = 1.0

    counts = X_final.sum(axis=0)
    composition = counts / max(float(counts.sum()), 1.0)
    frac_coords = np.asarray(atoms_base.get_scaled_positions(), dtype=float)
    lattice_vectors = np.asarray(atoms_base.cell.array, dtype=float)
    db_label = " ".join(f"{el}{composition[i]:.3f}" for i, el in enumerate(selected_elements))

    result = OptimizationResult(
        history={},
        all_history={
            "schema_version": 3,
            "seed_parent_xsf": str(seed_xsf_path),
            "population": {
                "enabled": True,
                "generation": 0,
                "run_index": -1,
                "parent_run_index": None,
                "child_index": None,
                "selection_metric": "seed_parent",
            },
            "final_symbols": list(seed_symbols),
            "final_site_occupancies_mix": np.ones(len(seed_symbols), dtype=float).tolist(),
            "final_selected_composition": composition.tolist(),
            "final_energy_evaluations": 0,
        },
        anneal_boundaries=[],
        final_energy_total_ev=0.0,
        final_formation_energy_ev_per_atom=0.0,
        final_energy_above_hull_ev_per_atom=0.0,
        selected_elements=list(selected_elements),
        X_final=X_final,
        site_occupancies=np.ones(len(seed_symbols), dtype=float),
        composition_selected=composition,
        frac_coords=frac_coords,
        lattice_vectors=lattice_vectors,
        symbols=list(seed_symbols),
        db_label=db_label,
        energy_evaluations=0,
        cum_global_iter=0,
        aborted=False,
        abort_reason="",
    )

    return {
        "result": result,
        "score": 0.0,
        "generation": 0,
        "run_index": -1,
        "parent_run_index": None,
        "child_index": None,
        "seed_parent": True,
        "composition_tag": composition_tag,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=str, required=True, help="YAML config file")
    args = parser.parse_args()

    infile = Path(args.infile).expanduser().resolve()
    if not infile.exists():
        raise FileNotFoundError(f"infile not found: {infile}")

    inp = load_inputs(infile)

    if inp.seed is not None:
        np.random.seed(inp.seed)
        torch.manual_seed(inp.seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(inp.seed)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log = setup_logging(inp.out_dir)
    log.info(f"Output dir: {inp.out_dir}")
    log.info(f"Device: {device}")
    log.info(f"Infile: {infile}")
    log.info(f"Reference energies: {inp.reference_json_path}")
    log.info(f"Target compositions: {len(inp.compositions)}")
    if inp.seed is not None:
        log.info(f"Seed: {inp.seed}")
    histories_dir = inp.out_dir / "histories"
    histories_dir.mkdir(parents=True, exist_ok=True)
    progress_status_path = inp.out_dir / "logs" / "progress_status.log"
    progress_status_path.write_text(
        "High-level EAT progress status\n"
        "Energy evaluations count calls to objfun.get_E_F_S(), i.e. FairChem get_potential_energy().\n"
        "Includes main syntropization, shock/escape trials, stage-end checks, and final relaxation.\n"
        "\n",
        encoding="utf-8",
    )
    log.info(f"High-level progress status log: {progress_status_path}")

    rng_state = save_rng_state()
    predictor = pretrained_mlip.get_predict_unit("uma-s-1")
    calc = FAIRChemCalculator(predictor, task_name="omat")
    load_rng_state(rng_state)

    inp.db_path.parent.mkdir(parents=True, exist_ok=True)
    if not inp.db_path.exists():
        inp.db_path.write_text("{}", encoding="utf-8")

    atoms_base = build_base_atoms(
        inp.frac_coords_init,
        inp.lattice_vectors,
        inp.supercell_size,
        null_default_for_mix=inp.null_default_for_mix,
    )

    null_site_count = int(np.asarray(atoms_base.arrays["site_is_null"], dtype=int).sum())
    mix_site_count = int(np.asarray(atoms_base.arrays["site_is_mix"], dtype=int).sum())
    log.info(f"Variable sites: {mix_site_count}")
    log.info(f"Null sites: {null_site_count}")
    log.info(f"Null minimum distance: {inp.null_min_distance}")
    if inp.null_initial_occupancy_sum is not None:
        log.info(f"Initial null occupancy sum target: {inp.null_initial_occupancy_sum}")
    log.info(f"Composition penalty multiplier: {inp.composition_penalty_multiplier}")
    log.info(
        "Absolute cell bounds: "
        f"length=[{inp.cell_bound_length_min}, {inp.cell_bound_length_max}] Angstrom, "
        f"angle=[{inp.cell_bound_angle_min_deg}, {inp.cell_bound_angle_max_deg}] deg"
    )
    if inp.abort_energy_above_hull is not None:
        log.info(f"Abort stage if energy above hull exceeds: {inp.abort_energy_above_hull} eV/atom")
    log.info(f"First syntropization stage max PBFGS iterations: {inp.first_stage_maxit}")
    log.info(f"Syntropization stage max PBFGS iterations: {inp.stage_maxit}")
    log.info(f"PBFGS convergence tolerance: {inp.pbfgs_tol}")
    log.info(f"Null occupancy threshold: {inp.null_threshold}")
    log.info(
        "Shock trigger: objective_no_entropy stagnation over "
        f"{inp.escape_stagnation_window} recent stage iterations"
    )
    log.info(f"Shock maximum frequency: once every {inp.escape_min_shock_interval} annealing stages")
    log.info(
        "Shock within-stage repeats: "
        f"max_attempts_per_stage={inp.escape_max_attempts_per_stage}, "
        f"accept_tol={inp.escape_accept_tol} on objective_no_entropy"
    )
    log.info("Escape null reseeding: random, always enabled")
    log.info(
        "Shock perturbation sigmas: "
        f"p_empty={inp.escape_empty_prob_sigma}, "
        f"p_species={inp.escape_species_prob_sigma}"
    )
    log.info("Shock acceptance: accept the lowest objective_no_entropy trial only if it improves enough")
    if inp.population_enabled:
        log.info(
            "Population basin hopping enabled: "
            f"initial_runs={inp.population_initial_runs}, "
            f"max_generations={inp.population_max_generations}, "
            f"parents/generation={inp.population_parents_per_generation}, "
            f"children/parent={inp.population_children_per_parent}, "
            f"archive={inp.population_keep_archive}, "
            f"metric={inp.population_selection_metric}, "
            f"memory={inp.population_memory_strength}"
        )
        if inp.population_initial_lattice_vectors is not None:
            log.info("Population generation-0 initial lattice override enabled.")
        if inp.population_seed_parent_xsf is not None:
            log.info(f"Population seed parent XSF: {inp.population_seed_parent_xsf}")

    DB_JSON_PATH = str(inp.db_path)
    REFERENCE_JSON_PATH = str(inp.reference_json_path)

    ensure_endmembers_in_db_from_yaml(inp.db_path, inp.selected_elements, inp.endmembers)

    cum_iter = 0
    cum_energy_evals = 0
    best_formation_energy = float("inf")
    best_energy_eval = 0
    best_run_label = "none"
    for composition in inp.compositions:
        if len(composition) != len(inp.selected_elements):
            raise ValueError(
                f"Composition length {len(composition)} != number of elements {len(inp.selected_elements)} "
                f"for composition={composition} elements={inp.selected_elements}"
            )

        hullcalc = HullEnergyCalculator(DB_JSON_PATH, inp.selected_elements)
        mu = torch.tensor(composition, device=device)
        composition_tag = _composition_tag(inp.selected_elements, composition)
        run_counter = 0

        def run_one(
            *,
            generation: int,
            run_index: int,
            parent_entry: dict[str, Any] | None = None,
            child_index: int | None = None,
        ) -> dict[str, Any]:
            nonlocal cum_iter, cum_energy_evals, best_formation_energy, best_energy_eval, best_run_label
            parent_result = None if parent_entry is None else parent_entry["result"]
            parent_run_index = None if parent_entry is None else int(parent_entry["run_index"])
            run_role = "initial" if parent_entry is None else "child"
            child_text = "none" if child_index is None else str(child_index)
            parent_text = "none" if parent_run_index is None else str(parent_run_index)
            _append_progress_status(
                progress_status_path,
                (
                    f"START target={composition_tag} gen={generation} run={run_index} "
                    f"role={run_role} parent={parent_text} child={child_text} | "
                    f"current_evals={cum_energy_evals} | "
                    f"best_E_form={_fmt_energy(best_formation_energy)} eV/atom "
                    f"at_eval={best_energy_eval} best_run={best_run_label}"
                ),
            )

            syn = syntropizer(
                selected_elements=inp.selected_elements,
                atoms_base=atoms_base,
                calc=calc,
                log=log,
                mu=mu,
                fix=inp.fix,
                alpha=inp.alpha,
                Tsallis_q=inp.Tsallis_q,
                tau_init=inp.tau_init,
                tau_mult=inp.tau_mult,
                first_stage_maxit=inp.first_stage_maxit,
                stage_maxit=inp.stage_maxit,
                pbfgs_tol=inp.pbfgs_tol,
                DB_JSON_PATH=DB_JSON_PATH,
                REFERENCE_JSON_PATH=REFERENCE_JSON_PATH,
                null_min_distance=inp.null_min_distance,
                null_normalization_eps=inp.null_normalization_eps,
                null_threshold=inp.null_threshold,
                initial_occupancy_sum=inp.null_initial_occupancy_sum,
                cell_init_mode=inp.cell_init_mode,
                random_cell_length_min=inp.random_cell_length_min,
                random_cell_length_max=inp.random_cell_length_max,
                random_cell_angle_min_deg=inp.random_cell_angle_min_deg,
                random_cell_angle_max_deg=inp.random_cell_angle_max_deg,
                cell_bound_length_min=inp.cell_bound_length_min,
                cell_bound_length_max=inp.cell_bound_length_max,
                cell_bound_angle_min_deg=inp.cell_bound_angle_min_deg,
                cell_bound_angle_max_deg=inp.cell_bound_angle_max_deg,
                composition_penalty_multiplier=inp.composition_penalty_multiplier,
                escape_enabled=inp.escape_enabled,
                escape_trials=inp.escape_trials,
                escape_pos_sigma=inp.escape_pos_sigma,
                escape_empty_prob_sigma=inp.escape_empty_prob_sigma,
                escape_species_prob_sigma=inp.escape_species_prob_sigma,
                escape_cell_length_sigma=inp.escape_cell_length_sigma,
                escape_cell_angle_sigma_deg=inp.escape_cell_angle_sigma_deg,
                escape_trial_maxit=inp.escape_trial_maxit,
                escape_stagnation_window=inp.escape_stagnation_window,
                escape_min_shock_interval=inp.escape_min_shock_interval,
                escape_max_attempts_per_stage=inp.escape_max_attempts_per_stage,
                escape_accept_tol=inp.escape_accept_tol,
                restart_parent_result=parent_result,
                restart_memory_strength=inp.population_memory_strength,
                initial_lattice_vectors=inp.population_initial_lattice_vectors,
                abort_energy_above_hull=inp.abort_energy_above_hull,
            )

            result: OptimizationResult = syn.run()
            result.all_history["population"] = {
                "enabled": bool(inp.population_enabled),
                "generation": int(generation),
                "run_index": int(run_index),
                "parent_run_index": parent_run_index,
                "child_index": None if child_index is None else int(child_index),
                "selection_metric": inp.population_selection_metric,
            }

            label = _db_label_from_result(result)
            fractions = np.asarray(result.composition_selected, dtype=float).tolist()
            syn_iter = result.cum_global_iter
            cum_iter += syn_iter
            run_energy_evals = int(getattr(result, "energy_evaluations", 0))
            cum_energy_evals += run_energy_evals
            score = _score_result(result, inp.population_selection_metric)
            result.all_history["population"]["energy_evaluations_run"] = int(run_energy_evals)
            result.all_history["population"]["energy_evaluations_cumulative"] = int(cum_energy_evals)
            if (not result.aborted) and np.isfinite(float(result.final_formation_energy_ev_per_atom)):
                if float(result.final_formation_energy_ev_per_atom) < best_formation_energy:
                    best_formation_energy = float(result.final_formation_energy_ev_per_atom)
                    best_energy_eval = int(cum_energy_evals)
                    best_run_label = f"target={composition_tag},gen={generation},run={run_index}"

            _append_progress_status(
                progress_status_path,
                (
                    f"DONE  target={composition_tag} gen={generation} run={run_index} "
                    f"role={run_role} parent={parent_text} child={child_text} | "
                    f"status={'ABORTED' if result.aborted else 'OK'} | "
                    f"run_evals={run_energy_evals} current_evals={cum_energy_evals} | "
                    f"E_form={_fmt_energy(result.final_formation_energy_ev_per_atom)} eV/atom | "
                    f"E_above_hull={_fmt_energy(result.final_energy_above_hull_ev_per_atom)} eV/atom | "
                    f"current_best_E_form={_fmt_energy(best_formation_energy)} eV/atom "
                    f"best_at_eval={best_energy_eval} best_run={best_run_label}"
                ),
            )

            if result.aborted:
                log.warning(
                    f"Skipping DB update for target={composition_tag} gen={generation} run={run_index} "
                    f"because the run aborted: {result.abort_reason}"
                )
            else:
                record = PhaseRecord(
                    elements=result.selected_elements,
                    fractions=fractions,
                    formation_energy_per_atom=result.final_formation_energy_ev_per_atom,
                    total_energy_ev=result.final_energy_total_ev,
                    n_atoms=len(result.symbols),
                    symbols=list(result.symbols),
                    frac_coords=result.frac_coords.tolist(),
                    lattice_vectors=result.lattice_vectors.tolist(),
                    cell_pbc=[True, True, True],
                    supercell_size=list(inp.supercell_size),
                    metadata={
                        "source": "syntropizer_null_population" if inp.population_enabled else "syntropizer_null",
                        "site_occupancies": result.site_occupancies.tolist(),
                        "composition_selected": fractions,
                        "population_generation": int(generation),
                        "population_run_index": int(run_index),
                        "population_parent_run_index": parent_run_index,
                        "population_child_index": None if child_index is None else int(child_index),
                        "population_selection_metric": inp.population_selection_metric,
                        "population_score": float(score),
                        "energy_evaluations_run": int(run_energy_evals),
                        "energy_evaluations_cumulative": int(cum_energy_evals),
                    },
                )

                hullcalc.add_phase_to_db(label=label, record=record)
                log.info(f"Added/updated hull DB entry: {label}")

            if inp.population_enabled:
                parent_tag = "none" if parent_run_index is None else str(parent_run_index)
                child_tag = "none" if child_index is None else str(child_index)
                history_name = (
                    f"all_history_target={composition_tag}_gen={generation}_run={run_index}"
                    f"_parent={parent_tag}_child={child_tag}.json"
                )
            else:
                history_name = f"all_history_target={composition_tag}_i={run_index}.json"
            history_path = histories_dir / history_name
            with open(history_path, "w", encoding="utf-8") as f:
                json.dump(result.all_history, f)

            log.info(
                f"Completed target={composition_tag} gen={generation} run={run_index} "
                f"parent={parent_run_index} child={child_index} | "
                f"evals(run/cum)={run_energy_evals}/{cum_energy_evals} | "
                f"score({inp.population_selection_metric})={score:+.6f} | "
                f"E_form={result.final_formation_energy_ev_per_atom:+.6f} eV/atom | "
                f"E_above_hull={result.final_energy_above_hull_ev_per_atom:+.6f} eV/atom"
            )

            if (not result.aborted) and len(inp.selected_elements) <= 3:
                plot_path = plot_phase_diagram_from_db(
                    db_json_path=DB_JSON_PATH,
                    design_space=inp.selected_elements,
                    out_dir=inp.out_dir,
                    fname="phase_diagram.png",
                    annotate=False,
                    energy_evaluations=cum_energy_evals,
                )
                log.info(f"Saved phase diagram to {plot_path}")
            elif len(inp.selected_elements) > 3:
                log.info("Skipping phase diagram plot (only works for <= 3 elements)")

            return {
                "result": result,
                "score": float(score),
                "generation": int(generation),
                "run_index": int(run_index),
                "parent_run_index": parent_run_index,
                "child_index": None if child_index is None else int(child_index),
            }

        if not inp.population_enabled:
            for i in range(inp.n_trials):
                run_one(generation=0, run_index=i)
            continue

        seed_parent_entry = None
        if inp.population_seed_parent_xsf is not None:
            if inp.population_seed_parent_symbols is None:
                raise RuntimeError("Internal error: seed parent symbols were not loaded")
            seed_parent_entry = _seed_parent_entry_from_atoms(
                atoms_base=atoms_base,
                seed_symbols=inp.population_seed_parent_symbols,
                selected_elements=inp.selected_elements,
                seed_xsf_path=inp.population_seed_parent_xsf,
                composition_tag=composition_tag,
            )
            population_pool = [seed_parent_entry]
            log.info(
                f"[population target={composition_tag}] using XSF seed as the only generation-0 parent "
                f"(atoms={len(inp.population_seed_parent_symbols)}). Skipping random generation-0 runs."
            )
        else:
            population_pool: list[dict[str, Any]] = []
            for _ in range(inp.population_initial_runs):
                entry = run_one(generation=0, run_index=run_counter)
                run_counter += 1
                population_pool.append(entry)

            population_pool = sorted(population_pool, key=_population_entry_sort_key)[: inp.population_keep_archive]
            if population_pool:
                best = population_pool[0]
                log.info(
                    f"[population target={composition_tag}] generation 0 best run={best['run_index']} "
                    f"score={best['score']:+.6f}"
                )

        for generation in range(1, inp.population_max_generations + 1):
            if seed_parent_entry is not None:
                parents = [seed_parent_entry]
            else:
                selectable_pool = [entry for entry in population_pool if np.isfinite(float(entry["score"]))]
                parents = selectable_pool[: inp.population_parents_per_generation]
            if not parents:
                log.warning(f"[population target={composition_tag}] no valid parents remain; stopping.")
                break

            log.info("")
            log.info("=" * 80)
            log.info(
                f"POPULATION GENERATION {generation}/{inp.population_max_generations} | "
                f"parents={len(parents)} children/parent={inp.population_children_per_parent}"
            )
            log.info("=" * 80)

            children = []
            for parent_entry in parents:
                for child_idx in range(inp.population_children_per_parent):
                    entry = run_one(
                        generation=generation,
                        run_index=run_counter,
                        parent_entry=parent_entry,
                        child_index=child_idx,
                    )
                    run_counter += 1
                    children.append(entry)

            population_pool = sorted(population_pool + children, key=_population_entry_sort_key)[
                : inp.population_keep_archive
            ]
            best = population_pool[0]
            log.info(
                f"[population target={composition_tag}] generation {generation} best run={best['run_index']} "
                f"score={best['score']:+.6f} | archive_size={len(population_pool)}"
            )


if __name__ == "__main__":
    main()
