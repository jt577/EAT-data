#!/usr/bin/env python3
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import json
import numpy as np
import torch

import yaml  # pip install pyyaml

# FairChem
from fairchem.core import FAIRChemCalculator, pretrained_mlip
from modules.atoms_helpers import build_base_atoms
from modules.syntropizer import syntropizer, OptimizationResult
from modules.energy_above_hull import HullEnergyCalculator, PhaseRecord
from modules.logging_utils import setup_logging
from modules.plotting import plot_phase_diagram_from_db
from modules.seeding import save_rng_state, load_rng_state


def _db_label_from_result(result: OptimizationResult) -> str:
    if getattr(result, "db_label", None):
        return result.db_label
    comp = np.mean(result.X_final, axis=0)
    parts = [f"{el}{comp[i]:.3f}" for i, el in enumerate(result.selected_elements)]
    return " ".join(parts)


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



@dataclass(frozen=True)
class RunInputs:
    out_dir: Path
    n_trials: int
    compositions: list[list[float]]
    fix: bool
    selected_elements: list[str]
    lattice_vectors: list[list[float]]
    frac_coords_init: Any
    supercell_size: tuple[int, int, int]
    db_path: Path
    reference_json_path: Path
    endmembers: list[dict[str, Any]]   # NEW


def _load_yaml(path: Path) -> dict[str, Any]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"YAML root must be a mapping/object. Got: {type(data)}")
    return data


def _subst_a(obj: Any, a: float) -> Any:
    """Recursively replace the string '{a}' with the float a inside lists/dicts."""
    if isinstance(obj, str):
        if obj.strip() == "{a}":
            return float(a)
        return obj
    if isinstance(obj, list):
        return [_subst_a(x, a) for x in obj]
    if isinstance(obj, dict):
        return {k: _subst_a(v, a) for k, v in obj.items()}
    return obj


def load_inputs(infile: Path) -> RunInputs:
    cfg = _load_yaml(infile)

    # main.py directory (reference point for repo-relative data paths)
    main_dir = Path(__file__).resolve().parent

    outdir = cfg.get("outdir", "out")
    out_dir = Path(outdir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    n_trials = int(cfg["n_trials"])
    compositions = cfg["compositions"]
    if not isinstance(compositions, list) or not compositions:
        raise ValueError("config: compositions must be a non-empty list of lists")

    fix = bool(cfg.get("fix", True))
    selected_elements = list(cfg["selected_elements"])

    supercell_list = cfg["supercell_size"]
    if len(supercell_list) != 3:
        raise ValueError("config: supercell_size must have length 3")
    supercell_size = (int(supercell_list[0]), int(supercell_list[1]), int(supercell_list[2]))

    frac_coords_init = cfg["frac_coords_init"]
    lattice_vectors_raw = cfg["lattice_vectors"]

    # Optional "{a}" substitution
    if "lattice_a" in cfg:
        a = float(cfg["lattice_a"])
        lattice_vectors_raw = _subst_a(lattice_vectors_raw, a)

    lattice_vectors = [[float(x) for x in row] for row in lattice_vectors_raw]

    # DB path lives in out_dir by default, but name can be changed
    hull_db_filename = cfg.get("hull_db_filename", "hull_db.json")
    db_path = (out_dir / hull_db_filename).resolve()

    # Reference energies JSON is ALWAYS relative to main.py location (unless you *explicitly* override)
    # Your requirement: "data/reference_energies.json where the reference is where main.py code is."
    ref_override = cfg.get("reference_json_path", None)
    if ref_override is None:
        reference_json_path = (main_dir / "data" / "reference_energies.json").resolve()
    else:
        # If user provides one, interpret relative to infile location (common expectation),
        # or allow absolute.
        p = Path(ref_override).expanduser()
        reference_json_path = (p if p.is_absolute() else (infile.parent / p)).resolve()

    if not reference_json_path.exists():
        raise FileNotFoundError(f"Reference energies JSON not found: {reference_json_path}")

    endmembers = cfg.get("endmembers", None)
    if endmembers is None or not isinstance(endmembers, list) or len(endmembers) == 0:
        raise ValueError("config: endmembers must be a non-empty list (fractions + formation_energy_per_atom).")

    # validate endmembers
    for em in endmembers:
        if not isinstance(em, dict):
            raise ValueError(f"config: each endmember must be a dict, got {type(em)}")
        if "fractions" not in em or "formation_energy_per_atom" not in em:
            raise ValueError("config: each endmember needs 'fractions' and 'formation_energy_per_atom'")
        fr = em["fractions"]
        if not isinstance(fr, list) or len(fr) != len(selected_elements):
            raise ValueError(
                f"config: endmember fractions must have length {len(selected_elements)} (selected_elements order). "
                f"Got: {fr}"
            )
        # optional: sanity check for "pure" endpoints
        # (allows floats too)
        s = float(sum(fr))
        if abs(s - 1.0) > 1e-6:
            raise ValueError(f"config: endmember fractions must sum to 1. Got sum={s} fr={fr}")


    return RunInputs(
        out_dir=out_dir,
        n_trials=n_trials,
        compositions=[[float(v) for v in comp] for comp in compositions],
        fix=fix,
        selected_elements=selected_elements,
        lattice_vectors=lattice_vectors,
        frac_coords_init=frac_coords_init,
        supercell_size=supercell_size,
        db_path=db_path,
        reference_json_path=reference_json_path,
        endmembers=endmembers,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--infile", type=str, required=True, help="YAML config file")
    args = parser.parse_args()

    infile = Path(args.infile).expanduser().resolve()
    if not infile.exists():
        raise FileNotFoundError(f"infile not found: {infile}")

    inp = load_inputs(infile)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log = setup_logging(inp.out_dir)
    log.info(f"Output dir: {inp.out_dir}")
    log.info(f"Device: {device}")
    log.info(f"Infile: {infile}")
    log.info(f"Reference energies: {inp.reference_json_path}")

    # Save RNG state before FairChem touches seeds
    rng_state = save_rng_state()

    # -----------------------------
    # FAIRCHEM CALCULATOR
    # -----------------------------
    predictor = pretrained_mlip.get_predict_unit("uma-s-1")
    calc = FAIRChemCalculator(predictor, task_name="omat")

    # Restore RNG state so your random Dirichlet init isn’t “stuck”
    load_rng_state(rng_state)

    # -----------------------------
    # DB init
    # -----------------------------
    inp.db_path.parent.mkdir(parents=True, exist_ok=True)
    if not inp.db_path.exists():
        inp.db_path.write_text("{}", encoding="utf-8")

    # atoms template (same for all trials)
    atoms_base = build_base_atoms(inp.frac_coords_init, inp.lattice_vectors, inp.supercell_size)

    DB_JSON_PATH = str(inp.db_path)
    REFERENCE_JSON_PATH = str(inp.reference_json_path)

    # ensure endmembers exist BEFORE hull calculator is built
    ensure_endmembers_in_db_from_yaml(
        inp.db_path,
        inp.selected_elements,
        inp.endmembers,
    )

    # -----------------------------
    # OPTIMIZE
    # -----------------------------
    for composition in inp.compositions:
        if len(composition) != len(inp.selected_elements):
            raise ValueError(
                f"Composition length {len(composition)} != number of elements {len(inp.selected_elements)} "
                f"for composition={composition} elements={inp.selected_elements}"
            )

        # Construct hullcalc per composition loop (or move outside if you prefer; either works)
        hullcalc = HullEnergyCalculator(DB_JSON_PATH, inp.selected_elements)

        for i in range(inp.n_trials):
            mu = torch.tensor(composition, device=device)

            syn = syntropizer(
                selected_elements=inp.selected_elements,
                atoms_base=atoms_base,
                calc=calc,
                log=log,
                mu=mu,
                fix=inp.fix,
                DB_JSON_PATH=DB_JSON_PATH,
                REFERENCE_JSON_PATH=REFERENCE_JSON_PATH,
            )

            result: OptimizationResult = syn.run()

            label = _db_label_from_result(result)
            fractions = np.mean(result.X_final, axis=0).tolist()

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
                metadata={"source": "syntropizer"},
            )

            hullcalc.add_phase_to_db(label=label, record=record)
            log.info(f"✅ Added/updated hull DB entry: {label}")

            # Save artifacts (optionally include trial index to avoid overwriting)
            np.savez(
                inp.out_dir / f"result_trial_{i:04d}.npz",
                X_final=result.X_final,
                frac_coords=result.frac_coords,
                lattice_vectors=result.lattice_vectors,
                energy_total_ev=result.final_energy_total_ev,
                formation_energy_ev_per_atom=result.final_formation_energy_ev_per_atom,
                energy_above_hull_ev_per_atom=result.final_energy_above_hull_ev_per_atom,
            )

            if len(inp.selected_elements) <= 2:
                plot_path = plot_phase_diagram_from_db(
                    db_json_path=DB_JSON_PATH,
                    design_space=inp.selected_elements,
                    out_dir=inp.out_dir,
                    fname="phase_diagram.png",
                    annotate=False,
                )
                log.info(f"📈 Saved phase diagram to {plot_path}")
            else:
                log.info("ℹ️ Skipping phase diagram plot (only works for <= 2 elements)")


if __name__ == "__main__":
    main()
