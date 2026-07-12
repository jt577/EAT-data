# plotting.py
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Tuple, Optional, Iterable

import numpy as np
import matplotlib.pyplot as plt


# -------------------------
# Helpers
# -------------------------
def _extract_x_for_element0(entry: Dict[str, Any], design_space: List[str]) -> Optional[float]:
    """
    Return x = fraction of design_space[0] for a DB entry/record, or None if not compatible.

    Expects:
      entry["elements"]  : list[str]
      entry["fractions"] : list[float]
    """
    elems = entry.get("elements", None)
    fracs = entry.get("fractions", None)
    if not isinstance(elems, list) or not isinstance(fracs, list) or len(elems) != len(fracs):
        return None

    ds = set(design_space)
    if not set(elems).issubset(ds):
        return None

    f = np.asarray(fracs, dtype=float)
    s = float(f.sum())
    if s <= 0:
        return None
    f = f / s

    frac_map = {e: float(fi) for e, fi in zip(elems, f)}
    return float(frac_map.get(design_space[0], 0.0))


def _iter_record_dicts(label: str, entry: Any) -> Iterable[Tuple[str, Dict[str, Any]]]:
    """
    Yield (label_variant, record_dict) pairs from a DB value.

    Supports your current schemas:
      A) entry is a single dict with "formation_energy_per_atom"
      B) entry is a dict with "entries": [ {...}, {...}, ... ]   <-- YOUR DB
      C) entry is a dict with "records": [ {...}, {...}, ... ]
      D) entry is a dict with "formation_energies_per_atom": [y1, y2, ...] (energies-only list)

    Notes:
      - If "entries"/"records" exists, we primarily plot those.
      - If those lists are empty or missing, we fall back to plotting the top-level entry.
    """
    if not isinstance(entry, dict):
        return

    # Prefer multi-entry payload if present (your DB uses "entries")
    for list_key in ("entries", "records"):
        rec_list = entry.get(list_key, None)
        if isinstance(rec_list, list) and len(rec_list) > 0:
            for j, rec in enumerate(rec_list):
                if isinstance(rec, dict):
                    yield f"{label}::{list_key}[{j}]", rec
            return  # do not also yield top-level unless list was empty

    # Energies-only list schema
    ys = entry.get("formation_energies_per_atom", None)
    if isinstance(ys, list) and len(ys) > 0:
        for j, y in enumerate(ys):
            rec = dict(entry)
            rec["formation_energy_per_atom"] = y
            yield f"{label}::Elist[{j}]", rec
        return

    # Fallback: single record
    yield label, entry


def _lower_hull(points: List[Tuple[float, float]]) -> List[Tuple[float, float]]:
    """
    Lower convex hull for 2D points (x, y) with x in [0,1].
    Returns hull vertices in increasing x.

    Assumes points are already "one per x" (or at least no pathological duplicates).
    """
    pts = sorted(points, key=lambda p: (p[0], p[1]))

    # De-duplicate exact duplicates
    uniq: List[Tuple[float, float]] = []
    for p in pts:
        if not uniq or (p[0] != uniq[-1][0] or p[1] != uniq[-1][1]):
            uniq.append(p)
    pts = uniq

    def cross(o, a, b):
        # z-component of (a-o) x (b-o)
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower: List[Tuple[float, float]] = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)

    return lower


def _xbin(x: float, tol: float) -> float:
    """
    Quantize x to a tolerance bin so near-equal compositions group together.
    Returns a representative x (the bin center on the tol grid).
    """
    if tol <= 0:
        raise ValueError("composition_tol must be > 0")
    k = int(np.round(x / tol))
    return float(k * tol)


# -------------------------
# Main plotter
# -------------------------
def plot_phase_diagram_from_db(
    db_json_path: str | Path,
    design_space: List[str],
    out_dir: str | Path,
    *,
    fname: str = "phase_diagram.png",
    annotate: bool = False,
    composition_tol: float = 1e-6,
    cum_iter: Optional[int] = None,
) -> Path:
    """
    Reads the DB, scatters ALL formation energies vs x=frac(design_space[0]),
    and overlays the convex hull constructed from the MINIMUM energy at each
    (tolerance-grouped) composition.

    Assumes BINARY design space: len(design_space)==2.
    """
    design_space = list(design_space)
    if len(design_space) != 2:
        raise ValueError("This plotter supports binary only (len(design_space)==2).")

    db_json_path = Path(db_json_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not db_json_path.exists():
        raise FileNotFoundError(f"DB json not found: {db_json_path}")

    db = json.loads(db_json_path.read_text(encoding="utf-8") or "{}")

    # All points (scatter)
    xs_all: List[float] = []
    ys_all: List[float] = []
    labels_all: List[str] = []

    # For hull: min energy per composition bin
    min_energy_by_x: Dict[float, float] = {}

    for label, entry in db.items():
        for lab2, rec in _iter_record_dicts(label, entry):
            y = rec.get("formation_energy_per_atom", None)
            if y is None:
                continue
            try:
                y = float(y)
            except Exception:
                continue

            x = _extract_x_for_element0(rec, design_space)
            if x is None:
                continue
            x = float(x)

            # Scatter: ALL points
            xs_all.append(x)
            ys_all.append(y)
            labels_all.append(lab2)

            # Hull: min per x-bin
            xb = _xbin(x, composition_tol)
            prev = min_energy_by_x.get(xb, None)
            if prev is None or y < prev:
                min_energy_by_x[xb] = y

    if len(xs_all) == 0:
        fig = plt.figure(figsize=(8, 4))
        plt.title("No compatible entries found in DB for design space")
        out_path = out_dir / fname
        fig.savefig(out_path, dpi=200)
        plt.close(fig)
        return out_path

    xs_np = np.asarray(xs_all, dtype=float)
    ys_np = np.asarray(ys_all, dtype=float)

    # Hull input: ensure ONE y per x (already min per bin)
    hull_input = sorted([(xk, yk) for xk, yk in min_energy_by_x.items()], key=lambda p: p[0])

    # Need at least 2 points to draw a line, but hull algo likes >= 3
    if len(hull_input) >= 3:
        hull_pts = _lower_hull(hull_input)
    else:
        hull_pts = hull_input

    hull_x = np.array([p[0] for p in hull_pts], dtype=float) if hull_pts else np.array([])
    hull_y = np.array([p[1] for p in hull_pts], dtype=float) if hull_pts else np.array([])

    fig = plt.figure(figsize=(9, 5))
    ax = plt.gca()

    ax.scatter(
        xs_np, ys_np,
        s=50, alpha=0.35, edgecolors="k", linewidths=0.3,
        label="All structures"
    )

    if len(hull_x) >= 2:
        ax.plot(hull_x, hull_y, linewidth=1.8, color="red", label="Convex hull (min per composition)")
        ax.scatter(hull_x, hull_y, s=40, color="red", zorder=5)
    elif len(hull_x) == 1:
        ax.scatter(hull_x, hull_y, s=40, color="red", zorder=5, label="Min per composition")

    if annotate:
        for x, y, lab in zip(xs_np, ys_np, labels_all):
            ax.text(x, y, lab, fontsize=7, alpha=0.8)

    ax.set_xlabel(f"x = fraction({design_space[0]})")
    ax.set_ylabel("Formation energy (eV/atom)")
    ax.set_title(f"Binary phase diagram: {design_space[0]}–{design_space[1]}")
    ax.grid(True, alpha=0.25)
    ax.set_xlim(-0.02, 1.02)
    ax.legend()

    if cum_iter is not None:
        ax.text(
            0.98, 0.02,
            f"Cumulative iteration: {cum_iter}",
            horizontalalignment='right',
            verticalalignment='top',
            transform=ax.transAxes,
            fontsize=14,
            color='gray',
            alpha=1.0
        )


    out_path = out_dir / fname
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path
