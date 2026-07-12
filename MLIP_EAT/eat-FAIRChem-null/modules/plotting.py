# plotting.py
from __future__ import annotations

import json
import itertools
import math
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


def _extract_composition(entry: Dict[str, Any], design_space: List[str]) -> Optional[np.ndarray]:
    """
    Return a normalized composition vector in design_space order, or None if incompatible.
    """
    elems = entry.get("elements", None)
    fracs = entry.get("fractions", None)
    if not isinstance(elems, list) or not isinstance(fracs, list) or len(elems) != len(fracs):
        return None

    ds = set(design_space)
    if not set(elems).issubset(ds):
        return None

    f = np.asarray(fracs, dtype=float)
    if f.ndim != 1 or not np.all(np.isfinite(f)):
        return None
    s = float(f.sum())
    if s <= 0:
        return None
    f = f / s

    frac_map = {e: float(fi) for e, fi in zip(elems, f)}
    comp = np.asarray([frac_map.get(el, 0.0) for el in design_space], dtype=float)
    comp_sum = float(comp.sum())
    if comp_sum <= 0:
        return None
    return comp / comp_sum


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


def _composition_bin_key(comp: np.ndarray, tol: float) -> tuple[int, ...]:
    if tol <= 0:
        raise ValueError("composition_tol must be > 0")
    return tuple(int(v) for v in np.round(np.asarray(comp, dtype=float) / tol))


def _composition_from_bin_key(key: tuple[int, ...], tol: float) -> np.ndarray:
    comp = np.asarray(key, dtype=float) * float(tol)
    comp = np.clip(comp, 0.0, None)
    s = float(comp.sum())
    if s <= 0:
        return comp
    return comp / s


def _ternary_to_xy(compositions: np.ndarray) -> np.ndarray:
    comps = np.asarray(compositions, dtype=float)
    if comps.ndim == 1:
        comps = comps.reshape(1, -1)
    if comps.shape[1] != 3:
        raise ValueError("ternary coordinates require composition shape (N, 3)")
    h = math.sqrt(3.0) / 2.0
    x = comps[:, 1] + 0.5 * comps[:, 2]
    y = h * comps[:, 2]
    return np.column_stack([x, y])


def _draw_ternary_frame(ax, design_space: List[str]) -> None:
    h = math.sqrt(3.0) / 2.0
    tri = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.5, h], [0.0, 0.0]], dtype=float)
    ax.plot(tri[:, 0], tri[:, 1], color="black", linewidth=1.2)

    ticks = np.linspace(0.2, 0.8, 4)
    for t in ticks:
        # Constant fraction of design_space[0], [1], [2], respectively.
        p0 = _ternary_to_xy(np.asarray([[t, 1.0 - t, 0.0], [t, 0.0, 1.0 - t]]))
        p1 = _ternary_to_xy(np.asarray([[1.0 - t, t, 0.0], [0.0, t, 1.0 - t]]))
        p2 = _ternary_to_xy(np.asarray([[1.0 - t, 0.0, t], [0.0, 1.0 - t, t]]))
        for p in (p0, p1, p2):
            ax.plot(p[:, 0], p[:, 1], color="0.82", linewidth=0.5, zorder=0)

    ax.text(-0.045, -0.045, design_space[0], ha="right", va="top", fontsize=12)
    ax.text(1.045, -0.045, design_space[1], ha="left", va="top", fontsize=12)
    ax.text(0.5, h + 0.045, design_space[2], ha="center", va="bottom", fontsize=12)

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(-0.08, 1.08)
    ax.set_ylim(-0.08, h + 0.08)
    ax.axis("off")


def _lower_hull_triangles_bruteforce(
    points_3d: np.ndarray,
    *,
    max_points: int = 120,
    tol: float = 1e-9,
) -> np.ndarray:
    n = len(points_3d)
    if n < 4 or n > max_points:
        return np.empty((0, 3), dtype=int)

    triangles: list[tuple[int, int, int]] = []
    pts = np.asarray(points_3d, dtype=float)
    for i, j, k in itertools.combinations(range(n), 3):
        normal = np.cross(pts[j] - pts[i], pts[k] - pts[i])
        norm = float(np.linalg.norm(normal))
        if norm <= tol or abs(float(normal[2])) <= tol:
            continue
        normal = normal / norm
        if normal[2] > 0.0:
            normal = -normal
        signed = (pts - pts[i]) @ normal
        if np.all(signed <= tol):
            triangles.append(tuple(sorted((i, j, k))))

    if not triangles:
        return np.empty((0, 3), dtype=int)
    return np.asarray(sorted(set(triangles)), dtype=int)


def _lower_hull_triangles(compositions: np.ndarray, energies: np.ndarray) -> np.ndarray:
    comps = np.asarray(compositions, dtype=float)
    energies = np.asarray(energies, dtype=float)
    if len(comps) < 4:
        return np.empty((0, 3), dtype=int)

    points_3d = np.column_stack([comps[:, 0], comps[:, 1], energies])
    if np.linalg.matrix_rank(points_3d - points_3d.mean(axis=0)) < 3:
        return np.empty((0, 3), dtype=int)

    try:
        from scipy.spatial import ConvexHull

        hull = ConvexHull(points_3d)
        lower = []
        for simplex, equation in zip(hull.simplices, hull.equations):
            if float(equation[2]) < -1e-10:
                lower.append(tuple(int(i) for i in simplex))
        return np.asarray(lower, dtype=int) if lower else np.empty((0, 3), dtype=int)
    except Exception:
        return _lower_hull_triangles_bruteforce(points_3d)


def _plot_no_compatible_entries(out_dir: Path, fname: str, title: str) -> Path:
    fig = plt.figure(figsize=(8, 4))
    plt.title(title)
    out_path = out_dir / fname
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def _plot_binary_phase_diagram_from_db(
    db: Dict[str, Any],
    design_space: List[str],
    out_dir: Path,
    *,
    fname: str,
    annotate: bool,
    composition_tol: float,
    cum_iter: Optional[int],
    energy_evaluations: Optional[int],
) -> Path:
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

            xs_all.append(x)
            ys_all.append(y)
            labels_all.append(lab2)

            xb = _xbin(x, composition_tol)
            prev = min_energy_by_x.get(xb, None)
            if prev is None or y < prev:
                min_energy_by_x[xb] = y

    if len(xs_all) == 0:
        return _plot_no_compatible_entries(out_dir, fname, "No compatible entries found in DB for design space")

    xs_np = np.asarray(xs_all, dtype=float)
    ys_np = np.asarray(ys_all, dtype=float)

    hull_input = sorted([(xk, yk) for xk, yk in min_energy_by_x.items()], key=lambda p: p[0])
    hull_pts = _lower_hull(hull_input) if len(hull_input) >= 3 else hull_input
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
    ax.set_title(f"Binary phase diagram: {design_space[0]}-{design_space[1]}")
    ax.grid(True, alpha=0.25)
    ax.set_xlim(-0.02, 1.02)
    ax.legend()

    if energy_evaluations is not None:
        status_text = f"Energy evaluations: {energy_evaluations}"
    elif cum_iter is not None:
        status_text = f"Cumulative iteration: {cum_iter}"
    else:
        status_text = None

    if status_text is not None:
        ax.text(
            0.98, 0.02,
            status_text,
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
            fontsize=14,
            color="gray",
            alpha=1.0,
        )

    out_path = out_dir / fname
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path


def _plot_ternary_phase_diagram_from_db(
    db: Dict[str, Any],
    design_space: List[str],
    out_dir: Path,
    *,
    fname: str,
    annotate: bool,
    composition_tol: float,
    cum_iter: Optional[int],
    energy_evaluations: Optional[int],
) -> Path:
    comps_all: list[np.ndarray] = []
    ys_all: list[float] = []
    labels_all: list[str] = []
    min_by_comp: dict[tuple[int, ...], tuple[np.ndarray, float]] = {}

    for label, entry in db.items():
        for lab2, rec in _iter_record_dicts(label, entry):
            y = rec.get("formation_energy_per_atom", None)
            if y is None:
                continue
            try:
                y = float(y)
            except Exception:
                continue

            comp = _extract_composition(rec, design_space)
            if comp is None:
                continue

            comps_all.append(comp)
            ys_all.append(y)
            labels_all.append(lab2)

            key = _composition_bin_key(comp, composition_tol)
            prev = min_by_comp.get(key, None)
            if prev is None or y < prev[1]:
                min_by_comp[key] = (_composition_from_bin_key(key, composition_tol), y)

    if len(comps_all) == 0:
        return _plot_no_compatible_entries(out_dir, fname, "No compatible ternary entries found in DB for design space")

    comps_np = np.asarray(comps_all, dtype=float)
    ys_np = np.asarray(ys_all, dtype=float)
    xy_all = _ternary_to_xy(comps_np)

    hull_records = list(min_by_comp.values())
    hull_comps = np.asarray([rec[0] for rec in hull_records], dtype=float)
    hull_energies = np.asarray([rec[1] for rec in hull_records], dtype=float)
    lower_triangles = _lower_hull_triangles(hull_comps, hull_energies)

    fig = plt.figure(figsize=(8, 7))
    ax = plt.gca()
    _draw_ternary_frame(ax, design_space)

    if len(lower_triangles) > 0:
        stable_indices = sorted(set(int(i) for tri in lower_triangles for i in tri))
        for tri in lower_triangles:
            xy_tri = _ternary_to_xy(hull_comps[np.asarray(tri, dtype=int)])
            ax.fill(
                xy_tri[:, 0],
                xy_tri[:, 1],
                facecolor="red",
                edgecolor="red",
                linewidth=0.8,
                alpha=0.10,
                zorder=1,
            )
        stable_xy = _ternary_to_xy(hull_comps[np.asarray(stable_indices, dtype=int)])
        ax.scatter(
            stable_xy[:, 0],
            stable_xy[:, 1],
            s=85,
            facecolors="none",
            edgecolors="red",
            linewidths=1.2,
            label="Lower-hull vertices",
            zorder=4,
        )

    scatter = ax.scatter(
        xy_all[:, 0],
        xy_all[:, 1],
        c=ys_np,
        cmap="viridis",
        s=46,
        alpha=0.78,
        edgecolors="black",
        linewidths=0.25,
        label="All structures",
        zorder=3,
    )
    cbar = fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Formation energy (eV/atom)")

    if annotate:
        for (x, y), lab in zip(xy_all, labels_all):
            ax.text(x, y, lab, fontsize=6, alpha=0.75)

    ax.set_title(f"Ternary phase diagram: {design_space[0]}-{design_space[1]}-{design_space[2]}")
    if len(lower_triangles) > 0:
        ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.02), frameon=False, ncol=2)

    if energy_evaluations is not None:
        status_text = f"Energy evaluations: {energy_evaluations}"
    elif cum_iter is not None:
        status_text = f"Cumulative iteration: {cum_iter}"
    else:
        status_text = None

    if status_text is not None:
        ax.text(
            0.98,
            0.02,
            status_text,
            horizontalalignment="right",
            verticalalignment="bottom",
            transform=ax.transAxes,
            fontsize=12,
            color="gray",
            alpha=1.0,
        )

    out_path = out_dir / fname
    fig.tight_layout()
    fig.savefig(out_path, dpi=300)
    plt.close(fig)
    return out_path


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
    energy_evaluations: Optional[int] = None,
) -> Path:
    """
    Plot binary or ternary phase diagrams from the DB.

    Binary: formation energy vs fraction(design_space[0]), with lower hull line.
    Ternary: simplex scatter colored by formation energy, with projected lower
    convex-hull facets from the minimum energy at each composition bin.
    """
    design_space = list(design_space)
    if len(design_space) not in (2, 3):
        raise ValueError("This plotter supports binary or ternary systems only (len(design_space) in {2, 3}).")

    db_json_path = Path(db_json_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not db_json_path.exists():
        raise FileNotFoundError(f"DB json not found: {db_json_path}")

    db = json.loads(db_json_path.read_text(encoding="utf-8") or "{}")

    if len(design_space) == 2:
        return _plot_binary_phase_diagram_from_db(
            db,
            design_space,
            out_dir,
            fname=fname,
            annotate=annotate,
            composition_tol=composition_tol,
            cum_iter=cum_iter,
            energy_evaluations=energy_evaluations,
        )

    return _plot_ternary_phase_diagram_from_db(
        db,
        design_space,
        out_dir,
        fname=fname,
        annotate=annotate,
        composition_tol=composition_tol,
        cum_iter=cum_iter,
        energy_evaluations=energy_evaluations,
    )
