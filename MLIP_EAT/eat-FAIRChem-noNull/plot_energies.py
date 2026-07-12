import numpy as np
import matplotlib.pyplot as plt
import json

DB_PATH = "/data2/jt577/2025/December_2025/Phase_Diagram/fairchem_phase_diagram/out3/hull_db.json"
db = json.load(open(DB_PATH, "r"))

def lower_convex_hull_1d(x, y):
    """
    Return (x_unique, y_unique, hull_idx) for the LOWER convex hull of (x,y),
    assuming x is 1D composition. Uses monotone chain on sorted x.

    Important:
      - If multiple points have same x, we keep only the lowest y (for hull input).
      - Hull is computed on those unique-x min-y points.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    # sort by x then y
    order = np.lexsort((y, x))
    x = x[order]
    y = y[order]

    # If multiple points have same x, keep only the lowest y
    x_unique = []
    y_unique = []
    for xi in np.unique(x):
        mask = (x == xi)
        x_unique.append(xi)
        y_unique.append(y[mask].min())
    x = np.array(x_unique, dtype=float)
    y = np.array(y_unique, dtype=float)

    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])

    hull = []
    for i in range(len(x)):
        p = (x[i], y[i])
        while len(hull) >= 2:
            o = (x[hull[-2]], y[hull[-2]])
            a = (x[hull[-1]], y[hull[-1]])
            # For LOWER hull: pop if we would make a non-left turn (clockwise or collinear)
            if cross(o, a, p) <= 0:
                hull.pop()
            else:
                break
        hull.append(i)

    return x, y, np.array(hull, dtype=int)

def _get_frac_map(elements, fractions):
    f = np.asarray(fractions, dtype=float)
    s = float(f.sum())
    if s <= 0:
        return None
    f = f / s
    return {el: float(fi) for el, fi in zip(elements, f)}

def _iter_records_from_db_value(val):
    """
    Supports:
      - new schema: val has "entries": [record, record, ...]
      - old/single schema: val itself is a record with "formation_energy_per_atom"
    """
    if not isinstance(val, dict):
        return
    if isinstance(val.get("entries", None), list):
        for rec in val["entries"]:
            if isinstance(rec, dict):
                yield rec
    else:
        # single-record style
        yield val

# ---------------------------
# Collect ALL scatter points
# and min-energy-per-composition for hull
# ---------------------------
composition_tol = 1e-6  # your grid is exact (0, 1/16, ...), so this is fine
if composition_tol <= 0:
    raise ValueError("composition_tol must be > 0")
decimals = int(max(0, np.ceil(-np.log10(composition_tol))))

xs_all = []
ys_all = []

# map: xkey -> min energy at that x
min_energy_by_xkey = {}

for label, val in db.items():
    for rec in _iter_records_from_db_value(val):
        y = rec.get("formation_energy_per_atom", None)
        if y is None:
            continue
        try:
            y = float(y)
        except Exception:
            continue
        if not np.isfinite(y):
            continue

        elems = rec.get("elements", None)
        fracs = rec.get("fractions", None)
        if not isinstance(elems, list) or not isinstance(fracs, list) or len(elems) != len(fracs):
            continue

        frac_map = _get_frac_map(elems, fracs)
        if frac_map is None:
            continue

        # x = Cu fraction (change "Cu" to design_space[0] if you want generic)
        x = float(frac_map.get("Cu", 0.0))

        xs_all.append(x)
        ys_all.append(y)

        xkey = round(x, decimals)
        prev = min_energy_by_xkey.get(xkey, None)
        if prev is None or y < prev:
            min_energy_by_xkey[xkey] = y

xs_all = np.asarray(xs_all, dtype=float)
ys_all = np.asarray(ys_all, dtype=float)

if xs_all.size == 0:
    raise RuntimeError("No valid (x, energy) points found in DB. Check DB schema / keys.")

# ---------------------------
# Hull from min per composition
# ---------------------------
x_min = np.array(sorted(min_energy_by_xkey.keys()), dtype=float)
e_min = np.array([min_energy_by_xkey[xk] for xk in x_min], dtype=float)

x_u, emin_u, hull_idx = lower_convex_hull_1d(x_min, e_min)
xh = x_u[hull_idx]
Eh = emin_u[hull_idx]

# ---------------------------
# Plot (keep your style)
# ---------------------------
plt.figure(figsize=(8, 6))

# all points (faint)
plt.scatter(xs_all, ys_all, color="blue", alpha=0.08, s=18)

# hull (line through hull vertices)
plt.plot(xh, Eh, "-o", linewidth=2.5, markersize=6, label="Convex hull (lower envelope)")

plt.xlabel("Cu Concentration", fontsize=18)
plt.ylabel("Formation Energy (eV/atom)", fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend()
plt.tight_layout()
plt.savefig("hull.png", dpi=300)
plt.show()
