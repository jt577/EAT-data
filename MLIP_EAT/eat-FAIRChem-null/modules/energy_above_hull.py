#!/usr/bin/env python3
"""
HullEnergyCalculator: energy on the convex hull + gradient wrt composition
for a fixed design space, using a Materials Project-derived dictionary.

Expected JSON format (mp_materials.json):

{
  "Fe1O1": {
    "elements": ["Fe", "O"],
    "fractions": [0.5, 0.5],
    "formation_energy_per_atom": -1.23
  },
  ...
}

Usage pattern:

    from hull_energy_calculator import HullEnergyCalculator

    design_space = ["H", "O", "C"]
    calc = HullEnergyCalculator("mp_materials.json", design_space)

    c_star = [0.2, 0.3, 0.5]
    E_hull, grad_c = calc.get_energy_and_gradient(c_star)

    # E_hull: scalar hull energy at c_star (eV/atom)
    # grad_c: np.array(m,), dE_hull/dc_i in order of design_space
"""
from __future__ import annotations
import json
from pathlib import Path
from typing import List, Tuple, Optional, Union, Any, Dict

import numpy as np
import torch
from scipy.optimize import linprog
from dataclasses import dataclass, asdict
import time
import uuid

ArrayLike = Union[List[float], np.ndarray]

@dataclass
class PhaseRecord:
    # required for hull
    elements: List[str]
    fractions: List[float]
    formation_energy_per_atom: float

    # reproducibility payload (recommended)
    total_energy_ev: Optional[float] = None
    n_atoms: Optional[int] = None
    symbols: Optional[List[str]] = None
    frac_coords: Optional[List[List[float]]] = None  # (n_atoms, 3)
    lattice_vectors: Optional[List[List[float]]] = None  # (3, 3)
    cell_pbc: Optional[List[bool]] = None
    supercell_size: Optional[List[int]] = None

    # optional provenance
    metadata: Optional[Dict[str, Any]] = None


class HullEnergyCalculator:
    """
    Precomputes the subset of Materials Project phases relevant to a given
    design space and provides hull energy + gradient wrt composition.

    Attributes
    ----------
    design_space : list of str
        Elements defining the compositional basis, e.g. ["H", "O", "C"].
    m : int
        Number of elements in the design space.
    C : np.ndarray, shape (N, m)
        Composition matrix of all valid phases in the design space.
    E : np.ndarray, shape (N,)
        Formation energies per atom for those phases.
    """

    def __init__(
        self,
        db_json_path: Union[str, Path],
        design_space: List[str],
    ):
        """
        Initialize the calculator for a specific design space.

        Parameters
        ----------
        db_json_path : str or Path
            Path to the full MP-derived JSON database (e.g. "mp_materials.json").
        design_space : list of str
            Elements defining the compositional basis, in the desired order.
        """
        self.db_json_path = Path(db_json_path)
        self.design_space = list(design_space)
        self.m = len(self.design_space)

        # Load full database
        db = self._load_phase_database(self.db_json_path)

        # Build sub-database for this design space
        sub_db = self._build_subdatabase_for_design_space(db, self.design_space)
        if len(sub_db) == 0:
            raise RuntimeError(
                f"No phases found in database for design space {self.design_space}."
            )

        # Convert sub-db to phase list, then build C and E
        phases = self._phases_from_subdatabase(sub_db)
        self.C, self.E = self._build_composition_matrix(phases, self.design_space)

    # ======================================================================
    # Public API
    # ======================================================================

    def get_energy_and_gradient(
        self,
        c_star: ArrayLike,
    ) -> Tuple[Optional[float], Optional[np.ndarray]]:
        """
        Compute hull energy and gradient wrt composition for a given composition.

        Parameters
        ----------
        c_star : array-like of shape (m,)
            Composition vector in the same order as self.design_space.
            Must sum to 1.

        Returns
        -------
        E_hull : float or None
            Hull energy (eV/atom) at c_star, or None on failure.
        grad_c : np.ndarray of shape (m,) or None
            Gradient dE_hull/dc_i in the same order as self.design_space,
            or None if gradient unavailable or solve failed.
        """
        C = self.C
        E = self.E
        c_star = np.asarray(c_star, dtype=float)

        if c_star.shape != (self.m,):
            raise ValueError(
                f"c_star must have shape ({self.m},), got {c_star.shape}"
            )

        if not np.isclose(c_star.sum(), 1.0, atol=1e-6):
            raise ValueError(
                f"c_star must sum to 1.0; got sum={c_star.sum()}."
            )

        E_hull, grad_c = self._hull_energy_and_gradient_from_matrices(C, E, c_star)
        return E_hull, grad_c

    # ======================================================================
    # Internal helpers: DB + phase handling
    # ======================================================================

    @staticmethod
    def _load_phase_database(path: Path):
        """Load the full materials dictionary from JSON."""
        if not path.exists():
            raise FileNotFoundError(f"Database JSON not found at: {path}")
        with path.open("r") as f:
            data = json.load(f)
        return data

    @staticmethod
    def _build_subdatabase_for_design_space(db_dict, design_space: List[str]):
        """
        Build a smaller dictionary that only contains entries whose elements
        are a subset of the given design_space.
        """
        design_set = set(design_space)
        sub_db = {}

        for label, entry in db_dict.items():
            elems = entry["elements"]
            if set(elems).issubset(design_set):
                sub_db[label] = entry

        return sub_db

    @staticmethod
    def _phases_from_subdatabase(sub_db):
        """
        Return one phase per label for hull-building.
        If label has multiple entries, use the MIN formation energy among them.
        """
        phases = []
        for label, bucket in sub_db.items():
            elems = bucket["elements"]
            fracs = bucket["fractions"]

            if "entries" in bucket:
                energies = [float(e["formation_energy_per_atom"]) for e in bucket["entries"]]
                E_min = float(np.min(energies)) if len(energies) > 0 else float("inf")
            else:
                # backward compatibility (old schema)
                E_min = float(bucket["formation_energy_per_atom"])

            phases.append(
                {
                    "label": label,
                    "elements": elems,
                    "fractions": fracs,
                    "E_form": E_min,
                }
            )
        return phases

    @staticmethod
    def _normalize_fractions(elements: List[str], fractions: List[float]) -> List[float]:
        f = np.asarray(fractions, dtype=float)
        if f.ndim != 1 or len(f) != len(elements):
            raise ValueError("fractions must be length == elements")
        s = float(f.sum())
        if s <= 0:
            raise ValueError("fractions sum must be > 0")
        f = f / s
        return f.tolist()


    def add_phase_to_db(
        self,
        label: str,
        record: PhaseRecord,
    ) -> None:
        """
        Append a phase record under `label` (store ALL entries).
        Hull matrices will still be built from the min-energy entry per label.
        """
        db = self._load_phase_database(self.db_json_path)

        # normalize fractions
        record.fractions = self._normalize_fractions(record.elements, record.fractions)

        entry = asdict(record)
        entry["_uuid"] = str(uuid.uuid4())
        entry["_timestamp"] = time.time()

        if label not in db:
            # initialize a bucket for this composition label
            db[label] = {
                "elements": list(record.elements),
                "fractions": list(record.fractions),
                "entries": [],
            }
        else:
            # if the DB had the *old* single-record format, upgrade it on the fly
            if "entries" not in db[label]:
                prev = dict(db[label])
                db[label] = {
                    "elements": prev.get("elements", list(record.elements)),
                    "fractions": prev.get("fractions", list(record.fractions)),
                    "entries": [],
                }
                # if previous looked like a record, preserve it as one entry
                if "formation_energy_per_atom" in prev:
                    prev_entry = dict(prev)
                    prev_entry["_uuid"] = str(uuid.uuid4())
                    prev_entry["_timestamp"] = time.time()
                    db[label]["entries"].append(prev_entry)

        db[label]["entries"].append(entry)

        # write
        self.db_json_path.parent.mkdir(parents=True, exist_ok=True)
        with self.db_json_path.open("w") as f:
            json.dump(db, f, indent=2)

        # refresh hull matrices (min energy per label!)
        sub_db = self._build_subdatabase_for_design_space(db, self.design_space)
        phases = self._phases_from_subdatabase(sub_db)  # updated below
        self.C, self.E = self._build_composition_matrix(phases, self.design_space)



    @staticmethod
    def _build_composition_matrix(phases, design_space: List[str]):
        """
        From a list of phases + design_space, build:

            C: (N, m) matrix of phase compositions in design-space basis
            E: (N,)   vector of formation energies

        For each phase, embed its composition into the design-space basis:
          - columns correspond to design_space elements
          - entries are atomic fractions (sum to 1 for each row).
        """
        m = len(design_space)
        elem_to_idx = {el: i for i, el in enumerate(design_space)}

        N = len(phases)
        C = np.zeros((N, m), dtype=float)
        E = np.zeros(N, dtype=float)

        for i, ph in enumerate(phases):
            E[i] = float(ph["E_form"])
            elems = ph["elements"]
            fracs = ph["fractions"]

            for el, f in zip(elems, fracs):
                C[i, elem_to_idx[el]] += float(f)

        return C, E

    # ======================================================================
    # Internal core: LP solve + gradient from C, E
    # ======================================================================

    @staticmethod
    def _hull_energy_and_gradient_from_matrices(
        C: np.ndarray,
        E: np.ndarray,
        c_star: np.ndarray,
    ) -> Tuple[Optional[float], Optional[np.ndarray]]:
        """
        Compute hull energy and gradient wrt composition from precomputed
        composition matrix C and energy vector E.

        Args
        ----
        C : (N, m) ndarray
            Composition matrix in design-space basis.
        E : (N,) ndarray
            Formation energies per atom.
        c_star : (m,) ndarray
            Target composition, must sum to 1.

        Returns
        -------
        E_hull : float or None
            Hull energy at c_star (eV/atom).
        grad_c : np.ndarray of shape (m,) or None
            Gradient dE_hull/dc_star, or None if unavailable.
        """
        C = np.asarray(C, float)
        E = np.asarray(E, float)
        c_star = np.asarray(c_star, float)
        c_star = c_star / c_star.sum()  # ensure normalization

        N, m = C.shape
        if c_star.shape != (m,):
            raise ValueError(f"c_star must have shape ({m},), got {c_star.shape}")

        if not np.isclose(c_star.sum(), 1.0, atol=1e-6):
            raise ValueError(f"c_star must sum to 1.0; got sum={c_star.sum()}")

        # Objective: minimize E^T λ
        c = E.copy()

        # Equality constraints:
        #   A_eq @ λ = b_eq
        #   where:
        #       first m rows encode C^T λ = c_star
        #       last row encodes sum(λ) = 1
        A_eq = np.vstack([C.T, np.ones((1, N), dtype=float)])
        b_eq = np.concatenate([c_star, [1.0]])

        # Bounds: λ_i >= 0
        bounds = [(0.0, None)] * N

        res = linprog(
            c=c,
            A_eq=A_eq,
            b_eq=b_eq,
            bounds=bounds,
            method="highs"
        )

        if not res.success:
            print(f"[WARNING] linprog failed: {res.message}")
            return None, None

        lam = res.x
        min_lam = lam.min()

        # Thresholds for "how negative is too negative"
        soft_tol = -1e-6   # small negatives ≈ numerical noise
        hard_tol = -1e-2   # more negative than this is probably a real problem

        if min_lam < hard_tol:
            # Truly bad solution – treat as invalid
            print(f"[WARNING] Strongly negative λ encountered (min={min_lam:.3e}); treating as invalid.")
            return None, None

        if min_lam < soft_tol:
            # Mildly negative – likely numerical noise; clip to zero and renormalize
            print(f"[INFO] Mildly negative λ encountered (min={min_lam:.3e}); clipping and renormalizing.")
            lam = np.clip(lam, 0.0, None)
            s = lam.sum()
            if s <= 0:
                print("[WARNING] λ sum non-positive after clipping; treating as invalid.")
                return None, None
            lam /= s

        # Recompute hull energy from (possibly adjusted) λ
        E_hull = float(lam @ E)

        # Composition reconstruction sanity check (after any clipping)
        c_recon = C.T @ lam
        if not np.allclose(c_recon, c_star, atol=1e-4):
            print("[WARNING] Composition mismatch after optimization.")
            return None, None


        # Duals / gradient extraction:
        # For method="highs", SciPy exposes:
        #   res.eqlin.marginals: dual variables for equality constraints.
        #
        # Our A_eq is stacked as:
        #   [C^T; ones_row] @ λ = [c_star; 1]
        #
        # So:
        #   duals_eq[:m] correspond to C^T λ = c_star,
        #   duals_eq[m]   corresponds to sum(λ) = 1.
        #
        # The marginals are the derivative of the optimal objective value with
        # respect to the RHS entries in b_eq, so duals_eq[:m] = dE_hull/dc_star.
        try:
            duals_eq = res.eqlin.marginals  # shape: (m + 1,)
        except AttributeError:
            print("[WARNING] SciPy version does not expose eqlin.marginals; gradient unavailable.")
            return E_hull, None

        if duals_eq is None or len(duals_eq) != m + 1:
            print("[WARNING] Unexpected duals shape; gradient unavailable.")
            return E_hull, None

        grad_c = np.asarray(duals_eq[:m], dtype=float)

        return E_hull, grad_c

def calc_hull(x, hullcalculator, null_mask=None, fixed_counts=None, min_total=1e-8):
    counts = torch.sum(x, dim=0)
    if fixed_counts is not None:
        counts = counts + fixed_counts.to(device=x.device, dtype=x.dtype)

    total = torch.sum(counts)
    total_safe = torch.clamp(total, min=float(min_total))
    composition = counts / total_safe

    composition_np = composition.detach().cpu().numpy()
    E_hull, grad_c = hullcalculator.get_energy_and_gradient(composition_np)
    if E_hull is None or grad_c is None:
        return E_hull, None

    grad_c_t = torch.tensor(grad_c, dtype=x.dtype, device=x.device)
    grad_rows = grad_c_t.unsqueeze(0).expand_as(x) / total_safe

    if null_mask is not None:
        null_mask = torch.as_tensor(null_mask, device=x.device, dtype=torch.bool).view(-1)
        if null_mask.numel() != x.shape[0]:
            raise ValueError("null_mask must have one entry per row in x")
        if torch.any(null_mask):
            offset = torch.dot(grad_c_t, composition) / total_safe
            grad_rows = grad_rows.clone()
            grad_rows[null_mask] -= offset

    return E_hull, grad_rows.flatten()


# =============================================================================
# Optional: simple self-test when run as a script
# =============================================================================

if __name__ == "__main__":
    # Example quick test
    DB_JSON_PATH = "alexandria/alexandria_materials.json"  # Adjust path as needed
    design_space = ["Mg", "B", "O"]
    calc = HullEnergyCalculator(DB_JSON_PATH, design_space)

    compositions_to_test = [
        [2/6, 1/6, 1/2],
    ]

    print(f"Design space: {design_space}")
    for c_star in compositions_to_test:
        E_hull, grad_c = calc.get_energy_and_gradient(c_star)
        print("\n---------------------------------")
        print(f"Composition:  {c_star}")
        if E_hull is None:
            print("Hull energy not available.")
            continue
        print(f"Hull energy:  {E_hull:.6f} eV/atom")
        if grad_c is None:
            print("Gradient not available.")
        else:
            print("Gradient wrt composition:")
            for el, g in zip(design_space, grad_c):
                print(f"  dE_hull/dc_{el:>2s} = {g: .6f} eV/atom")
