from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class NullRepairResult:
    frac_all: np.ndarray
    moved_null_local_indices: list[int]
    deactivated_null_local_indices: list[int]
    remaining_close_pairs: int


def periodic_pairwise_dist2_frac(frac_a: np.ndarray, frac_b: np.ndarray, cell_np: np.ndarray) -> np.ndarray:
    frac_a = np.asarray(frac_a, dtype=float)
    frac_b = np.asarray(frac_b, dtype=float)
    if frac_a.size == 0 or frac_b.size == 0:
        return np.empty((len(frac_a), len(frac_b)), dtype=float)
    ds = frac_a[:, None, :] - frac_b[None, :, :]
    ds = ds - np.round(ds)
    dr = ds @ np.asarray(cell_np, dtype=float)
    return np.sum(dr**2, axis=-1)


def perturb_fractional_positions(
    frac_np: np.ndarray,
    cell_np: np.ndarray,
    sigma: float,
    movable_mask: np.ndarray | None = None,
) -> np.ndarray:
    frac_np = np.asarray(frac_np, dtype=float)
    if sigma <= 0.0 or frac_np.size == 0:
        return frac_np.copy()

    xyz = frac_np @ np.asarray(cell_np, dtype=float)
    if movable_mask is None:
        xyz = xyz + float(sigma) * np.random.randn(*xyz.shape)
    else:
        movable_mask = np.asarray(movable_mask, dtype=bool).reshape(-1)
        xyz = xyz.copy()
        xyz[movable_mask] = xyz[movable_mask] + float(sigma) * np.random.randn(int(np.count_nonzero(movable_mask)), 3)

    frac_out = np.linalg.solve(np.asarray(cell_np, dtype=float).T, xyz.T).T
    return np.mod(frac_out, 1.0)


def _count_close_pairs(frac_all: np.ndarray, cell_np: np.ndarray, cutoff2: float) -> int:
    if len(frac_all) < 2:
        return 0
    dist2 = periodic_pairwise_dist2_frac(frac_all, frac_all, cell_np)
    iu = np.triu_indices(len(frac_all), k=1)
    return int(np.count_nonzero(dist2[iu] < cutoff2))


def repair_null_overlaps(
    frac_all: np.ndarray,
    cell_np: np.ndarray,
    null_indices: np.ndarray,
    candidate_frac: np.ndarray,
    min_distance: float,
    *,
    priority_scores: np.ndarray | None = None,
    max_passes: int = 4,
) -> NullRepairResult:
    frac_all = np.asarray(frac_all, dtype=float).copy()
    null_indices = np.asarray(null_indices, dtype=int).reshape(-1)
    candidate_frac = np.asarray(candidate_frac, dtype=float)
    if min_distance <= 0.0 or null_indices.size == 0 or len(frac_all) < 2:
        return NullRepairResult(
            frac_all=frac_all,
            moved_null_local_indices=[],
            deactivated_null_local_indices=[],
            remaining_close_pairs=0,
        )

    if priority_scores is None:
        priority_scores = np.zeros(len(null_indices), dtype=float)
    else:
        priority_scores = np.asarray(priority_scores, dtype=float).reshape(-1)
        if priority_scores.size != null_indices.size:
            raise ValueError("priority_scores must match the number of null indices")

    null_local_by_atom = {int(atom_idx): local_idx for local_idx, atom_idx in enumerate(null_indices.tolist())}
    cutoff2 = float(min_distance) ** 2
    moved: set[int] = set()
    deactivated: set[int] = set()

    for _ in range(max(1, int(max_passes))):
        dist2 = periodic_pairwise_dist2_frac(frac_all, frac_all, cell_np)
        iu = np.triu_indices(len(frac_all), k=1)
        close_pairs = np.column_stack(iu)[dist2[iu] < cutoff2]
        if close_pairs.size == 0:
            break

        atoms_to_move: set[int] = set()
        for pair in close_pairs:
            i = int(pair[0])
            j = int(pair[1])
            i_local = null_local_by_atom.get(i)
            j_local = null_local_by_atom.get(j)

            if i_local is None and j_local is None:
                continue
            if i_local is not None and j_local is None:
                atoms_to_move.add(i)
                continue
            if i_local is None and j_local is not None:
                atoms_to_move.add(j)
                continue

            pri_i = float(priority_scores[i_local])
            pri_j = float(priority_scores[j_local])
            if pri_i < pri_j:
                atoms_to_move.add(i)
            elif pri_j < pri_i:
                atoms_to_move.add(j)
            else:
                atoms_to_move.add(max(i, j))

        if not atoms_to_move:
            break

        ordered_to_move = sorted(
            atoms_to_move,
            key=lambda atom_idx: (float(priority_scores[null_local_by_atom[int(atom_idx)]]), int(atom_idx)),
        )
        keep_mask = np.ones(len(frac_all), dtype=bool)
        keep_mask[np.asarray(ordered_to_move, dtype=int)] = False
        placed_frac = []

        for atom_idx in ordered_to_move:
            blocked_frac = frac_all[keep_mask]
            if placed_frac:
                blocked_frac = np.vstack([blocked_frac, np.asarray(placed_frac, dtype=float)])

            min_dist2 = np.full(len(candidate_frac), np.inf, dtype=float)
            if blocked_frac.size > 0:
                min_dist2 = periodic_pairwise_dist2_frac(candidate_frac, blocked_frac, cell_np).min(axis=1)

            feasible = np.where(min_dist2 >= cutoff2 - 1e-12)[0]
            if feasible.size > 0:
                candidate_idx = int(feasible[int(np.argmax(min_dist2[feasible]))])
            else:
                candidate_idx = int(np.argmax(min_dist2))
                deactivated.add(int(null_local_by_atom[int(atom_idx)]))

            frac_all[int(atom_idx)] = candidate_frac[candidate_idx]
            placed_frac.append(frac_all[int(atom_idx)].copy())
            keep_mask[int(atom_idx)] = True
            moved.add(int(null_local_by_atom[int(atom_idx)]))

    return NullRepairResult(
        frac_all=frac_all,
        moved_null_local_indices=sorted(moved),
        deactivated_null_local_indices=sorted(deactivated),
        remaining_close_pairs=_count_close_pairs(frac_all, cell_np, cutoff2),
    )
