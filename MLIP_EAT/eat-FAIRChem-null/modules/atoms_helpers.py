# modules/atoms_helpers.py
from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers


def build_base_atoms(
    frac_coords_init,
    lattice_vectors,
    supercell_size,
    *,
    mix_placeholder: str = "H",
    null_default_for_mix: bool = False,
) -> Atoms:
    """
    frac_coords_init can be either:
      - old format: [[x, y, z], ...]
      - labeled format: [[label, x, y, z], ...]

    Supported labels:
      - "mix": standard EAT row with per-site weights summing to 1
      - "null": null row with per-site weights summing to <= 1
      - element symbols like "O", "Cu", ...
    """
    cell = np.array(lattice_vectors, dtype=float)

    symbols = []
    frac = []
    site_is_mix = []
    site_is_null = []
    site_fixed_Z = []

    for row in frac_coords_init:
        if len(row) == 3:
            label = "mix"
            x, y, z = row
        elif len(row) == 4:
            label, x, y, z = row
        else:
            raise ValueError(f"Bad frac_coords_init row: {row}")

        frac.append([float(x), float(y), float(z)])
        label_lower = str(label).lower()

        if label_lower in {"mix", "null"}:
            is_null = (label_lower == "null") or bool(null_default_for_mix)
            symbols.append(mix_placeholder)
            site_is_mix.append(1)
            site_is_null.append(1 if is_null else 0)
            site_fixed_Z.append(0)
        else:
            el = str(label)
            symbols.append(el)
            site_is_mix.append(0)
            site_is_null.append(0)
            site_fixed_Z.append(int(atomic_numbers[el]))

    atoms = Atoms(symbols=symbols, cell=cell, pbc=True)
    atoms.set_scaled_positions(np.array(frac, dtype=float))

    atoms.new_array("site_is_mix", np.array(site_is_mix, dtype=np.int32))
    atoms.new_array("site_is_null", np.array(site_is_null, dtype=np.int32))
    atoms.new_array("site_fixed_Z", np.array(site_fixed_Z, dtype=np.int32))

    sx, sy, sz = map(int, supercell_size)
    atoms = atoms.repeat((sx, sy, sz))
    return atoms
