# modules/atoms_helpers.py
from __future__ import annotations
import numpy as np
import torch
from ase import Atoms
from .lattice_helpers import cell_to_ucparams
from ase.data import atomic_numbers

def build_base_atoms(frac_coords_init, lattice_vectors, supercell_size, *, mix_placeholder: str = "H") -> Atoms:
    """
    frac_coords_init can be either:
      - old format: [[x,y,z], ...]
      - new format: [[label,x,y,z], ...] where label is "mix" or an element like "O"
    """
    cell = np.array(lattice_vectors, dtype=float)

    symbols = []
    frac = []
    site_is_mix = []
    site_fixed_Z = []

    for row in frac_coords_init:
        if len(row) == 3:
            # old format: treat as mix
            label = "mix"
            x, y, z = row
        elif len(row) == 4:
            label, x, y, z = row
        else:
            raise ValueError(f"Bad frac_coords_init row: {row}")

        frac.append([float(x), float(y), float(z)])

        if str(label).lower() == "mix":
            symbols.append(mix_placeholder)      # placeholder; EAT weights will define chemistry
            site_is_mix.append(1)
            site_fixed_Z.append(0)
        else:
            el = str(label)
            symbols.append(el)
            site_is_mix.append(0)
            site_fixed_Z.append(int(atomic_numbers[el]))

    atoms = Atoms(symbols=symbols, cell=cell, pbc=True)
    atoms.set_scaled_positions(np.array(frac, dtype=float))

    atoms.new_array("site_is_mix", np.array(site_is_mix, dtype=np.int32))
    atoms.new_array("site_fixed_Z", np.array(site_fixed_Z, dtype=np.int32))

    # supercell replication (arrays replicate too)
    sx, sy, sz = map(int, supercell_size)
    atoms = atoms.repeat((sx, sy, sz))
    return atoms
