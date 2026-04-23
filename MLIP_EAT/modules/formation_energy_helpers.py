#!/usr/bin/env python
"""
Formation energy calculation helpers.

Usage:
    python optimize_tc_with_formation_energy.py
"""

import sys
from pathlib import Path
import json
import numpy as np
import torch
from ase import Atoms

# Add utils to path
sys.path.insert(0, str(Path(__file__).parent / "utils"))
# Prefer Justin's patched fairchem if available
patched_fairchem = Path("/Users/justint/Library/CloudStorage/OneDrive-Personal/Desktop/Academic Stuff/Arias Research/Tc_optimization/December_2025/fairchem")
if patched_fairchem.exists():
    sys.path.insert(0, str(patched_fairchem))


# Import Justin's modified fairchem
try:
    from fairchem.core import FAIRChemCalculator, pretrained_mlip
    FAIRCHEM_AVAILABLE = True
except ImportError:
    FAIRCHEM_AVAILABLE = False
    FAIRChemCalculator = None  # For type hints
    pretrained_mlip = None  # For type hints

# Expose bulk for convenience
__all__ = ['bulk', 'FAIRCHEM_AVAILABLE', 'FAIRChemCalculator', 'pretrained_mlip']

def set_enable_eat_grad(atoms: Atoms, enabled: bool) -> Atoms:
    """
    Control whether EAT gradient computation is enabled for an Atoms object.
    
    When enabled, gradients with respect to EAT weights can be computed,
    which is useful for composition optimization.
    
    Args:
        atoms: ASE Atoms object to modify
        enabled: Whether to enable EAT gradient computation
    
    Returns:
        The modified Atoms object (same as input, for chaining)
    """
    atoms.info["enable_eat_grad"] = bool(enabled)
    return atoms

def create_eat_atoms_for_fairchem_multi(
    atoms: Atoms,
    elements: list[str],
    fractions: list[float],
    enable_grad: bool = True,
    dummy_element: str = "H",
) -> Atoms:
    """
    Create an ASE Atoms object configured for FAIRChem EAT (multi-element system).
    
    Args:
        atoms: Source atoms (for positions, cell, pbc)
        elements: List of element symbols
        fractions: List of fractions for each element (must sum to 1.0)
        enable_grad: Whether to enable EAT gradient computation
        dummy_element: Dummy element symbol to use (default: "H")
    
    Returns:
        New Atoms object configured for FAIRChem EAT
    """
    from pymatgen.core import Element
    
    if len(elements) != len(fractions):
        raise ValueError(f"elements and fractions must have same length")
    if abs(sum(fractions) - 1.0) > 1e-6:
        raise ValueError(f"fractions must sum to 1.0, got {sum(fractions)}")
    
    # Create atoms with dummy element
    atoms_eat = Atoms(
        symbols=[dummy_element] * len(atoms),
        positions=atoms.get_positions(),
        cell=atoms.get_cell(),
        pbc=atoms.get_pbc(),
    )
    
    # Get atomic numbers and create EAT weights
    MAX_NUM_ELEMENTS = 100
    n_atoms = len(atoms_eat)
    ew = np.zeros((n_atoms, MAX_NUM_ELEMENTS), dtype=np.float32)
    
    for elem, frac in zip(elements, fractions):
        z = Element(elem).Z
        ew[:, z] = frac
    
    atoms_eat.new_array("eat_weights", ew)
    
    # Set gradient flag
    set_enable_eat_grad(atoms_eat, enable_grad)
    
    return atoms_eat




class FairchemEnergyForceFunction(torch.autograd.Function):
    """
    Wrap a FAIRChem energy/force call as a torch autograd Function.
    Backward provides dE/dR = -forces and dE/d(eat_weights) from eat_grad.
    Supports arbitrary active species mapped by element_types.
    """
    @staticmethod
    def forward(ctx, positions, eat_weights, cell_tensor, atoms_template, active_elements, active_indices, element_types, fairchem_calc):
        # positions: (n_atoms, 3) tensor
        # eat_weights: (n_atoms, n_elem_model) tensor
        # cell_tensor: (3, 3) tensor
        from pymatgen.core import Element
        from ase.stress import voigt_6_to_full_3x3_stress
        
        pos_np = positions.detach().cpu().numpy()
        ew_np = eat_weights.detach().cpu().numpy()
        cell_np = cell_tensor.detach().cpu().numpy()
        n_atoms = ew_np.shape[0]

        # Build FAIRChem atoms with dummy symbols and full eat_weights indexed by Z
        MAX_NUM_ELEMENTS = 100
        eat_full = np.zeros((n_atoms, MAX_NUM_ELEMENTS), dtype=np.float32)
        for idx, elem in enumerate(element_types):
            if idx >= ew_np.shape[1]:
                continue
            z = Element(elem).Z
            eat_full[:, z] = ew_np[:, idx]

        from ase import Atoms
        atoms = Atoms(
            symbols=["H"] * n_atoms,
            positions=pos_np,
            cell=cell_np,
            pbc=atoms_template.get_pbc(),
        )
        atoms.new_array("eat_weights", eat_full)
        set_enable_eat_grad(atoms, True)
        atoms.calc = fairchem_calc

        # Forward FAIRChem call
        total_energy = atoms.get_potential_energy()
        forces = atoms.get_forces()  # numpy, shape (n_atoms, 3)
        eat_grad = atoms.calc.results.get("eat_grad", None)  # (n_atoms, MAX_NUM_ELEMENTS)
        stress_voigt = atoms.calc.results.get("stress", None)  # Voigt (6,)
        stress_full = None
        if stress_voigt is not None:
            stress_full = voigt_6_to_full_3x3_stress(stress_voigt)  # (3,3)

        # Map eat_grad (by Z) back to the model's eat_weight indices
        grad_eat_weights = np.zeros_like(ew_np, dtype=np.float32)
        if eat_grad is not None:
            for idx, elem in enumerate(element_types):
                if idx >= ew_np.shape[1]:
                    continue
                z = Element(elem).Z
                grad_eat_weights[:, idx] = eat_grad[:, z]

        forces_t = torch.tensor(forces, device=positions.device, dtype=positions.dtype)
        grad_eat_t = torch.tensor(grad_eat_weights, device=eat_weights.device, dtype=eat_weights.dtype)
        stress_t = torch.zeros((3, 3), device=positions.device, dtype=positions.dtype)
        if stress_full is not None:
            stress_t = torch.tensor(stress_full, device=positions.device, dtype=positions.dtype)
        volume_t = torch.tensor(atoms.get_volume(), device=positions.device, dtype=positions.dtype)

        ctx.save_for_backward(forces_t, grad_eat_t, stress_t, cell_tensor, volume_t)
        ctx.stress_available = stress_full is not None
        return torch.tensor(total_energy, device=positions.device, dtype=positions.dtype), forces_t, grad_eat_t, stress_t

    @staticmethod
    def backward(ctx, grad_energy, grad_forces, grad_eatgrad, grad_stress):
        forces_t, grad_eat_t, stress_t, cell_t, volume_t = ctx.saved_tensors

        # dE/dR = -forces; chain with upstream grad
        grad_positions = grad_energy * (-forces_t)
        grad_eat_weights = grad_energy * grad_eat_t
        grad_cell = None
        if cell_t.requires_grad and ctx.stress_available:
            # Gradient of energy w.r.t. cell matrix H using stress:
            # dE/dH = -V * sigma * H^{-T}
            cell_detached = cell_t.detach()
            cell_inv_t = torch.inverse(cell_detached).transpose(0, 1)
            grad_cell = grad_energy * (-volume_t) * torch.matmul(stress_t, cell_inv_t)

        # None for non-tensor inputs
        return grad_positions, grad_eat_weights, grad_cell, None, None, None, None, None


class FormationEnergyCalculator:
    """Calculate formation energy using reference energies."""
    
    def __init__(
        self, 
        reference_energies_path: str | Path | None = None, 
        reference_energies: dict | None = None,
    ):
        """
        Initialize with reference energies.
        
        Args:
            reference_energies_path: Path to JSON file with reference energies (optional)
            reference_energies: Dictionary of reference energies (optional, takes precedence)
            fallback_to_file: If True and reference_energies is provided, fall back to file
                             for missing elements. If False, only use provided reference_energies.
        """
        if reference_energies is not None:
            self.ref_energies = reference_energies.copy()
            # If fallback is enabled and a file path is provided, load file and merge
            if reference_energies_path is not None:
                try:
                    with open(reference_energies_path, 'r') as f:
                        file_ref_energies = json.load(f)
                    # Add any missing elements from file
                    for elem, energy in file_ref_energies.items():
                        if elem not in self.ref_energies:
                            self.ref_energies[elem] = energy
                except Exception as e:
                    print(f"Warning: Could not load fallback reference energies from {reference_energies_path}: {e}")
        elif reference_energies_path is not None:
            with open(reference_energies_path, 'r') as f:
                self.ref_energies = json.load(f)
        else:
            self.ref_energies = {}
    
    def compute_formation_energy(
        self,
        total_energy: float | torch.Tensor,
        composition: dict[str, float],
        num_atoms: int,
    ) -> float | torch.Tensor:
        """
        Compute formation energy per atom.
        
        Formation energy = (E_total - sum(n_i * E_ref_i)) / n_atoms
        
        For EAT (fractional compositions):
        Formation energy = (E_total - sum(f_i * E_ref_i * n_atoms)) / n_atoms
                         = E_total/n_atoms - sum(f_i * E_ref_i)
        
        Args:
            total_energy: Total energy from calculator (eV)
            composition: Dictionary mapping element symbols to fractions (must sum to 1.0)
            num_atoms: Number of atoms in the structure
        
        Returns:
            Formation energy per atom (eV/atom)
        """
        # Compute reference energy contribution
        ref_energy_contribution = 0.0
        is_torch = torch.is_tensor(total_energy)
        
        for element, fraction in composition.items():
            if element not in self.ref_energies:
                available = sorted(self.ref_energies.keys())
                raise ValueError(
                    f"Reference energy not found for element '{element}'. "
                    f"Available elements: {available[:20]}{'...' if len(available) > 20 else ''}"
                )
            
            ref_energy = self.ref_energies[element]
            if is_torch:
                ref_energy = torch.tensor(ref_energy, dtype=total_energy.dtype, device=total_energy.device)
            
            ref_energy_contribution += fraction * ref_energy
        
        # Formation energy per atom
        energy_per_atom = total_energy / num_atoms
        formation_energy = energy_per_atom - ref_energy_contribution
        
        return formation_energy
    
    def compute_formation_energy_from_eat_weights(
        self,
        total_energy: float | torch.Tensor,
        eat_weights: torch.Tensor,
        element_types: tuple[str, ...],
        num_atoms: int,
    ) -> float | torch.Tensor:
        """
        Compute formation energy from EAT weights tensor.
        
        Args:
            total_energy: Total energy from calculator (eV)
            eat_weights: Tensor of shape (num_atoms, len(element_types)) with EAT weights
            element_types: Tuple of element symbols matching the model's element_types
            num_atoms: Number of atoms in the structure
        
        Returns:
            Formation energy per atom (eV/atom)
        """
        # Average EAT weights across atoms to get composition fractions
        # eat_weights shape: (num_atoms, len(element_types))
        composition_fractions = eat_weights.mean(dim=0)  # (len(element_types),)
        
        # Convert to dict
        composition = {
            element_types[i]: composition_fractions[i].item() if torch.is_tensor(composition_fractions[i]) else composition_fractions[i]
            for i in range(len(element_types))
            if composition_fractions[i] > 1e-6  # Only include significant fractions
        }
        
        return self.compute_formation_energy(total_energy, composition, num_atoms)