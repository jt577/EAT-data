###################################################################
# Functions that read files from jdftx outputs and obtain energy and gradients wrt weights
###################################################################

# imports
import os
import numpy as np

# read energy file and extract free energy value
def read_energy(adsorbate, folder):
    prefix = f'{adsorbate}'
    subfolder = os.path.join(folder, prefix)
    energy_file = os.path.join(subfolder, f'{prefix}.Ecomponents')
    with open(energy_file, 'r') as file:
        for line_number, line in enumerate(file, 1):
            stripped_line = line.strip()
            if 'F' in stripped_line:
                split_line = stripped_line.split()
                energy = split_line[2]
                break

    return float(energy)

# gradient of energy wrt weights corr to atom #i
def gradient_per_atom(file_name, i, elements):
    num_species = len(elements)
    grad = np.zeros(num_species)
    with open(file_name, "r") as file:
        for line_number, line in enumerate(file, 1):
            stripped_line = line.strip()
            if f'mix{i}' in stripped_line:
                split_line = stripped_line.split()
                for j in range(num_species):
                    grad[j] = float(split_line[2+j])
                break
    return grad


# gradient of energy wrt weights
def read_gradient(w, adsorbate, elements, folder):
    prefix = f'{adsorbate}'
    num_species = len(elements)
    weights = w.reshape(-1, num_species)
    # extract gradient
    num_atoms = weights.shape[0]
    grad_w = np.zeros((num_atoms, num_species))
    grad_file = os.path.join(folder, prefix, f'{prefix}.mixgrad')
    for i in range(num_atoms):
        grad_w[i, :] = gradient_per_atom(grad_file, i+1, elements)
    # vectorize
    res = grad_w.flatten()

    return res
