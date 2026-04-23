import os
import argparse
import shutil
import numpy as np
import torch
import time
import subprocess
from modules.syntropizer import syntropizer_EAT
from modules.run_dft import run_dft
import modules.logging as logging
import modules.folder_ops as folder_ops
import pickle

# Note: the lattice constant for fixed ref in spline calc is crucial to know for setting the pure ref calc in EAT opt --- here it is 7.10 Bohr
with open("/data2/jt577/2025/November_2025/Phase_Diagram/EAT/spline_E.pkl", "rb") as f:
    spline_E = pickle.load(f)

with open("/data2/jt577/2025/November_2025/Phase_Diagram/EAT/spline_dE_Cu.pkl", "rb") as f:
    spline_dE_Cu = pickle.load(f)

with open("/data2/jt577/2025/November_2025/Phase_Diagram/EAT/spline_dE_Pd.pkl", "rb") as f:
    spline_dE_Pd = pickle.load(f)

def obj_fun(x_full, energies, energy_grads, dft_names):
    """
    Objective function to be maximized.
    x_full: (n_atoms, n_species) tensor (stoichiometry matrix)
    """
    # First species must be Cu and second must be Pd
    h2eV = 27.211386245988
    n_atoms = x_full.shape[0]

    x_np = x_full.detach().numpy()
    n_species = x_np.shape[1]

    energy_alloy = energies[dft_names[0]]
    energy_pures = sum([energies[dft_name] for dft_name in dft_names[1:]]) / 4 # div by 4 since 4 atoms per pure cell

    energy_grad_alloy = energy_grads[dft_names[0]].reshape(n_atoms, -1)
    energy_grad_pures = [energy_grads[dft_name] for dft_name in dft_names[1:]]
    energy_grad_pures = np.array(energy_grad_pures)

    # Evaluate energy difference with spline, treating 'x' in spline as the concentration of first species (Cu). These are already per atom.
    x_Cu = x_np[:,0]
    dE_Cu = spline_dE_Cu(x_Cu)
    dE_Pd = spline_dE_Pd(x_Cu)
    grad_delta_E = np.array([dE_Cu, dE_Pd]).T  # Shape (n_atoms, n_species)
    delta_E = spline_E(x_Cu)

    # Get true energy and gradient using the spline estimates and the energy_pures
    energy_pures_real = energy_pures + np.sum(delta_E)
    energy_grad_pures_real = grad_delta_E + energy_grad_pures

    energy = (energy_alloy - energy_pures_real)
    energy_grad = (energy_grad_alloy - energy_grad_pures_real).flatten()

    return -h2eV * (energy), -h2eV * (energy_grad)


# Obtain paths of input file and output directory from command line arguments
parser = argparse.ArgumentParser(
    description="Run effective atom theory (EAT)"
)
parser.add_argument(
    "infile",
    help="path to the input file"
)
parser.add_argument(
    "out_dir",
    help="directory where results will be written"
)
args = parser.parse_args()
infile  = args.infile
out_dir  = args.out_dir
if os.path.exists(out_dir):
    raise FileExistsError(f"The output directory already exists. Please choose a different name.")

# read input_file
with open(infile, 'r') as f:
    lines = f.readlines()
n_alpha = None
x0 = None
# Parse input file
for i, line in enumerate(lines):
    stripped_line = line.split('#', 1)[0].rstrip()
    if stripped_line:
        split_line = stripped_line.split()
        if split_line[0] == 'dft_files':
            dft_infiles = []
            dft_names = []
            dft_infiles.append(lines[i+1].strip())  # First DFT infile is alloy
            dft_names.append('dft_alloy')
            dft_infiles.append(lines[i+2].strip())  # Second DFT infile is pure
            dft_names.append('dft_pure')
        elif split_line[0] == 'dft_parallel':
            dft_parallel = int(split_line[1])
        elif split_line[0] == 'n_iter':
            n_iter = int(split_line[1])
        elif split_line[0] == 'n_atoms':
            n_atoms = int(split_line[1])
        elif split_line[0] == 'n_species':
            n_species = int(split_line[1])
        elif split_line[0] == 'n_alpha':
            n_alpha = [int(val) for val in split_line[1:]]
            if len(n_alpha) != n_species:
                raise ValueError("Length of n_alpha must be equal to n_species.")
            if sum(n_alpha) != n_atoms:
                raise ValueError("Sum of n_alpha must be equal to n_atoms.")
        elif split_line[0] == 'q':
            q = float(split_line[1])
        elif split_line[0] == 'eta_init':
            eta_init = float(split_line[1])
        elif split_line[0] == 'n_anneals':
            n_anneals = int(split_line[1])
        elif split_line[0] == 'anneal_mult':
            anneal_mult = float(split_line[1])
        elif split_line[0] == 'BFGS_hist':
            BFGS_hist = int(split_line[1])
        elif split_line[0] == 'tol':
            tol = float(split_line[1])
        elif split_line[0] == 'max_line_search':
            max_line_search = int(split_line[1])
        elif split_line[0] == 'jitter':
            jitter = float(split_line[1])
        elif split_line[0] == 'x0_path':
            x0_path = split_line[1]
            x0 = torch.load(x0_path)

# Create output directory if it doesn't exist
plot_dir = folder_ops.create_subfolder(out_dir, 'plots')  # Create new output directory
outfile = os.path.join(out_dir, 'output.txt')
d = n_atoms * n_species

# Log input file into output file
with open(infile, 'r') as f:
    infile_content = f.read()
logging.log("Input file\n", log_file=outfile)
logging.log(infile_content, log_file=outfile)
logging.log("--------------------------------------------------------------------------------------------------------------\n", log_file=outfile)


# Set up logging
logging.log(f"Effective Atom Theory\n", log_file=outfile)

# Syntropize
syntropizer_instance = syntropizer_EAT(func=obj_fun, N=n_atoms, S=n_species, n_alpha=n_alpha, log_file=outfile, out_dir=out_dir,
                                        dft_names=dft_names, dft_infiles=dft_infiles, dft_parallel=dft_parallel, 
                                        q=q, eta_init=eta_init, n_anneals=n_anneals, anneal_mult=anneal_mult, 
                                        BFGS_hist=BFGS_hist, n_iter=n_iter, tol=tol, max_line_search=max_line_search,
                                        jitter=jitter, x0=x0)
f_val, x = syntropizer_instance.objfun_min()