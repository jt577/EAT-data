##################################################################################
# Main file for EAT O, OH adsorption energy minimization to target
##################################################################################

import os
import numpy as np
import modules.bsruncalc as bs
import modules.egread as eg
import modules.minimize as mn
import modules.tqentropy as tq
import modules.cdmfolders as cdm

# User defined parameters
elements = ['Ti', 'Cu', 'Ni', 'Zn']
target = -0.014699721440278707 # target binding energy (Hartree) (1.6 eV)
q = 2 # Tsallis q param
num_epochs = 1 # number of epochs
maxit = 1 # max num iterations per epoch

# Input files for 
def generic_inputs_init(adsorbate):
    res = (
    f"initial-state {adsorbate}.$VAR\n"
    f"ion-species SG15/$ID_ONCV_PBE.upf\n"
    f"elec-smearing Fermi 0.001\n"
    f"kpoint 0.5 0.5 0.5  1\n"
    f"kpoint-folding 3 3 1\n"
    f"elec-cutoff 30\n"
    f"elec-ex-corr gga-x-rpbe gga-c-pbe\n"
    f"symmetries none\n"
    f"# Relaxation commands\n"
    f"electronic-minimize energyDiffThreshold 1e-6 nIterations 1000 nEnergyDiff 5\n"
    f"ionic-minimize energyDiffThreshold 1e-4 nIterations 1000 nEnergyDiff 3\n"
    f"\n"
    f"dump-name {adsorbate}.$VAR\n"
    f"dump End Lattice IonicPositions Ecomponents\n"
    f"\n")

    return res

# Input files
def generic_inputs(adsorbate):
    res = (
    f"initial-state {adsorbate}.$VAR\n"
    f"ion-species SG15/$ID_ONCV_PBE.upf\n"
    f"elec-smearing Fermi 0.001\n"
    f"kpoint 0.5 0.5 0.5  1\n"
    f"kpoint-folding 3 3 1\n"
    f"elec-cutoff 30\n"
    f"elec-ex-corr gga-x-rpbe gga-c-pbe\n"
    f"symmetries none\n"
    f"# Relaxation commands\n"
    f"electronic-minimize energyDiffThreshold 1e-6 nIterations 1000 nEnergyDiff 5\n"
    f"\n"
    f"dump-name {adsorbate}.$VAR\n"
    f"dump End Lattice IonicPositions State MixGrad Ecomponents\n"
    f"\n")

    return res

# specify spatial structure (here we want rock salt)
lattice_constant = 7.88 
super_x = lattice_constant * np.array([1, 1, 0])
super_y = lattice_constant * np.array([1, -1, 0])
super_z = np.array([0, 0, 3/2 * 1/0.3 * lattice_constant])
# bounds of cell - columns are superlattice vectors
xbox = f'{super_x[0]} {super_y[0]} {super_z[0]}'
ybox = f'{super_x[1]} {super_y[1]} {super_z[1]}'
zbox = f'{super_x[2]} {super_y[2]} {super_z[2]}'
# define lattice explicitly
lattice = (
f"lattice \\ \n"
f"{xbox}\\ \n"
f"{ybox}\\ \n"
f"{zbox} \n")

# positions of atoms
positions_file = 'positions.txt'
with open(positions_file, 'r') as file:
    position_lines = file.readlines()

# Initialize a counter for lines containing 'mix'
mix_count = 0
# Iterate through each line and check for 'mix'
mixed_lines = []
nonmixed_lines = []
for line in position_lines:
    if 'mix' in line:
        mixed_lines.append(line)
        mix_count += 1
    else:
        if 'ion' in line:
            nonmixed_lines.append(line)
# order positions so that all mixed are first, then nonmixed
positions_ordered = mixed_lines + nonmixed_lines 
# positions of adsorbates
ads_O_pos = [f'ion O   0.0   0.0   0.23 1']
ads_OH_pos = [f'ion O   0.0   0.0   0.23 1']
ads_OH_pos += [f'ion H   0.1   0.1   0.25 1']


#################################################################################
# Functions relevant to cost minimization
#################################################################################

# cost function is the square of difference between current binding energy and platinum surface binding energy
def binding(w, folder, elements, target, eta, q):
    S = len(elements)
    # get binding energy of HEA
    MixBinding = eg.read_energy(adsorbate='O', folder=folder) - eg.read_energy(adsorbate='OH', folder=folder) + 1/2 *  -1.1781008671071755 - eV2h(0.36)
    # return cost fn and binding energy
    return MixBinding

# cost function is the square of difference between current binding energy and platinum surface binding energy
def cost(w, folder, elements, target, eta, q, first_iteration, generic_inputs, lattice, positions_ordered, ads_O_pos, ads_OH_pos):
    S = len(elements)
    # get binding energy of HEA
    MixBinding = eg.read_energy(adsorbate='O', folder=folder) - eg.read_energy(adsorbate='OH', folder=folder) + 1/2 *  -1.1781008671071755 - eV2h(0.36)
    # return cost fn and binding energy
    return 1/2 * (h2eV(MixBinding - target))**2 + eta * tq.Tsallis(w, S, q)

# gradient of cost function
def grad_cost(w, folder, elements, target, eta, q, first_iteration, generic_inputs, lattice, positions_ordered, ads_O_pos, ads_OH_pos):
    S = len(elements)
    # get binding energy of HEA
    MixBinding = eg.read_energy(adsorbate='O', folder=folder) - eg.read_energy(adsorbate='OH', folder=folder) + 1/2 * -1.1781008671071755 - eV2h(0.36)
    # get gradient of binding energy
    grad_binding = eg.read_gradient(w, adsorbate='O', elements=elements, folder=folder) - eg.read_gradient(w, adsorbate='OH', elements=elements, folder=folder)

    res = h2eV(MixBinding - target) * h2eV(grad_binding) + eta * tq.grad_Tsallis(w, S, q)
    # return grad of cost fn
    return res

# convert units from Hartree to eV
def h2eV(energy):
    return energy*27.2114

# convert units from eV to Hartree
def eV2h(energy):
    return energy/27.2114

# Custom callback function
def callback_func(x, folder):
    global target
    global iteration_counter
    global elements
    global eta
    global lattice_constant

    MixBinding = binding(x, folder, elements, target, eta=0, q=1)
    S = len(elements)
    weights = x.reshape(-1, S)
    entropy = tq.Tsallis(x, S, q=1)

    # Write the current solution and its objective value into a file
    with open('min_progress.txt', 'a') as file:
        file.write(f'\nIteration {iteration_counter}, eta = {eta}\n')
        for i in range(weights.shape[0]):
            weight_string = ''
            for j in range(S):
                weight_string += f'{weights[i,j]} {elements[j]} '
            file.write(f'Atom {i+1} weights: {weight_string}\n')

        file.write(f'HEA binding free energy: {h2eV(MixBinding)} eV, Target: {h2eV(target)} eV\n')
        file.write(f'Entropy: {entropy}\n')
    # update lattice and ion position in position files
    cdm.update_pos(iteration_counter)
    # incremement the iteration
    iteration_counter += 1


###################################################################
# Main minimization loop
###################################################################

# delete progress files
progress_file_name = 'min_progress.txt'
details_file_name = 'min_details.txt'
cdm.delete_progress(progress_file_name, details_file_name)

# Initialize variables
num_species = len(elements)
eta = 1e-5 # Set entropy constant eta to small value

# w_full = np.zeros((mix_count, num_species))
# for i in range(mix_count):
#     w_full[i,:] = np.random.uniform(0,1,num_species)
#     w_full[i,:] = w_full[i,:]/np.sum(w_full[i,:])
# w = w_full.flatten()
# x = w

# # Read in initial weights
# with open('initial_weights.txt') as weights_file:
#     weights_lines = weights_file.readlines()
# w_full = np.zeros((mix_count, num_species))
# for i, weights_line in enumerate(weights_lines):
#     split_weights_line = weights_line.strip().split()
#     w_full[i, :] = np.array([float(split_weights_line[2]), float(split_weights_line[3]), float(split_weights_line[4]), float(split_weights_line[5])])
# w = w_full.flatten()
# x = w.copy()

w = np.zeros((1, 4))
x = w.copy()

quasi_Hessian = None
iteration_counter = 1
unique_folder = None
progress_file_path = os.path.join(os.getcwd(), progress_file_name)
details_file_path = os.path.join(os.getcwd(), details_file_name)

# initialize minimization progress file
element_string = ''
for i in range(len(elements)):
    element_string += f'{elements[i]} '
with open(progress_file_name, 'a') as file:
    file.write(f'Initilizing minimization for species {element_string}\n\n')
    file.write(f'Target: {h2eV(target)} eV\n')

with open(progress_file_name, 'a') as file:
    file.write(f'Performing first surface calculations...\n')

first_iteration = True
# Perform first surface calculation 
unique_folder = bs.perform_calc(w, elements, target, eta, q, first_iteration, generic_inputs_init, lattice, positions_ordered, ads_O_pos, ads_OH_pos)
first_iteration = False

with open(progress_file_name, 'a') as file:
    file.write(f'Surface calculations successfully finished.\n')

with open(progress_file_name, 'a') as file:
    file.write(f'Initial entropy parameter eta set to {eta}\n')
    file.write(f'_____________________________________________________________________\nBeginning minimzation\n\n')

# Initilialize min progress files
with open(details_file_name, 'a') as file:
    file.write(f"Beginning minimization")

# Main loop
for i in range(num_epochs):
    with open(details_file_name, 'a') as file:
        file.write("\n__________________________________________________________________________________\n")
        file.write(f"Epoch {i+1}: Loss = {cost(x, unique_folder, elements, target, 0, q, first_iteration, generic_inputs, lattice, positions_ordered, ads_O_pos, ads_OH_pos)}, Entropy = {tq.Tsallis(x, S=num_species, q=1)}, eta = {eta}\n")
    x, _, unique_folder = mn.PBFGS(fun=cost, x0=x, jac=grad_cost, args=(elements, target, eta, q, first_iteration, generic_inputs, lattice, positions_ordered, ads_O_pos, ads_OH_pos), S=num_species, B=quasi_Hessian, folder=unique_folder, progress_file_path=progress_file_path, details_file_path=details_file_path, maxit=maxit, tol=1e-2, callback=callback_func)
    eta *= 10
    if tq.Tsallis(x, S=num_species, q=1) < 1e-2:
        break
with open(details_file_name, 'a') as file:
    file.write("\n__________________________________________________________________________________\n")
    file.write(f"FINAL: Loss = {cost(x, unique_folder, elements, target, 0, q, first_iteration, generic_inputs, lattice, positions_ordered, ads_O_pos, ads_OH_pos)}, Entropy = {tq.Tsallis(x, S=num_species, q=1)}, eta = {eta}\n")
