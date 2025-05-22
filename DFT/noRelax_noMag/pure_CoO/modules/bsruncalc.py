###############################################################################
# Functions related to bulk and surface calcualtions and extracting energy and energy gradient
###############################################################################

# imports
import subprocess
import os
import numpy as np
import shutil
import time
import modules.cdmfolders as cdm

# writes the bulk input file for certain specified elements and weights
def write_input_surface(weights, adsorbate, elements, folder, first_iteration, generic_inputs, lattice, positions_ordered, ads_O_pos, ads_OH_pos):
    mix_count = 0
    num_species = 4
    prefix = f'{adsorbate}'

    # generic input params
    input_file_0 = generic_inputs(adsorbate)

    # initialize
    input_file_2 = []

    if first_iteration == True:

        # positions of mixed atoms
        mixed_pos = positions_ordered[:mix_count]
        nonmixed_pos = positions_ordered[mix_count:]
    
        for i in range(mix_count):
            weight_string = ''
            for j in range(num_species):
                weight_string += f'{weights[i,j]} '
            curr_line = mixed_pos[i]
            split_line = curr_line.split()
            curr_mixed_pos = split_line[2:]
            curr_mixed_pos_str = ''
            for entry in curr_mixed_pos:
                curr_mixed_pos_str += entry + '    '
            input_file_2 += [f'ion mix{i+1}   {weight_string}    {curr_mixed_pos_str}']
        for entry in nonmixed_pos:
            input_file_2 += [entry]

        # Oxygen adatom
        if adsorbate == 'O':
            input_file_2 += ads_O_pos 
        elif adsorbate == 'OH':
            input_file_2 += ads_OH_pos

        input_file_0 += lattice

    else:
        mix_count = weights.shape[0]
        num_species = weights.shape[1]

        ionpos_path = os.path.join('runs', 'positions', f'{prefix}.ionpos')
        # extract ion positions and lattice
        with open(ionpos_path) as file:
            # Iterate through each line and check for 'mix'
            mixed_lines = []
            nonmixed_lines = []
            for line in file:
                if 'mix' in line:
                    mixed_lines.append(line)
                else:
                    if 'ion' in line and '#' not in line:
                        nonmixed_lines.append(line)
        # order positions so that all mixed are first, then nonmixed
        positions_ordered = mixed_lines + nonmixed_lines 
        lattice_path = os.path.join('runs', 'positions', f'{prefix}.lattice')
        with open(lattice_path) as file:
            lattice = file.readlines()

        # add lattice
        input_file_0 += ''.join(lattice)
        
        for i in range(mix_count):  
            weight_string = ''
            for j in range(num_species):
                weight_string += f'{weights[i,j]} '
            curr_line = mixed_lines[i]
            split_line = curr_line.split()
            curr_mixed_pos = split_line[2:]
            curr_mixed_pos_str = ''
            for entry in curr_mixed_pos:
                curr_mixed_pos_str += entry + '    '
            input_file_2 += [f'ion mix{i+1}   {weight_string}    {curr_mixed_pos_str}']
        for entry in nonmixed_lines:
            input_file_2 += [entry]

    # add mixed ions
    input_file_1 = []
    element_string = ''
    for i in range(num_species):
        element_string += f'{elements[i]} '

    for i in range(mix_count):
        input_file_1 += [f'add-mix mix{i+1} {element_string}']

    # combine lines
    input_file = input_file_0 + '\n' + '\n'.join(input_file_1) + '\n\n' + '\n'.join(input_file_2)
    
    # write the file
    working_folder = os.path.join(folder, prefix)
    input_file_path = os.path.join(working_folder, f'{prefix}.in')
    try:
        with open(input_file_path, "w") as f:
            f.write(input_file)
    except Exception as e:
        print(f"An error occurred: {e}")

def run_jdftx(input_file, output_file):
    cwd = os.getcwd()
    # The --oversubscribe option is sometimes needed if Open MPI finds fewer slots than requested.
    command = (
        f"mpirun --oversubscribe --bind-to none -n 9 "
        f"/data2/jt577/jdftx_eat_withLibXC/build/jdftx -i {cwd}/{input_file} | tee {cwd}/{output_file}"
    )
    subprocess.run(command, shell=True, check=True)


def energy_surface(w, adsorbate, elements, folder, first_iteration, generic_inputs, lattice, positions_ordered, ads_O_pos, ads_OH_pos):
    prefix = f'{adsorbate}'
    num_species = len(elements)
    weights = w.reshape(-1, num_species)

    write_input_surface(weights, adsorbate, elements, folder, first_iteration, generic_inputs, lattice, positions_ordered, ads_O_pos, ads_OH_pos)
    original_cwd = os.getcwd()
    try:
        working_folder = os.path.join(folder, prefix)
        os.chdir(working_folder)
        # Run the DFT calc; this will block until complete.
        run_jdftx(f'{adsorbate}.in', f'{adsorbate}.out')
    finally:
        os.chdir(original_cwd)
    # No job_id to return since the call was synchronous.
    return None


def perform_calc(w, elements, target, eta, q, first_iteration, generic_inputs, lattice, positions_ordered, ads_O_pos, ads_OH_pos):
    base_path = os.getcwd()
    run_folder = cdm.create_subfolder(folder_path=base_path, subfolder_name='runs')
    unique_folder = cdm.create_unique_folder(base_path=run_folder, prefix='calc_')

    O_path = cdm.create_subfolder(folder_path=unique_folder, subfolder_name='O')
    OH_path = cdm.create_subfolder(folder_path=unique_folder, subfolder_name='OH')
    wfns_path = cdm.create_subfolder(folder_path=run_folder, subfolder_name='wavefunctions')
    cdm.mv_wfns_to_unique(source_folder=wfns_path, destination_folders=[O_path, OH_path])

    # Each energy calculation now blocks until finished.
    energy_surface(w, adsorbate='O', elements=elements, folder=unique_folder, first_iteration=first_iteration, generic_inputs=generic_inputs, lattice=lattice, positions_ordered=positions_ordered, ads_O_pos=ads_O_pos, ads_OH_pos=ads_OH_pos)
    energy_surface(w, adsorbate='OH', elements=elements, folder=unique_folder, first_iteration=first_iteration, generic_inputs=generic_inputs, lattice=lattice, positions_ordered=positions_ordered, ads_O_pos=ads_O_pos, ads_OH_pos=ads_OH_pos)

    wfns_path = cdm.create_subfolder(folder_path=run_folder, subfolder_name='wavefunctions')
    cdm.mv_wfns_from_unique(source_folder=unique_folder, destination_folder=wfns_path)

    pos_path = cdm.create_subfolder(folder_path=run_folder, subfolder_name='positions')
    cdm.mv_pos_surface(source_folder=unique_folder, destination_folder=pos_path)

    return unique_folder