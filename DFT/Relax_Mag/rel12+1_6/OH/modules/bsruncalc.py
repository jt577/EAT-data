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
def write_input_surface(weights, elements, folder, generic_inputs, lattice, mixed_positions, nonmixed_input, magmom, common_inputs):
    mix_count = 0
    num_species = 0
    
    # Initialize
    generated_lines = []  # List to store generated lines
    input_file_2 = []  # Final list of all lines

    # Generic input params
    if common_inputs is not None:
        prefix = 'bandstruct'
        input_file_0 = common_inputs
    else:
        prefix = 'run'
        input_file_0 = generic_inputs(prefix, magmom).splitlines(True)

    # Generate mixed ion lines
    for i in range(mix_count):
        weight_string = ''
        one_index = np.argmax(weights[i, :])
        curr_element = elements[one_index]
        line = (
            f'ion {curr_element} '
            f'{mixed_positions[i, 0]} {mixed_positions[i, 1]} '
            f'{mixed_positions[i, 2]} {mixed_positions[i, 3]}\n'
        )

        input_file_2 += line

    # Add lattice information to input_file_0
    lattice_lines = lattice.splitlines(True)
    input_file_0 += lattice_lines

    nonmixed_str = ''.join(nonmixed_input)

    # Combine all lines
    all_lines = input_file_0 + input_file_2 + nonmixed_input


    # add mixed ions
    input_file_1 = []
    element_string = ' '.join(elements) + ' '

    for i in range(mix_count):
        input_file_1.append(f'add-mix mix{i+1} {element_string}\n')

    # combine lines
    input_file = input_file_1 + all_lines
    input_file_combined = ''.join(input_file)
    
    # write the file
    input_file_path = os.path.join(folder, f'{prefix}.in')
    try:
        with open(input_file_path, "w") as f:
            f.write(input_file_combined)
    except Exception as e:
        print(f"An error occurred: {e}")


# run the dft calc as slurm job under reservation
def run_jdftx(input_file, output_file):
    cwd = os.getcwd()
    slurm_script_path = os.path.join(cwd, f'slurm_{input_file}')

    with open(slurm_script_path, 'w') as file:
        file.write(f"""#!/bin/bash
#SBATCH -N 1               # Request 1 node
#SBATCH -n 9 # Request 1 tasks (processes) on the node
#SBATCH -c 10              # Request 1 CPUs per task
#SBATCH -p hi_mem4         # Use the hi_mem partition
#SBATCH -J rel12+1_6           # Job name
#SBATCH -o slurm_out.o%j   # Output file
#SBATCH --time=999:00:00   # Max run time
module load mpi/openmpi-x86_64
module load intel-mkl
mpirun --bind-to none -n 9 /data2/jt577/jdftx_eat_withLibXC/build/jdftx -i {cwd}/{input_file} | tee {cwd}/{output_file}""")

    # Submit the job and get the job ID
    job_id = submit_job(slurm_script_path)

    # Return the job_id to wait for the job in a different function
    return job_id




def energy_surface(w, elements, folder, generic_inputs, lattice, mixed_positions, nonmixed_input, magmom, common_in):
    
    if common_in is not None:
        prefix = 'bandstruct'
    else:
        prefix = 'run'
    
    # determine input weights
    num_species = len(elements)
    weights = []

    # write input file and run jdftx
    write_input_surface(weights, elements, folder, generic_inputs, lattice, mixed_positions, nonmixed_input, magmom, common_in)
    # Save the current working directory
    original_cwd = os.getcwd()
    try:
        # Change the working directory to the subfolder
        os.chdir(folder)
        job_id = run_jdftx(f'{prefix}.in',f'{prefix}.out')
    finally:
        # Change back to the original directory
        os.chdir(original_cwd)

    # return job_id so can wait for job in a different function
    return job_id


# make folders, move files, run energy calculations. Must do before reading from files when computing cost or grad_cost
def perform_calc(w, elements, generic_inputs, lattice, mixed_positions, nonmixed_input, magmom, common_in):
    base_path = os.getcwd()
    # create folder to store runs if it doesn't exist
    run_folder = cdm.create_subfolder(folder_path = base_path, subfolder_name = 'runs')

    # create new unique folder
    unique_folder = cdm.create_unique_folder(base_path=run_folder, prefix='calc_')


    # do jdftx calcs
    job_id = energy_surface(w, elements, unique_folder, generic_inputs, lattice, mixed_positions, nonmixed_input, magmom, None)

    time.sleep(1)

    # return unique_folder so can read from files within it
    return unique_folder


# slurm running and waiting functions
def submit_job(script_path):
    """Submits a job to SLURM and returns the job ID."""
    result = subprocess.run(['sbatch', script_path], stdout=subprocess.PIPE, universal_newlines=True, check=True)
    # Output is something like 'Submitted batch job 12345'
    job_id = result.stdout.strip().split()[-1]
    return job_id