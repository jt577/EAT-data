#####################################################################
# Functions that deal with creating, deleting, or moving folders
#####################################################################
# imports
import os
import shutil
import time

# create a subfolder if it doesn't exist already
def create_subfolder(folder_path, subfolder_name):
    subfolder_path = os.path.join(folder_path, subfolder_name)
    try:
        os.makedirs(subfolder_path, exist_ok=True)  # exist_ok=True to avoid error if the directory already exists
    except PermissionError as e:
        raise PermissionError(f"Permission denied when creating {subfolder_path}") from e
    
    return subfolder_path

# creates unique folder based on current timestamp
def create_unique_folder(base_path, prefix="calc_"):
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    folder_name = f"{prefix}"
    unique_folder_path = os.path.join(base_path, folder_name)
    return unique_folder_path

# removes folder
def remove_folder(folder):
    # Check if the directory exists and then delete it
    if os.path.exists(folder) and os.path.isdir(folder):
        shutil.rmtree(folder)
        print(f"Directory '{folder}' has been deleted.")
    else:
        print(f"Directory '{folder}' does not exist, so it cannot be deleted.")

# removes file
def remove_file(file):
    # Check if the file exists and then delete it
    if os.path.exists(file):
        os.remove(file)
        print(f"File '{file}' has been deleted.")
    else:
        print(f"File '{file}' does not exist, so it cannot be deleted.")

# delete progress files
def delete_progress(progress_file_name, details_file_name):
    # Remove progress files before restarting minimization
    remove_file(progress_file_name)
    remove_file(details_file_name)
    remove_file(f'O_pos.txt')
    remove_file(f'OH_pos.txt')
    remove_folder('runs')

# move wavefunctions from wavefn folder to unique folder
def mv_wfns_to_unique(source_folder, destination_folders):   
    # Define the wavefunction file paths
    wf_files = [
        os.path.join(source_folder, f'O.wfns'),
        os.path.join(source_folder, f'OH.wfns')
    ]

    # Move each wavefunction file to the destination folder
    numit = 0
    for file_path in wf_files:
        if os.path.exists(file_path):
            # Construct the destination path
            dest_path = os.path.join(destination_folders[numit], os.path.basename(file_path))
            # move the file
            shutil.move(file_path, dest_path)
            print(f"File '{file_path}' has been copied to '{dest_path}'.")
            numit += 1
        else:
            print(f"File '{file_path}' does not exist and cannot be copied.")
            numit += 1

# move wavefns from unique folder to wavefn folder
def mv_wfns_from_unique(source_folder, destination_folder):   
    # Ensure the destination folder exists, if not, create it
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
        print(f"Destination folder '{destination_folder}' has been created.")
    
    # Define the wavefunction file paths
    wf_files = [
        os.path.join(source_folder, f'O', f'O.wfns'),
        os.path.join(source_folder, f'OH', f'OH.wfns')
    ]

    # Move each wavefunction file to the destination folder
    for file_path in wf_files:
        if os.path.exists(file_path):
            # Construct the destination path
            dest_path = os.path.join(destination_folder, os.path.basename(file_path))
            # copy the file
            shutil.move(file_path, dest_path)
            print(f"File '{file_path}' has been copied to '{dest_path}'.")
        else:
            print(f"File '{file_path}' does not exist and cannot be copied.")


# move position files from unique timestamped folder to positions folder
def mv_pos_surface(source_folder, destination_folder):
    # Ensure the destination folder exists, if not, create it
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
        print(f"Destination folder '{destination_folder}' has been created.")
    
    # Define the position file paths
    pos_files = [
        os.path.join(source_folder, f'O', f'O.ionpos'), 
        os.path.join(source_folder, f'O', f'O.lattice'), 
        os.path.join(source_folder, f'OH', f'OH.ionpos'), 
        os.path.join(source_folder, f'OH', f'OH.lattice')
    ]

    # copy each position file to the destination folder
    for file_path in pos_files:
        if os.path.exists(file_path):
            # Construct the destination path
            dest_path = os.path.join(destination_folder, os.path.basename(file_path))
            # Move the file
            shutil.copy(file_path, dest_path)
            print(f"File '{file_path}' has been copied to '{dest_path}'.")
        else:
            print(f"File '{file_path}' does not exist and cannot be copied.")

# update lattice and ion position documents
def update_pos(iteration_counter):
    # Read ion positions from the file, skipping the first line
    ionpos_path = os.path.join('runs', 'positions', f'O.ionpos')
    with open(ionpos_path, 'r') as pos_file:
        pos_lines = pos_file.readlines()[1:]  # Skip the first line
    # Read lattice information from the file, reformatting each line except the first one
    latt_lines = []
    latt_path = os.path.join('runs', 'positions', f'O.lattice')
    with open(latt_path, 'r') as latt_file:
        for line_number, line in enumerate(latt_file):  # Use enumerate for line numbers
            if line_number != 0:  # Skip the first line
                parts = line.strip().split()
                if len(parts) >= 3:  # Ensure line has enough parts to avoid index errors
                    latt_lines.append(f'{parts[0]} {parts[1]} {parts[2]}')
    # Write to the destination file
    with open('O_pos.txt', 'a') as dest_file:
        dest_file.write(f'Iteration {iteration_counter}\n')
        dest_file.write('Lattice\n')
        dest_file.write('\n'.join(latt_lines) + '\n')  # Join the list into a single string
        dest_file.write('Positions\n')
        dest_file.writelines(pos_lines)  # Write each position line
        dest_file.write('\n')

    # Read ion positions from the file, skipping the first line
    ionpos_path = os.path.join('runs', 'positions', f'OH.ionpos')
    with open(ionpos_path, 'r') as pos_file:
        pos_lines = pos_file.readlines()[1:]  # Skip the first line
    # Read lattice information from the file, reformatting each line except the first one
    latt_lines = []
    latt_path = os.path.join('runs', 'positions', f'OH.lattice')
    with open(latt_path, 'r') as latt_file:
        for line_number, line in enumerate(latt_file):  # Use enumerate for line numbers
            if line_number != 0:  # Skip the first line
                parts = line.strip().split()
                if len(parts) >= 3:  # Ensure line has enough parts to avoid index errors
                    latt_lines.append(f'{parts[0]} {parts[1]} {parts[2]}')
    # Write to the destination file
    with open(f'OH_pos.txt', 'a') as dest_file:
        dest_file.write(f'Iteration {iteration_counter}\n')
        dest_file.write('Lattice\n')
        dest_file.write('\n'.join(latt_lines) + '\n')  # Join the list into a single string
        dest_file.write('Positions\n')
        dest_file.writelines(pos_lines)  # Write each position line
        dest_file.write('\n')