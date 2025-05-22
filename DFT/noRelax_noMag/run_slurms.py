import os

# Get the current working directory
cwd = os.getcwd()

# Iterate over each folder in the current working directory
for folder in os.listdir(cwd):
    # Check if the item is a directory
    if os.path.isdir(os.path.join(cwd, folder)):
        # Change the working directory to the current folder
        os.chdir(os.path.join(cwd, folder))
        
        # Run sbatch slurm.sh
        os.system('sbatch slurm_runfile.sh')
        
        # Change the working directory back to the original directory
        os.chdir(cwd)