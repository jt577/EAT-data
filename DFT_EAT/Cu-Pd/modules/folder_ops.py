# Functions that deal with creating, deleting, or moving files or folders
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