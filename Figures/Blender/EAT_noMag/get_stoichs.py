import os

# Get the current working directory
cwd = os.getcwd()

# Iterate over each folder in the current working directory
for folder in os.listdir(cwd):
    # Check if the item is a directory
    if os.path.isdir(os.path.join(cwd, folder)):
        # Change the working directory to the current folder
        os.chdir(os.path.join(cwd, folder))
        # Open min_progress.txt
        with open('min_progress.txt', 'r') as f:
            lines = f.readlines()
        final_weights_lines = lines[-18:-2]
        elementDict = {}
        elementDict['Ni'] = 0
        elementDict['Cr'] = 0
        elementDict['Co'] = 0
        elementDict['V'] = 0
        for line in final_weights_lines:
            # Split the line into parts
            parts = line.strip().split()
            elementDict['Ni'] += float(parts[3]) 
            elementDict['Cr'] += float(parts[5])
            elementDict['Co'] += float(parts[7])
            elementDict['V'] += float(parts[9])
        elementDict['Ni'] = round(elementDict['Ni']) / 16
        elementDict['Cr'] = round(elementDict['Cr']) / 16
        elementDict['Co'] = round(elementDict['Co']) / 16
        elementDict['V'] = round(elementDict['V']) / 16
        # Write the results to stoich.txt
        with open('stoich.txt', 'w') as f:
            f.write(f"Ni: {elementDict['Ni']}\n")
            f.write(f"Cr: {elementDict['Cr']}\n")
            f.write(f"Co: {elementDict['Co']}\n")
            f.write(f"V: {elementDict['V']}\n")
