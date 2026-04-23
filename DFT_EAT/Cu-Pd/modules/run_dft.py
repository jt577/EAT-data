import subprocess
import os
import numpy as np
import shutil
import time
import modules.folder_ops as folder_ops

class run_dft:
    def __init__(self, out_dir, run_name, dft_names, dft_infiles, dft_parallel, stoichiometry_matrix):
        self.out_dir = out_dir
        self.dft_parallel = dft_parallel
        self.stoichiometry_matrix = stoichiometry_matrix

        # Setup folder structure
        self.base_path = self.out_dir
        self.run_folder = folder_ops.create_subfolder(folder_path=self.base_path, subfolder_name='runs')
        self.unique_folder = folder_ops.create_subfolder(folder_path=self.run_folder, subfolder_name=run_name)
        self.dft_paths = []

        # Make a pure dft calc for each EAT atom in stoichiometry matrix
        n_atoms = self.stoichiometry_matrix.shape[0]
        self.dft_names = []
        self.dft_infiles = []
        self.dft_names.append(dft_names[0]) # alloy is first
        self.dft_infiles.append(dft_infiles[0]) # alloy is first
        for i in range(n_atoms):
            self.dft_names.append(f"{dft_names[-1]}_{i}") # pure is last
            self.dft_infiles.append(dft_infiles[-1]) # Use pure infile as template
        for dft_name in self.dft_names:
            self.dft_paths.append(folder_ops.create_subfolder(folder_path=self.unique_folder, subfolder_name=dft_name))
            

    def _write_input(self, dft_infile, dft_name, dft_path, dirUpdateScheme, linminMethod, is_alloy, pure_index):
        """
        Writes the DFT input file with the provided stoichiometry matrix.
        """
        # Read template without trailing newlines to avoid double-spacing
        with open(dft_infile, 'r') as file:
            lines = [ln.rstrip('\n') for ln in file]

        if is_alloy:
            stoich = self.stoichiometry_matrix.numpy()
            n_mix = stoich.shape[0]

            new_lines = []
            for line in lines:
                newline = line
                for i in range(n_mix):
                    split_line = line.strip().split()
                    if split_line:
                        if split_line[0] == 'ion' and split_line[1] == f'mix{i}':
                            weights = ' '.join(map(str, stoich[i, :]))
                            parts = line.split()
                            newline = ' '.join(parts[:2] + [weights] + parts[2:])
                            break
                if "electronic-minimize" in line:
                    split_line = line.strip().split()
                    split_line.append(f'dirUpdateScheme {dirUpdateScheme} linminMethod {linminMethod}')
                    newline = ' '.join(split_line)
                if "dump-name" in line:
                    split_line = line.strip().split()
                    split_line[1] = f"{dft_name}.$VAR"
                    newline = ' '.join(split_line)
                new_lines.append(newline)
        else:
            # Pure EAT atom case
            stoich_mtx = self.stoichiometry_matrix.numpy()
            stoich_vec = stoich_mtx[pure_index, :]
            new_lines = []
            for line in lines:
                newline = line
                tag = f'ion mix'
                if tag in line:
                    weights = ' '.join(map(str, stoich_vec))
                    parts = line.split()
                    newline = ' '.join(parts[:2] + [weights] + parts[2:])
                if "electronic-minimize" in line:
                    split_line = line.strip().split()
                    split_line.append(f'dirUpdateScheme {dirUpdateScheme} linminMethod {linminMethod}')
                    newline = ' '.join(split_line)
                if "dump-name" in line:
                    split_line = line.strip().split()
                    split_line[1] = f"{dft_name}.$VAR"
                    newline = ' '.join(split_line)
                new_lines.append(newline)

        # Write the modified file
        input_file_path = os.path.join(dft_path, f'{dft_name}.in')
        try:
            with open(input_file_path, "w") as f:
                f.write('\n'.join(new_lines) + '\n')
        except Exception as e:
            print(f"An error occurred: {e}")
    
    def _run_one(self, dft_path, dft_name):
        """
        Runs one DFT calculation using jdftx with MPI parallelization.
        """
        os.chdir(dft_path)
        try:
            command = (
                f"mpirun --bind-to none -n {self.dft_parallel} "
                f" /data2/jt577/jdftx_eat_withLibXC/build/jdftx -i {dft_path}/{dft_name}.in | tee {dft_path}/{dft_name}.out"
            )
            subprocess.run(command, shell=True, check=True)
        finally:
            os.chdir(self.base_path)

    def run(self):
        """
        Runs all DFT calculations specified in dft_names list.
        """
        # Multiple fallback schemes for electronic minimization to reduce probability of non-convergence
        dirUpdateSchemes = ["PolakRibiere", "L-BFGS"]
        linminMethods = ["DirUpdateRecommended", "CubicWolfe"]
        self.dft_names_detailed = []
        counter = 0
        for dft_infile, dft_path, dft_name in zip(self.dft_infiles, self.dft_paths, self.dft_names):
            converged = False
            if counter == 0:
                is_alloy = True
            else:
                is_alloy = False
            for dirUpdateScheme in dirUpdateSchemes:
                if converged == True:
                    break
                for linminMethod in linminMethods:
                    dft_name_detailed = f"{dft_name}_{dirUpdateScheme}_{linminMethod}"
                    self._write_input(dft_infile, dft_name_detailed, dft_path, dirUpdateScheme, linminMethod, is_alloy, counter - 1)
                    self._run_one(dft_path, dft_name_detailed)
                    converged = self._check_convergence(dft_path, dft_name_detailed)
                    if converged == True:
                        self.dft_names_detailed.append(dft_name_detailed)
                        break
            # If all fallbacks fail
            if converged == False:
                self.dft_names_detailed.append(dft_name_detailed) # Append the last attempted name even if not converged
            counter += 1
        
    def _check_convergence(self, dft_path, dft_name):
        """
        Checks if the DFT calculation converged by looking for 'ElecMinimize: Converged' in the output file.
        If not found, it waits and checks again until convergence is achieved.
        """
        out_file = os.path.join(dft_path, f'{dft_name}.out')
        if not os.path.exists(out_file):
            raise FileNotFoundError(f"Output file {out_file} does not exist.")
        with open(out_file, 'r') as file:
            lines = file.readlines()
            if not any('ElecMinimize: Converged' in line for line in lines):
                return False
            else:
                return True

    def read(self):
        """
        Reads the energy and energy gradient from the DFT output files.
        """
        skip_run = False # Flag to skip reading energies and gradients if electronic minimization did not converge
        # Check to make sure electronic minimization converged
        for dft_path, dft_name_detailed in zip(self.dft_paths, self.dft_names_detailed):
            out_file = os.path.join(dft_path, f'{dft_name_detailed}.out')
            if not os.path.exists(out_file):
                raise FileNotFoundError(f"Output file {out_file} does not exist.")
            with open(out_file, 'r') as file:
                lines = file.readlines()
                if not any('ElecMinimize: Converged' in line for line in lines):
                    skip_run = True
        if skip_run is False:
            energies = {}
            energy_grads = {}
            counter = 0
            for dft_path, dft_name, dft_name_detailed in zip(self.dft_paths, self.dft_names, self.dft_names_detailed):
                if counter == 0:
                    is_alloy = True
                else:
                    is_alloy = False
                # read energies
                energy_file = os.path.join(dft_path, f'{dft_name_detailed}.Ecomponents')
                energy = None
                with open(energy_file, 'r') as file:
                    for line_number, line in enumerate(file, 1):
                        stripped_line = line.strip()
                        if 'F' in stripped_line:
                            split_line = stripped_line.split()
                            energy = split_line[2]
                            break
                if energy is None: # if F not found (happens for some reason when JDFTx does not output F and just Etot), skip this run
                    energies = None
                    energy_grads = None
                    break
                else:
                    energies[dft_name] = float(energy)
                    # read gradients
                    if is_alloy:
                        grad = np.zeros_like(self.stoichiometry_matrix)
                        grad_file = os.path.join(dft_path, f'{dft_name_detailed}.mixgrad')
                        for i in range(grad.shape[0]):
                            with open(grad_file, "r") as file:
                                for line_number, line in enumerate(file, 1):
                                    stripped_line = line.strip()
                                    if f'mix{i}' in stripped_line:
                                        split_line = stripped_line.split()
                                        for j in range(grad.shape[1]):
                                            grad[i,j] = float(split_line[2+j])
                                        break
                        energy_grads[dft_name] = grad.flatten()
                    else:
                        grad = np.zeros(self.stoichiometry_matrix.shape[1])
                        grad_file = os.path.join(dft_path, f'{dft_name_detailed}.mixgrad')
                        for i in range(grad.shape[0]):
                            with open(grad_file, "r") as file:
                                for line_number, line in enumerate(file, 1):
                                    stripped_line = line.strip()
                                    if f'mix{i}' in stripped_line:
                                        split_line = stripped_line.split()
                                        for j in range(len(grad)):
                                            grad[j] = float(split_line[2+j])
                                        break
                        energy_grads[dft_name] = grad.flatten()
                
                counter += 1
        else:
            energies = None
            energy_grads = None
        
        return energies, energy_grads