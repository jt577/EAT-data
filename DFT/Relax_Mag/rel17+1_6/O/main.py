##################################################################################
# Main file for EAT line search
##################################################################################
import os
import numpy as np
import modules.bsruncalc as bs
import modules.egread as eg
import modules.minimize as mn
import modules.tqentropy as tq
import modules.cdmfolders as cdm
from collections import defaultdict


def generic_inputs(adsorbate, magmom):
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
    f"spintype z-spin\n"
    f"initial-magnetic-moments {magmom}\n"
    f"add-U Cr d 0.128623 V d 0.113923 Ni d 0.235196 Co d 0.121273\n"
    f"\n"
    f"dump-name {adsorbate}.$VAR\n"
    f"dump End Lattice IonicPositions Ecomponents\n"
    f"\n")

    return res

# Read in weights
with open('weights.txt') as weights_file:
    weights_lines = weights_file.readlines()
ions = []
for i, weights_line in enumerate(weights_lines):
    split_weights_line = weights_line.strip().split()
    w = np.array([float(split_weights_line[3]), float(split_weights_line[5]), float(split_weights_line[7]), float(split_weights_line[9])])
    elements = [split_weights_line[4], split_weights_line[6], split_weights_line[8], split_weights_line[10]]
    indmax = np.argmax(w)
    ions.append(elements[indmax])

# positions of atoms
positions_file = 'positions.txt'
with open(positions_file, 'r') as file:
    position_lines = file.readlines()

magmom_vals = [1, -1, 1, -1, 1, 1, -1, -1, -1 , 1 , -1 , 1 , -1 , -1 , 1 , 1]
# Iterate through each line and check for 'mix'
mixed_positions = []
nonmixed_input = []
mix_num = 0
magmom_dict = defaultdict(list)
for i, line in enumerate(position_lines):
    if 'Lattice' in line:
        lattice_lines = position_lines[i + 1:i + 4]
    if len(line.split()) > 1:
        if 'mix' in line:
            split_line = line.split()
            split_line[1] = ions[mix_num]
            if float(split_line[4]) > 0:
                split_line[5] = '1'
            line_new = ' '.join(split_line)
            line_new = line_new + '\n'
            nonmixed_input.append(line_new)
            # magnetic moments
            magmom_dict[ions[mix_num]].append(magmom_vals[mix_num])
            mix_num += 1
        elif line.split()[1] == 'O' or line.split()[1] == 'H':
            split_line = line.split()
            if float(split_line[4]) > 0:
                split_line[5] = '1'
            line_new = ' '.join(split_line)
            line_new = line_new + '\n'
            nonmixed_input.append(line_new)

# convert dict → "Co 1 1 -1 Fe -1 -1 1 …" string
magmom = ' '.join(
    f"{elem} " + ' '.join(str(m) for m in magmom_dict[elem])
    for elem in magmom_dict
)

# define lattice explicitly
lattice = (
f"lattice \\ \n"
f"{lattice_lines[0].strip()} \\ \n"
f"{lattice_lines[1].strip()} \\ \n"
f"{lattice_lines[2].strip()} \n\n")

cdm.create_subfolder('.', 'runs')

unique_folder = bs.perform_calc([], [], generic_inputs, lattice, mixed_positions, nonmixed_input, magmom, common_in=None)