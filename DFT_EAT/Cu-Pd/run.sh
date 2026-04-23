#!/bin/bash
#SBATCH -N 1 -n 96 -c 1 -p hi_mem4 -J EAT -o slurm_out.o%j
#SBATCH --time=999:00:00

source /data2/jt577/COSMO-RS/myenv/bin/activate

python /data2/jt577/2025/November_2025/Phase_Diagram/EAT/main.py /path/to/infile.txt /path/to/outfolder