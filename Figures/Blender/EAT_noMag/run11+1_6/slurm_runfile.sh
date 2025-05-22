#!/bin/bash
#SBATCH -N 1               # number of nodes
#SBATCH -n 9 # Request 1 tasks (processes) on the node
#SBATCH -c 4              # Request 1 CPUs per task
#SBATCH -p hi_mem2         # Use the hi_mem partition
#SBATCH -J r11+1_6           # Job name
#SBATCH -o slurm_out.o%j   # Output file
#SBATCH --time=999:00:00   # Max run time

# Activate virtual python environment
source /data2/jt577/COSMO-RS/myenv/bin/activate

# Run optimization script
python main.py