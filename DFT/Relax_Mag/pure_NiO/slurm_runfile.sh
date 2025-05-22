#!/bin/bash
#SBATCH -N 1               # number of nodes
#SBATCH -n 9 # Request 1 tasks (processes) on the node
#SBATCH -c 10              # Request 1 CPUs per task
#SBATCH -p hi_mem4         # Use the hi_mem partition
#SBATCH -J NiO           # Job name
#SBATCH -o slurm_out.o%j   # Output file
#SBATCH --time=0   # Max run time

# Activate virtual python environment
source /data2/jt577/COSMO-RS/myenv/bin/activate

# Run optimization script
python main.py