#!/bin/bash

# Slurm sbatsh options
#SBATCH -odof-scan-%a.log
#SBATCH -a 0-24
#SBATCH -c 16

source /etc/profile

module load anaconda/2022a

python -u dof-scan.py probe-s-rho $SLURM_ARRAY_TASK_ID
