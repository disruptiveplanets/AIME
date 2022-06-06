#!/bin/bash

# Slurm sbatsh options
#SBATCH -ounder-0.log
#SBATCH -c 16

source /etc/profile

module load anaconda/2022a

python -u compute-uncs-underdetermined.py probe-s-rho 0
