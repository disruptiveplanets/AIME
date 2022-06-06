#!/bin/bash

# Slurm sbatsh options
#SBATCH -ounder-6.log
#SBATCH -c 16

source /etc/profile

module load anaconda/2022a

python compute-uncs-underdetermined.py probe-s-rho 6
