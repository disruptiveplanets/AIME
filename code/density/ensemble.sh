#!/bin/bash

# Slurm sbatsh options
#SBATCH -o ensemble.log
#SBATCH -c 24

source /etc/profile

module load anaconda/2021a

python ensemble.py
