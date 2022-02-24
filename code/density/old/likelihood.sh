#!/bin/bash

# Slurm sbatsh options
#SBATCH -o likelihood-6.log
#SBATCH -c 24

source /etc/profile

module load anaconda/2021a

python likelihood.py
