#!/bin/bash

# Slurm sbatsh options
#SBATCH -o surface.log
#SBATCH -c 48

source /etc/profile

module load anaconda/2021a

python surface.py
