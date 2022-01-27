#!/bin/bash

# Slurm sbatsh options
#SBATCH -o surface.log
#SBATCH -c 24

source /etc/profile

module load anaconda/2021a

python surface.py
