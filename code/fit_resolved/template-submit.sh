#!/bin/bash

# Slurm sbatsh options
#SBATCH -o log-fit-%j.log
#SBATCH --exclusive
#ggSBATCH --gres=gpu:volta:1

source /etc/profile

module load anaconda/2021a

python fit.py 
