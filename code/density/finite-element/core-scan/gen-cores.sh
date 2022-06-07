#!/bin/bash

#SBATCH -ogen-cores-%a.log
#SBATCH -a 0-19
#SBATCH -c 4

source /etc/profile

module load anaconda/2022a

python gen-cores.py $SLURM_ARRAY_TASK_ID