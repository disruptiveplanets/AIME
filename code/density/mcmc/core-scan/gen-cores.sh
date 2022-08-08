#!/bin/bash

#SBATCH -ogen-cores-%a.log
#SBATCH -a 0-19
#SBATCH -c 8

source /etc/profile

module load anaconda/2022a

python gen-cores.py ori $SLURM_ARRAY_TASK_ID
python gen-cores.py vel $SLURM_ARRAY_TASK_ID
#python gen-cores.py sym $SLURM_ARRAY_TASK_ID
#python gen-cores.py double $SLURM_ARRAY_TASK_ID
#python gen-cores.py asym $SLURM_ARRAY_TASK_ID
#python gen-cores.py sph-3 $SLURM_ARRAY_TASK_ID
#python gen-cores.py sph-1.5 $SLURM_ARRAY_TASK_ID
#python gen-cores.py move-3 $SLURM_ARRAY_TASK_ID
#python gen-cores.py move-1.5 $SLURM_ARRAY_TASK_ID