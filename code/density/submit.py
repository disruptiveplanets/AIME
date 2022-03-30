import sys, os

names = [
    "sym-sph", 
    "asym-sph", 
    "sym-ell", 
    "asym-ell", 
    "tet", 
    "db", 
    "high", 
    "in", 
    "out", 
    "blob"
]

def run_text(asteroid):
    return """#!/bin/bash

# Slurm sbatsh options
#SBATCH -o {0}.log
#SBATCH -c 16

source /etc/profile

module load anaconda/2022a

python run.py {0} harmonic
python run.py {0} likelihood
python run.py {0} lumpy
""".format(asteroid)

def run(asteroid):
    with open(f"{asteroid}.sh", 'w') as f:
        f.write(run_text(asteroid))
    os.system(f"sbatch {asteroid}.sh")


if len(sys.argv) > 1:
    print(f"Running {sys.argv[1]}")
    run(sys.argv[1])

else:
    print("Running all")
    for n in names:
        print(n)
        run(n)