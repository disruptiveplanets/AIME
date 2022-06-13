import os
import numpy as np
import matplotlib.pyplot as plt

PULL = False
DIRECTORY = "fe/dof-scan"
DOFS = [9, 7, 5, 3, 2]
NUM_TRIALS = 5
SCAN_NAME = "cprobe-s-rho"

if PULL:
    for dof in DOFS:
        for trial_num in range(NUM_TRIALS):
            os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/thresholds/fe/{SCAN_NAME}-{dof}-{trial_num}.npy {DIRECTORY}/")

plt.figure()
plt.xlabel("Observational precision")
plt.ylabel("Uncertainty / local density")

for dof_index, dof in enumerate(DOFS):
    all_uncs = []
    for trial_num in range(NUM_TRIALS):
        fname = f"{SCAN_NAME}-{dof}-{trial_num}.npy"
        if not os.path.exists(DIRECTORY + "/" + fname):
            raise Exception(f"Directory {DIRECTORY + '/' + fname} did not exist")

        with open(DIRECTORY + "/" + fname, 'rb') as f:
            these_uncs = np.load(f)
            all_uncs.append(these_uncs)
    
    all_uncs = np.array(all_uncs)
    if len(fname[:-4]) < 8:
        print(fname[:-4], end='\t\t')
    else:
        print(fname[:-4], end='\t')
    xs = np.arange(len(all_uncs[0]))
    plt.scatter(xs, np.mean(all_uncs, axis=0), label=f"{dof} DOF", color=f"C{dof_index}", s=16)
    plt.plot(xs, np.mean(all_uncs, axis=0), color=f"C{dof_index}")
    for specific_uncs in all_uncs:
        plt.scatter(xs, specific_uncs, color=f"C{dof_index}", alpha = 0.3, s=8)

plt.legend()
plt.savefig("vary-num.png")
plt.show()