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
            os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/thresholds/fe/dof-scan/{SCAN_NAME}-{dof}-{trial_num}.npy {DIRECTORY}/")

fig_unc, ax_unc = plt.subplots()
fig_mean_devs, ax_mean_devs = plt.subplots()
fig_max_devs, ax_max_devs = plt.subplots()

ax_unc.set_xlabel("Observational uncertainty")
ax_unc.set_ylabel("Uncertainty / local density")

ax_mean_devs.set_xlabel("Observational uncertainty")
ax_mean_devs.set_ylabel("Mean |(density deviation) / (uncertainty)|")
ax_mean_devs.set_yscale("log")

ax_max_devs.set_xlabel("Observational uncertainty")
ax_max_devs.set_ylabel("Max |(density deviation) / (uncertainty)|")
ax_max_devs.set_yscale("log")

for dof_index, dof in enumerate(DOFS):
    all_mean_devs = []
    all_max_devs = []
    all_uncs = []
    for trial_num in range(NUM_TRIALS):
        fname = f"{SCAN_NAME}-{dof}-{trial_num}.npy"
        if not os.path.exists(DIRECTORY + "/" + fname):
            raise Exception(f"Directory {DIRECTORY + '/' + fname} did not exist")

        with open(DIRECTORY + "/" + fname, 'rb') as f:
            these_mean_devs, these_max_devs, these_uncs = np.load(f)
            all_mean_devs.append(these_mean_devs)
            all_max_devs.append(these_max_devs)
            all_uncs.append(these_uncs)
    
    all_mean_devs = np.array(all_mean_devs)
    all_max_devs = np.array(all_max_devs)
    all_uncs = np.array(all_uncs)
    if len(fname[:-4]) < 8:
        print(fname[:-4], end='\t\t')
    else:
        print(fname[:-4], end='\t')
    xs = np.arange(len(all_uncs[0]))

    ax_unc.scatter(xs, np.mean(all_uncs, axis=0), label=f"{dof} DOF", color=f"C{dof_index}", s=16)
    ax_unc.plot(xs, np.mean(all_uncs, axis=0), color=f"C{dof_index}")
    for specific_uncs in all_uncs:
        ax_unc.scatter(xs, specific_uncs, color=f"C{dof_index}", alpha = 0.3, s=8)

    ax_mean_devs.scatter(xs, np.mean(all_mean_devs, axis=0), label=f"{dof} DOF", color=f"C{dof_index}", s=16)
    ax_mean_devs.plot(xs, np.mean(all_mean_devs, axis=0), color=f"C{dof_index}")
    for specific in all_mean_devs:
        ax_mean_devs.scatter(xs, specific, color=f"C{dof_index}", alpha = 0.3, s=8)

    ax_max_devs.scatter(xs, np.mean(all_max_devs, axis=0), label=f"{dof} DOF", color=f"C{dof_index}", s=16)
    ax_max_devs.plot(xs, np.mean(all_max_devs, axis=0), color=f"C{dof_index}")
    for specific in all_max_devs:
        ax_max_devs.scatter(xs, specific, color=f"C{dof_index}", alpha = 0.3, s=8)

fig_unc.legend()
fig_unc.savefig("dof-unc.png")

fig_mean_devs.legend()
fig_mean_devs.savefig("dof-mean-devs.png")

fig_max_devs.legend()
fig_max_devs.savefig("dof-max-devs.png")
plt.show()