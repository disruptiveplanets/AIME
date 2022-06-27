import os
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('jcap')

PULL = True
DIRECTORY = "fe/dof-scan"
DOFS = [9, 7, 5, 3, 2]
NUM_TRIALS = 5
PRECISE = False
if PRECISE:
    SCAN_NAME = "probe-s-rho"
else:
    SCAN_NAME = "cprobe-s-rho"

if PULL:
    for dof in DOFS:
        for trial_num in range(NUM_TRIALS):
            os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/thresholds/fe/{SCAN_NAME}-{dof}-{trial_num}.npy {DIRECTORY}/")

fig, axs = plt.subplots(ncols=2, nrows=2, sharex=True)

ax_dev = axs[0,0]
ax_unc = axs[0,1]
ax_mean_dev_rats = axs[0,0]
ax_max_dev_rats = axs[0,1]

ax_max_dev_rats.set_xlabel("$\sigma_P$")
ax_max_dev_rats.set_xlabel("$\sigma_P$")

ax_unc.set_ylabel("$\sigma_\\rho / \\rho$")

ax_dev.set_ylabel("$\Delta \\rho / \\rho$")

ax_mean_dev_rats.set_ylabel("Mean $|\Delta \\rho / \sigma_\\rho|$")
ax_mean_dev_rats.set_yscale("log")

ax_max_dev_rats.set_ylabel("Max $|\Delta \\rho / \sigma_\\rho|$")
ax_max_dev_rats.set_yscale("log")

for dof_index, dof in enumerate(DOFS):
    all_mean_dev_rats = []
    all_max_dev_rats = []
    all_uncs = []
    all_devs = []
    for trial_num in range(NUM_TRIALS):
        fname = f"{SCAN_NAME}-{dof}-{trial_num}.npy"
        if not os.path.exists(DIRECTORY + "/" + fname):
            these_mean_devs = np.ones(5) * np.nan
            these_max_devs = np.ones(5) * np.nan
            these_uncs = np.ones(5) * np.nan
            #raise Exception(f"Directory {DIRECTORY + '/' + fname} did not exist")
        else:
            with open(DIRECTORY + "/" + fname, 'rb') as f:
                these_devs, these_devs, these_uncs = np.load(f)
        
        all_mean_dev_rats.append(np.mean(np.abs(these_devs / these_devs)))
        all_max_dev_rats.append(np.max(np.abs(these_devs / these_devs)))
        all_devs.append(these_devs)
        all_uncs.append(these_uncs)
    
    all_mean_dev_rats = np.array(all_mean_dev_rats)
    all_max_dev_rats = np.array(all_max_dev_rats)
    all_devs = np.array(all_devs)
    all_uncs = np.array(all_uncs)

    if len(fname[:-4]) < 8:
        print(fname[:-4], end='\t\t')
    else:
        print(fname[:-4], end='\t')
    xs = np.arange(len(all_uncs[0]))

    ax_unc.scatter(xs, np.nanmean(all_uncs, axis=0), label=f"{dof} DOF", color=f"C{dof_index}", s=16)
    ax_unc.plot(xs, np.nanmean(all_uncs, axis=0), color=f"C{dof_index}")
    for specific_uncs in all_uncs:
        ax_unc.scatter(xs, specific_uncs, color=f"C{dof_index}", alpha = 0.3, s=8)

    ax_mean_dev_rats.scatter(xs, np.nanmean(all_mean_dev_rats, axis=0), color=f"C{dof_index}", s=16)
    ax_mean_dev_rats.plot(xs, np.nanmean(all_mean_dev_rats, axis=0), color=f"C{dof_index}")
    for specific in all_mean_dev_rats:
        ax_mean_dev_rats.scatter(xs, specific, color=f"C{dof_index}", alpha = 0.3, s=8)

    ax_max_dev_rats.scatter(xs, np.nanmean(all_max_dev_rats, axis=0), color=f"C{dof_index}", s=16)
    ax_max_dev_rats.plot(xs, np.nanmean(all_max_dev_rats, axis=0), color=f"C{dof_index}")
    for specific in all_max_dev_rats:
        ax_max_dev_rats.scatter(xs, specific, color=f"C{dof_index}", alpha = 0.3, s=8)

    ax_dev.scatter(xs, np.nanmean(all_devs, axis=0), color=f"C{dof_index}", s=16)
    ax_dev.plot(xs, np.nanmean(all_devs, axis=0), color=f"C{dof_index}")
    for specific in all_devs:
        ax_dev.scatter(xs, specific, color=f"C{dof_index}", alpha = 0.3, s=8)

fig.legend(ncol=len(DOFS))
fig.tight_layout()
if PRECISE:
    fig.savefig("scan-dof-precise.png")
else:
    fig.savefig("scan-dof.png")

plt.show()