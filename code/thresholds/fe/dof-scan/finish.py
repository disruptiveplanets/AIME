import numpy as np
import matplotlib.pyplot as plt
import sys, os
sys.path.append("../../../density")
from core import UncertaintyTracker

N_RUNS = 48
DOFS = [9, 7, 5, 3, 3]
N_TRIALS = 5

def average(dof, run):
    unc_tracker = UncertaintyTracker()
    for trial_num in range(N_TRIALS):
        fname = f"cast-{dof}-{trial_num}-rho-{run}-map.npy"
        with open(fname, 'rb') as f:
            this_map = np.load(f)
            unc_tracker.update(this_map)
    density_map, uncertainty_map = unc_tracker.generate()
    true_map = np.ones_like(density_map)
    true_map[np.isnan(density_map)] = np.nan

    deviation_map = density_map - true_map
    ratio_map = deviation_map / uncertainty_map
    maps = [
        deviation_map / density_map,
        uncertainty_map / density_map,
        ratio_map,
    ]
    return [(
            np.nanpercentile(m, 100 - (100 - 68.27) / 2),
            np.nanpercentile(m, 100 - (100 - 95.45) / 2),
            np.nanpercentile(m, 50),
            np.nanpercentile(m, (100 - 95.45) / 2),
            np.nanpercentile(m, (100 - 68.27) / 2),
        ) for m in maps]

def single_plot(dof):
    global MIN_RUN
    avgs = []
    for run in range(N_RUNS):
        avg = average(dof, run)
        if avg is None:
            MIN_RUN = max(run, MIN_RUN)
            continue
        else:
            avgs.append(avg)
    return avgs

def scan(true_map):
    fig, axs = plt.subplots(ncols=0, nrows=3)
    axs[2].set_xlabel("Sigma rho")
    axs[0].set_ylabel("Deviation")
    axs[1].set_ylabel("Uncertainty")
    axs[2].set_ylabel("Ratio")
    
    results = np.array([single_plot(dof) for dof in DOFS])
    xs = range(MIN_RUN, N_RUNS)
    
    for i in range(3):
        for j, dof in enumerate(DOFS):
            axs[i].fill_between(xs, results[j, :,i,0], results[j, :,i,4], alpha=0.3, color=f"C{j}")
            axs[i].fill_between(xs, results[j, :,i,1], results[j, :,i,3], alpha=0.5, color=f"C{j}")
            axs[i].plot(xs, results[j, :,i,2], color=f"C{j}")
    
    fig.tight_layout()
    fig.savefig("scan.png")
    #fig.savefig("scan.pdf")
