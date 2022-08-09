import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("../../code/density")
import display

# Copied from average, which wrote the data used by this program
DIVISION = 49
MAX_RADIUS = 2000
grid_line = np.arange(-MAX_RADIUS, MAX_RADIUS, DIVISION)
# End copy

fe_vel_den = np.load("../../code/density/mcmc/core-scan/vel/avg-grid.npy")
fe_vel_unc = np.load("../../code/density/mcmc/core-scan/vel/unc-grid.npy")
fe_vel_rat = (fe_vel_unc / fe_vel_den)
fe_vel_rat_finite = fe_vel_rat.reshape(-1)[np.isfinite(fe_vel_rat.reshape(-1))]
fe_vel_den_finite = fe_vel_den.reshape(-1)[np.isfinite(fe_vel_den.reshape(-1))]

fe_ori_den = np.load("../../code/density/mcmc/core-scan/ori/avg-grid.npy")
fe_ori_unc = np.load("../../code/density/mcmc/core-scan/ori/unc-grid.npy")
fe_ori_rat = (fe_ori_unc / fe_ori_den)
fe_ori_rat_finite = fe_ori_rat.reshape(-1)[np.isfinite(fe_ori_rat.reshape(-1))]
fe_ori_den_finite = fe_ori_den.reshape(-1)[np.isfinite(fe_ori_den.reshape(-1))]

def compare_hist():
    bins = np.linspace(0.1, 0.35, 40)
    plt.hist(fe_vel_rat_finite, bins, density=True, histtype='step', label='Velocity')
    plt.hist(fe_ori_rat_finite, bins, density=True, histtype='step', label='Orientation')
    plt.xlabel("$\sigma_\\rho / \\rho$")
    plt.legend()
    plt.tight_layout()
    plt.savefig("hist.png")
    print("Velocity avg uncertainty", np.mean(fe_vel_rat_finite))
    print("Orientation avg uncertainty", np.mean(fe_ori_rat_finite))


def compare_gif():
    data = fe_ori_rat / fe_vel_rat - 1
    display.make_gif(data, grid_line, "$\sigma_\mathrm{orientation} / \sigma_\mathrm{velocity}$ - 1", "PRGn", "unc-ratio.gif", 2, percentile=99, balance=True)
    print("Uncertainty correlation", np.corrcoef(fe_vel_rat_finite, fe_ori_rat_finite))


    data = fe_ori_den / fe_vel_den - 1
    display.make_gif(data, grid_line, "$\\rho_\mathrm{orientation} / \\rho_\mathrm{velocity}$ - 1", "RdBu", "den-ratio.gif", 2, percentile=99, balance=True)
    print("Density correlation", np.corrcoef(fe_vel_den_finite, fe_ori_den_finite))


if __name__ == "__main__":
    compare_gif()
    compare_hist()