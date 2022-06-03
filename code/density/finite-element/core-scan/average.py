import os
import numpy as np
import warnings
os.path.append("../..")
from display import make_gif, make_slices

NAME_START = "avg-core"
NUM_DRAWS = 20
NUM_CHOOSE = 1000
ASTEROID_NAME = "avg-core"

grids = []

for i in range(NUM_DRAWS):
    with open(f"NAME_START-{i}-grids.npy", 'rb') as f:
        grids = np.load(f)
    with open(f"NAME_START-{i}-fe.npy", 'rb') as f:
        long_samples = np.load(f)
    
    densities = np.random.choose(long_samples, NUM_CHOOSE)
    for d in densities:
        grids.append(np.einsum("ijkl,i->jkl", grids, d))

mean_grid = np.mean(grids, axis=0)
unc_grid = np.std(grids, axis=0)

with open("avg-grid.npy", 'wb') as f:
    np.save(f, mean_grid)
with open("unc-grid.npy", 'wb') as f:
    np.save(f, unc_grid)

true_densities = 

unc_grid /= np.nanmean(mean_grid)
mean_grid /= np.nanmean(mean_grid)
true_densities /= np.nanmean(true_densities)


if not os.path.isdir(f"../figs/{ASTEROID_NAME}"):
    os.mkdir(f"../figs/{ASTEROID_NAME}")

warnings.filterwarnings("ignore")

ratios = (mean_grid - true_densities) / (unc_grid)
make_slices(ratios, grid_line, "$\\Delta\\sigma$", 'coolwarm', f"../figs/{ASTEROID_NAME}/fe-r", error, percentile=95, balance=True)
make_gif(ratios, grid_line, "$\\Delta\\sigma$", 'coolwarm', f"../figs/{ASTEROID_NAME}/fe-r.gif", duration=duration, percentile=95, balance=True)
difference = (true_densities - densities)

print("Plotting density")
make_slices(densities, grid_line, "$\\rho$", 'plasma', f"../figs/{ASTEROID_NAME}/fe-d", error)
make_gif(densities, grid_line, "$\\rho$", 'plasma', f"../figs/{ASTEROID_NAME}/fe-d.gif", duration)

print("Plotting uncertainty")
make_slices(uncertainty_ratios, grid_line, "$\\sigma_\\rho / \\rho$", 'Greys_r', f"../figs/{ASTEROID_NAME}/fe-u", error, 95)
make_gif(uncertainty_ratios, grid_line, "$\\sigma_\\rho / \\rho$", 'Greys_r', f"../figs/{ASTEROID_NAME}/fe-u.gif", duration, 95)

if true_densities is not None:
    print("Plotting differences")
    make_slices(difference, grid_line, "$\\Delta\\rho$", 'PuOr_r', f"../figs/{ASTEROID_NAME}/fe-s", error, 95, balance=True)
    make_gif(difference, grid_line, "$\\Delta\\rho$", 'PuOr_r', f"../figs/{ASTEROID_NAME}/fe-s.gif", duration, 95, balance=True)

warnings.filterwarnings("default")