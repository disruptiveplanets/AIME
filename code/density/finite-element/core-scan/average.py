import os, sys
import numpy as np
import warnings
sys.path.append("../..")
from display import make_gif, make_slices
from core import TrueShape

NAME_START = "den-core-sph"
NUM_DRAWS = 20
NUM_CHOOSE = 1000
ASTEROID_NAME = "avg-core"
MAX_RADIUS = 2000
DURATION = 5
GENERATE = True
DIVISION = 49
PULL = True

grid_line = np.arange(-MAX_RADIUS, MAX_RADIUS, DIVISION)

if PULL:
    os.system("scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/density/finite-element/core-scan/*.npy .")

if GENERATE:
    density_grid_sum = 0
    density_grid_square_sum = 0
    num_grids = 0

    for i in range(NUM_DRAWS):
        if not os.path.exists(f"{NAME_START}-{i}-fe.npy"):
            # Run did not successfully complete
            continue
        with open(f"{NAME_START}-{i}-grids.npy", 'rb') as f:
            grids = np.load(f)
        with open(f"{NAME_START}-{i}-fe.npy", 'rb') as f:
            long_samples = np.load(f)
        densities = long_samples[np.random.randint(0, len(long_samples), NUM_CHOOSE),:]
        for d in densities:
            density_grid = np.einsum("ijkl,i->jkl", grids, d)
            density_grid_sum += density_grid
            density_grid_square_sum += density_grid**2
            num_grids += 1

    print(num_grids, "grids were averaged over")
    mean_grid = density_grid_sum / num_grids

    # Insert nans
    indicator = np.sum(grids, axis=0)>0
    mean_grid[~indicator] = np.nan

    unc_grid = np.sqrt(density_grid_square_sum / num_grids - mean_grid**2)
    

    with open("avg-grid.npy", 'wb') as f:
        np.save(f, mean_grid)
    with open("unc-grid.npy", 'wb') as f:
        np.save(f, unc_grid)
else:
    with open("avg-grid.npy", 'rb') as f:
        mean_grid = np.load(f)
    with open("unc-grid.npy", 'rb') as f:
        unc_grid = np.load(f)

x,y,z = np.meshgrid(grid_line, grid_line, grid_line)
true_densities = TrueShape.core_sph(1.5, 500)(x,y,z)
true_densities = true_densities.astype(float)
true_densities[np.isnan(mean_grid)] = np.nan

unc_grid /= np.nanmean(mean_grid)
mean_grid /= np.nanmean(mean_grid)
true_densities /= np.nanmean(true_densities)


if not os.path.isdir(f"../../figs/{ASTEROID_NAME}"):
    os.mkdir(f"../../figs/{ASTEROID_NAME}")

warnings.filterwarnings("ignore")

difference = mean_grid - true_densities
ratios = difference / unc_grid
uncertainty_ratios = unc_grid / mean_grid
num_elements = np.sum(~np.isnan(mean_grid))
error = np.nansum(ratios**2) / num_elements # Reduced chi squared

print("Plotting density")
make_slices(mean_grid, grid_line, "$\\rho$", 'plasma', f"../../figs/{ASTEROID_NAME}/fe-d", error)
make_gif(mean_grid, grid_line, "$\\rho$", 'plasma', f"../../figs/{ASTEROID_NAME}/fe-d.gif", DURATION)

print("Plotting ratios")
make_slices(ratios, grid_line, "$\\Delta \\rho / \\sigma_\\rho$", 'coolwarm', f"../../figs/{ASTEROID_NAME}/fe-r", error, percentile=95, balance=True)
make_gif(ratios, grid_line, "$\\Delta \\rho / \\sigma_\\rho$", 'coolwarm', f"../../figs/{ASTEROID_NAME}/fe-r.gif", duration=DURATION, percentile=95, balance=True)

print("Plotting uncertainty")
make_slices(uncertainty_ratios, grid_line, "$\\sigma_\\rho / \\rho$", 'Greys_r', f"../../figs/{ASTEROID_NAME}/fe-u", error, 95)
make_gif(uncertainty_ratios, grid_line, "$\\sigma_\\rho / \\rho$", 'Greys_r', f"../../figs/{ASTEROID_NAME}/fe-u.gif", DURATION, 95)

if true_densities is not None:
    print("Plotting differences")
    make_slices(difference, grid_line, "$\\Delta\\rho$", 'PuOr_r', f"../../figs/{ASTEROID_NAME}/fe-s", error, 95, balance=True)
    make_gif(difference, grid_line, "$\\Delta\\rho$", 'PuOr_r', f"../../figs/{ASTEROID_NAME}/fe-s.gif", DURATION, 95, balance=True)

warnings.filterwarnings("default")