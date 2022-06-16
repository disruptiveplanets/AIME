import os, sys
import numpy as np
import warnings
sys.path.append("../..")
from display import make_gif, make_slices
from core import TrueShape

RUN_NAME = "sph-3"
NUM_DRAWS = 20
NUM_CHOOSE = 1000
ASTEROID_NAME = f"avg-{RUN_NAME}"
MAX_RADIUS = 2000
DURATION = 5
GENERATE = True
DIVISION = 49
PULL = False

SURFACE_AMS = {
    "sph-3": 1000,
    "sph-1.5": 1000,
    "move-3": 1002.0081758422925,
    "move-1.5": 1000.1281468600504,
}

BULK_AMS = {
    "sph-3": 922.9234884822591,
    "sph-1.5": 978.4541044108308,
    "move-3": 933.1648422811957,
    "move-1.5": 980.8811439828254,
}

grid_line = np.arange(-MAX_RADIUS, MAX_RADIUS, DIVISION)

def get_moments(densities):
    global asteroid
    zero_densities = densities.copy()
    zero_densities[np.isnan(zero_densities)] = 0
    import core
    k22, k20, surface_am = -0.05200629, -0.2021978, SURFACE_AMS[RUN_NAME] # For the shape
    bulk_am = BULK_AMS[RUN_NAME]
    asteroid = core.Asteroid("ast-test", f"../../samples/den-core-{RUN_NAME}-0-samples.npy", surface_am, DIVISION, MAX_RADIUS, core.Indicator.ell(surface_am, k22, k20), TrueShape.uniform(), bulk_am)
    moment_field = asteroid.moment_field(surface_am) * asteroid.indicator_map

    # Correct moment_field
    moment_field[0] *= (bulk_am / surface_am)**2
    moment_field[1:4] *= (bulk_am / surface_am)**1
    moment_field[9:16] *= (bulk_am / surface_am)**(-1)

    # Calculate klm
    unscaled_klm = np.einsum("iabc,abc->i", moment_field, zero_densities)
    radius_sqr = unscaled_klm[-1]
    klms = unscaled_klm / radius_sqr

    return klms


if PULL:
    os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/density/finite-element/core-scan/den-core-{RUN_NAME}*.npy {RUN_NAME}")
    os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/density/finite-element/core-scan/den-core-{RUN_NAME}*.png {RUN_NAME}")

if GENERATE:
    density_grid_sum = 0
    density_grid_square_sum = 0
    num_grids = 0

    # X, Y, Z = np.meshgrid(grid_line, grid_line, grid_line)
    # radius_map = (X**2 + Y**2 + Z**2).astype(float)
    # radius_map /= np.sum(radius_map)

    for i in range(NUM_DRAWS):
        if not os.path.exists(f"{RUN_NAME}/den-core-{RUN_NAME}-{i}-fe.npy"):
            # Run did not successfully complete
            continue
        with open(f"{RUN_NAME}/den-core-{RUN_NAME}-{i}-grids.npy", 'rb') as f:
            grids = np.load(f)
        with open(f"{RUN_NAME}/den-core-{RUN_NAME}-{i}-fe.npy", 'rb') as f:
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
    

    with open(f"{RUN_NAME}/avg-grid.npy", 'wb') as f:
        np.save(f, mean_grid)
    with open(f"{RUN_NAME}/unc-grid.npy", 'wb') as f:
        np.save(f, unc_grid)
else:
    with open(f"{RUN_NAME}/avg-grid.npy", 'rb') as f:
        mean_grid = np.load(f)
    with open(f"{RUN_NAME}/unc-grid.npy", 'rb') as f:
        unc_grid = np.load(f)

moments = get_moments(mean_grid)
for m, d in zip(moments, asteroid.data):
    print(f"Got {m}\t\t Wanted {d}")

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