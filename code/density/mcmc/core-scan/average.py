import os, sys
from re import T
import numpy as np
import warnings
sys.path.append("..")
sys.path.append("../..")
from mcmc_core import MCMCAsteroid, log_like
from fe import FiniteElement
from display import make_gif, make_slices
from core import TrueShape, Indicator

RUN_NAME = "asym"
PULL = True
GENERATE = True

NUM_DRAWS = 20
NUM_CHOOSE = 1000
ASTEROID_NAME = f"avg-{RUN_NAME}"
MAX_RADIUS = 2000
DURATION = 5
DIVISION = 49
k22a, k20a = -0.05200629, -0.2021978
ELLIPSOID_AM = 1000


SURFACE_AMS = {
    "sph-3": 1000,
    "sph-1.5": 1000,
    "move-3": 1002.0081758422925,
    "move-1.5": 1000.1281468600504,
    "sym": 1000,
    "asym": 1000,
    "double": 1000,
}

BULK_AMS = {
    "sym": 1000,
    "asym": 1000,
    "sph-3": 922.9234884822591,
    "sph-1.5": 978.4541044108308,
    "move-3": 933.1648422811957,
    "move-1.5": 980.8811439828254,
    "double": 970.4652599064898
}

a = np.sqrt(5/3) * ELLIPSOID_AM * np.sqrt(1 - 2 * k20a + 12 * k22a)
b = np.sqrt(5/3) * ELLIPSOID_AM * np.sqrt(1 - 2 * k20a - 12 * k22a)
c = np.sqrt(5/3) * ELLIPSOID_AM * np.sqrt(1 + 4 * k20a)
core_displacement = 300
core_rad = 500
core_vol = np.pi * 4 / 3 * core_rad**3
ellipsoid_vol = np.pi * 4 / 3 * a * b * c
density_factor_low = 0.5
density_factor_high = 2
core_shift_low = core_displacement * (core_vol * density_factor_low) / ellipsoid_vol
core_shift_high = core_displacement * (core_vol * density_factor_high) / ellipsoid_vol

blob_rad = 300
core_one = np.array([0, 500, 0])
core_two = np.array([0, -500, 0])

TRUE_SHAPES = {
    "sym": TrueShape.uniform(),
    "asym": TrueShape.uniform(),
    "sph-3": TrueShape.core_sph(3, 500),
    "sph-1.5": TrueShape.core_sph(1.5, 500),
    "move-3": TrueShape.core_shift(3, 500, core_displacement),
    "move-1.5": TrueShape.core_shift(1.5, 500, core_displacement),
    "double": TrueShape.two_core(3, blob_rad, core_one, 3, blob_rad, core_two),

}
INDICATORS = {
    "sym": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "asym": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "sph-3": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "sph-1.5": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "move-3": Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_high),
    "move-1.5": Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_low),
    "double": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
}
SAMPLE_NAME = {
    "sym": "den-sym",
    "asym": "den-asym",
    "sph-3": "den-core-sph-3",
    "sph-1.5": "den-core-sph-1.5",
    "move-3": "den-core-move-3",
    "move-1.5": "den-core-move-1.5",
    "double": "den-core-double",
}

grid_line = np.arange(-MAX_RADIUS, MAX_RADIUS, DIVISION)

def check_moments(densities):
    zero_densities = densities.copy()
    zero_densities[np.isnan(densities)] = 0
    
    mcmc_asteroid = MCMCAsteroid("ast-test", f"../../samples/{SAMPLE_NAME[RUN_NAME]}-0-samples.npy", INDICATORS[RUN_NAME],
        TRUE_SHAPES[RUN_NAME], SURFACE_AMS[RUN_NAME], DIVISION, MAX_RADIUS, 9, BULK_AMS[RUN_NAME])
    moment_field = mcmc_asteroid.asteroid.moment_field(SURFACE_AMS[RUN_NAME])

    # Calculate klm
    unscaled_klm = np.einsum("iabc,abc->i", moment_field, zero_densities) * DIVISION**3
    radius_sqr = unscaled_klm[-1].real
    klms = unscaled_klm / radius_sqr
    klms[0] *= radius_sqr

    k33 = mcmc_asteroid.data_storage.data[2] + mcmc_asteroid.data_storage.data[3] * 1j
    k32 = mcmc_asteroid.data_storage.data[4] + mcmc_asteroid.data_storage.data[5] * 1j
    k31 = mcmc_asteroid.data_storage.data[6] + mcmc_asteroid.data_storage.data[7] * 1j
    k30 = mcmc_asteroid.data_storage.data[8]
    complex_klms = [
        1, 
        0, 0, 0, 
        mcmc_asteroid.data_storage.data[0], 0, mcmc_asteroid.data_storage.data[1], 0, mcmc_asteroid.data_storage.data[0],
        -k33.conj(), k32.conj(), -k31.conj(), k30, k31, k32, k33
    ]

    for m, d in zip(klms, complex_klms):
        print(f"Got {m}\t\t Wanted {d}")
    print(f"Got {np.sqrt(radius_sqr)}\t\t Wanted {BULK_AMS[RUN_NAME]}")

    # Find likelihood
    free_real_klms = np.array([
        klms[4].real, # K22
        klms[6].real, # K20
        klms[15].real, # R K33
        klms[15].imag, # I K33
        klms[14].real, # R K32
        klms[14].imag, # I K32
        klms[13].real, # R K31
        klms[13].imag, # I K31
        klms[12].real, # K30
    ])

    likelihood = log_like(free_real_klms, mcmc_asteroid.data_storage)
    print(f"Total likelihood: {likelihood}")



if PULL:
    os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/density/mcmc/core-scan/den-core-{RUN_NAME}*.npy {RUN_NAME}")
    os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/density/mcmc/core-scan/den-core-{RUN_NAME}*.png {RUN_NAME}")

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
            density_grid /= np.sum(density_grid)
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

check_moments(mean_grid)

x,y,z = np.meshgrid(grid_line, grid_line, grid_line)
true_densities = TRUE_SHAPES[RUN_NAME](x,y,z)
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