import emcee, corner, sys, random, os, warnings
import numpy as np
from scipy.linalg import pinvh, inv
from matplotlib import pyplot as plt
from multiprocessing import Pool, Lock
import grids
sys.path.append("..")
from core import Asteroid, Indicator, TrueShape
from scipy.optimize import minimize
from threading import Thread
from display import make_gif, make_slices

N_WALKERS = 32
MAX_L = 3
N_FITTED_MOMENTS = 9
N_CONSTRAINED = 7
N_FREE = N_FITTED_MOMENTS
N_ALL = N_FREE + N_CONSTRAINED
MAX_N_STEPS = 100_000
DIVISION = 49
RLM_EPSILON = 1e-20
MAX_RADIUS = 2000
SMALL_SIGMA = 1e-5
CERTAIN_INDICES = [0, 1, 2, 3, 5, 7, -1]
MAX_DENSITY = 3 # Iron
MIN_DENSITY = 0.25
UNCERTAINTY_RATIO = 0.25
VERY_LARGE_SLOPE = 1e30 # Can be infinite for mcmc
NUM_THREADS = os.cpu_count()
MIN_LOG_LIKE = 1000# 100
MINIMIZATION_ATTEMPTS = 500

def get_cov(path):
    with open(path, 'rb') as f:
        flat_samples = np.load(f).reshape(-1, N_FITTED_MOMENTS + 1)[:, 1:]
    if flat_samples.shape[0] == 0:
        return None, None
    cov = np.cov(flat_samples.transpose())
    data = np.mean(flat_samples, axis=0)
    return data, cov

def get_theta_long(theta_short, info):
    # Get the densities consistent with making the mass 1 and com 0 and rotation
    return np.append(theta_short, info.rlm_fixed_inv[:,-1] - info.rlm_prod @ theta_short)

def get_klms(densities, info):
    unscaled_klms = np.einsum('iabc,abc->i', info.rlm_mat, densities)
    radius_sqr = np.sum(info.radius_vec * densities)
    #radius_sqr = 1000**2
    scaled_klms = np.array(unscaled_klms) / radius_sqr
    scaled_klms[-1] *= radius_sqr # Do not scale mass term
    return scaled_klms


def log_prior(theta_long, info):
    mask_too_small = theta_long < info.mean_density * MIN_DENSITY
    if np.any(mask_too_small):
        return VERY_LARGE_SLOPE * np.sum((theta_long / info.mean_density - MIN_DENSITY) * mask_too_small)
    mask_too_big = theta_long > info.mean_density * MAX_DENSITY
    if np.any(mask_too_big):
        return -VERY_LARGE_SLOPE * np.sum((np.max(theta_long) / info.mean_density - MAX_DENSITY) * mask_too_big)
    return 0.0
    
def log_like(theta_long, info):
    diff_klms = get_klms(theta_long, info)[:N_FITTED_MOMENTS] - info.data # Only need the unconstrained ones
    print(get_klms(theta_long, info))
    print(diff_klms)
    return -0.5 * diff_klms.transpose() @ info.data_inv_covs @ diff_klms

def log_probability(theta_short, info):
    theta_long = get_theta_long(theta_short, info)
    lp = log_prior(theta_long, info)
    ll = log_like(theta_long, info)
    return ll + lp

class AsteroidInfo:
    def __init__(self, mean_density, data, data_inv_covs, rlm_mat, radius_vec, rlm_fixed_inv, rlm_prod, masks):
        self.mean_density = mean_density
        self.data = data
        self.data_inv_covs = data_inv_covs
        self.rlm_mat = rlm_mat
        self.radius_vec = radius_vec
        self.rlm_fixed_inv = rlm_fixed_inv
        self.rlm_prod = rlm_prod
        self.masks = masks

def load(name, asteroid, surface_am, sample_path, division, generate=True):
    mean_density = 1 / (np.sum(asteroid.indicator_map) * division**3)
    data, cov = get_cov(sample_path)
    if cov is None:
        # Sample path was empty
        return None
    data_inv_covs = pinvh(cov)
    if generate:
        masks = grids.get_grids_centroid(N_ALL, asteroid.grid_line, asteroid.indicator_map, asteroid.indicator)[1]
        with open(name + "-grids.npy", 'wb') as f:
            np.save(f, masks)
    else:
        with open(name + "-grids.npy", 'rb') as f:
            masks = np.load(f)
    
    rlms = asteroid.moment_field(surface_am) * division**3

    rlm_mat_complex = rlms * division**3

    # Set order
    rlm_mat = np.zeros_like(rlm_mat_complex, dtype=float)
    # Unconstrained
    rlm_mat[0, :] = rlm_mat_complex[8,:].real # K22
    rlm_mat[1, :] = rlm_mat_complex[6,:].real # K20
    rlm_mat[2, :] = rlm_mat_complex[15,:].real # R K33
    rlm_mat[3, :] = rlm_mat_complex[15,:].imag # I K33
    rlm_mat[4, :] = rlm_mat_complex[14,:].real # R K32
    rlm_mat[5, :] = rlm_mat_complex[14,:].imag # I K32
    rlm_mat[6, :] = rlm_mat_complex[13,:].real # R K31
    rlm_mat[7, :] = rlm_mat_complex[13,:].imag # I K31
    rlm_mat[8, :] = rlm_mat_complex[12,:].real # K30
    # Constrained
    rlm_mat[9, :] = rlm_mat_complex[ 3,:].real # R K11
    rlm_mat[10, :] = rlm_mat_complex[3,:].imag # I K11
    rlm_mat[11, :] = rlm_mat_complex[2,:].real # K10
    rlm_mat[12, :] = rlm_mat_complex[7,:].real # R K21
    rlm_mat[13, :] = rlm_mat_complex[7,:].imag # I K21
    rlm_mat[14, :] = rlm_mat_complex[8,:].imag # I K22
    rlm_mat[15, :] = rlm_mat_complex[0,:].real # K00
    radius_vec = rlm_mat_complex[-1, :].real # radius

    rlm_fixed_inv = None # inv(rlm_mat[np.arange(N_FREE, N_ALL),:][:,np.arange(N_FREE, N_ALL)])
    rlm_cross = None # rlm_mat[np.arange(N_FREE, N_ALL),:][:,np.arange(0, N_FREE)]
    rlm_prod = None # rlm_fixed_inv @ rlm_cross
    return AsteroidInfo(mean_density, data, data_inv_covs, rlm_mat, radius_vec, rlm_fixed_inv, rlm_prod, masks)

def get_map(info, means, unc, asteroid):
    densities = np.einsum("ijkl,i->jkl", info.masks, means)
    unc_ratios = np.einsum("ijkl,i->jkl", info.masks, unc)
    densities = densities / np.nanmean(densities)

    densities[~asteroid.indicator_map] = np.nan
    unc_ratios[~asteroid.indicator_map] = np.nan

    return densities, unc_ratios

def pipeline(name, sample_path, indicator, surface_am, division, max_radius, map, used_bulk_am=None, generate=True):
    if used_bulk_am is None:
        used_bulk_am = surface_am
    
    asteroid = Asteroid(name, sample_path, surface_am, division, max_radius, indicator, TrueShape.uniform(), used_bulk_am)
    asteroid_info = load(name, asteroid, surface_am, sample_path, division, generate)
    
    with open("core-scan/avg-grid.npy", 'rb') as f:
        means = np.load(f)
    with open("core-scan/unc-grid.npy", 'rb') as f:
        unc = np.load(f)

    means[np.isnan(means)] = 0

    print(log_like(means, asteroid_info))

    return unc
    
if __name__ == "__main__":
    k22, k20, surface_am = -0.05200629, -0.2021978, 1000 # For the shape
    pipeline(f"den-core-sph", "../samples/den-core-sph-0-samples.npy", Indicator.ell(surface_am, k22, k20),
        surface_am, DIVISION, MAX_RADIUS, False, used_bulk_am=978.4541044108308)