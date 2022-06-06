import emcee, corner, sys, random, os, warnings
import numpy as np
from scipy.linalg import pinvh, inv
from matplotlib import pyplot as plt
from multiprocessing import Pool, Lock
from grids import get_grids_centroid
sys.path.append("..")
from core import Asteroid, Indicator, TrueShape
from scipy.optimize import minimize
from threading import Thread
from display import make_gif, make_slices

N_WALKERS = 32
MAX_L = 3
N_ALL_DIM = (MAX_L + 1)**2
N_FREE_DIM = N_ALL_DIM - 7
MAX_N_STEPS = 100_000
DIVISION = 99
RLM_EPSILON = 1e-20
MAX_RADIUS = 2000
SMALL_SIGMA = 1e-5
CERTAIN_INDICES = [0, 1, 2, 3, 5, 7, -1]
MAX_DENSITY = 9 # Iron
MIN_DENSITY = 0.5 
UNCERTAINTY_RATIO = 0.25
VERY_LARGE_SLOPE = 1e30 # Can be infinite for mcmc
NUM_THREADS = os.cpu_count()
MIN_LOG_LIKE = 1000# 100
MINIMIZATION_ATTEMPTS = 500

def get_cov(path):
    with open(path, 'rb') as f:
        flat_samples = np.load(f).reshape(-1, N_FREE_DIM + 1)[:, 1:]
    if flat_samples.shape[0] == 0:
        return None, None
    cov = np.cov(flat_samples.transpose())
    data = np.mean(flat_samples, axis=0)
    return data, cov

def get_theta_long(theta_short, info):
    # Get the densities consistent with making the mass 1 and com 0 and rotation
    return np.append(theta_short, info.rlm_fixed_inv[:,-1] - info.rlm_prod @ theta_short)

def get_klms(theta_long, info):
    unscaled_klms = info.rlm_mat @ theta_long
    radius_sqr = info.radius_vec @ theta_long
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
    diff_klms = get_klms(theta_long, info)[:N_FREE_DIM] - info.data
    return -0.5 * diff_klms.transpose() @ info.data_inv_covs @ diff_klms

def log_probability(theta_short, info):
    theta_long = get_theta_long(theta_short, info)
    lp = log_prior(theta_long, info)
    ll = log_like(theta_long, info)
    return ll + lp

def mcmc_fit(theta_start, output_name, info, generate=True):
    if generate:
        backend = emcee.backends.HDFBackend(output_name+".h5")
        backend.reset(N_WALKERS, N_FREE_DIM)
        old_tau = np.inf

        with Pool() as pool:
            sampler = emcee.EnsembleSampler(N_WALKERS, N_FREE_DIM, log_probability, args = (info,), backend=backend, pool=pool)

            pos = np.zeros((N_WALKERS, N_FREE_DIM))
            for i in range(N_FREE_DIM):
                pos[:,i] = np.random.randn(N_WALKERS) * UNCERTAINTY_RATIO * info.mean_density + theta_start[i]
            

            for sample in sampler.sample(pos, iterations=MAX_N_STEPS, progress=True):
                if sampler.iteration % 500 == 0:
                    # Compute the autocorrelation time so far
                    # Using tol=0 means that we'll always get an estimate even
                    # if it isn't trustworthy
                    tau = sampler.get_autocorr_time(tol=0)

                    # Check convergence
                    converged = np.all(tau * 100 < sampler.iteration)
                    converged &= np.all(np.abs(old_tau - tau) / tau < 0.01)
                    if converged:
                        print("Converged")
                        break
                    old_tau = tau
                    print(np.mean(sampler.get_last_sample().log_prob), tau)

            #sampler._previous_state = sample

    else:
        backend = emcee.backends.HDFBackend(output_name+".h5", read_only=True)
        sampler = emcee.EnsembleSampler(N_WALKERS, N_FREE_DIM, log_probability, args=(info,), backend=backend)

    return sampler

class MinResult:
    def __init__(self):
        self.lock = Lock()
        self.val = None
        self.attempts = 0
    
    def set(self, v):
        self.lock.acquire()
        self.val = v
        self.lock.release()

    def increment(self):
        self.lock.acquire()
        self.attempts += 1
        self.lock.release()

    def query(self, threshold):
        self.lock.acquire()
        greater = self.attempts > threshold
        self.lock.release()
        return greater

    def get(self):
        self.lock.acquire()
        v = self.val
        self.lock.release()
        return v

    def is_set(self):
        self.lock.acquire()
        result = False if self.val is None else True
        self.lock.release()
        return result

def minimize_func(theta, info, result):
    if result.is_set() or result.query(MINIMIZATION_ATTEMPTS):
        raise Exception()
    return -log_probability(theta, info)

def min_func_mcmc(seed, info, result):
    local_rng = random.Random()
    local_rng.seed(seed)
    while not result.is_set():
        val = [(local_rng.random() * 4 + 0.5) * info.mean_density for _ in range(N_FREE_DIM)]
        try:
            min_result = minimize(minimize_func, x0=val, method="Nelder-Mead", args=(info, result,), options = {"maxiter": 500 * len(val)})
        except Exception:
            print("Thread bailed")
            return
        if min_result.success and min_result.fun < MIN_LOG_LIKE:
            print(f"Thread successfully completed with log like {min_result.fun}")
            result.set(min_result.x)
            return
        else:
            print(f"Attempt failed with log like {min_result.fun}")
            result.increment()

def get_theta_start_mcmc(info):
    threads = []
    result = MinResult()
    print(f"Starting {NUM_THREADS} threads")
    for i in range(NUM_THREADS):
        seed = random.randint(0, 0xffff_ffff_ffff_ffff)
        t = Thread(target=min_func_mcmc, args=(seed, info, result))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    result = result.get()
    #print(Hessian(lambda theta: -log_probability(theta, info))(result))
    return result
    
def get_densities_mcmc(output_name, info, generate=True):
    if generate:
        theta_start = get_theta_start_mcmc(info)
        print("Theta start:", theta_start)
        if theta_start is None:
            return np.nan, np.nan
        sampler = mcmc_fit(theta_start, output_name, info, True)

        tau = sampler.get_autocorr_time()
        max_tau = int(np.max(tau))
        print("Autocorrelation time", tau)
        flat_samples = sampler.get_chain(discard=2 * max_tau, thin=max_tau // 2, flat=True)


        long_samples = np.array([get_theta_long(theta_short, info) for theta_short in flat_samples])

        with open(output_name + ".npy", 'wb') as f:
            np.save(f, long_samples)
    
    else:
        with open(output_name + ".npy", 'rb') as f:
            long_samples = np.load(f)

    long_means, unc = get_stats_from_long_samples(long_samples)

    fig = corner.corner(long_samples / info.mean_density, truths=np.ones(N_ALL_DIM))
    corner.overplot_lines(fig, long_means / info.mean_density, color='C1')
    fig.savefig(output_name + ".png")

    return long_means, unc / long_means

def get_stats_from_long_samples(long_samples):
    short_samples = long_samples[:, :N_FREE_DIM]

    # Don't compute the mean; compute the middle of the data set.
    least_dist = None
    for i, point in enumerate(short_samples):
        mean_dist = np.sum((short_samples - point)**2) / len(short_samples)
        if least_dist is None or least_dist > mean_dist:
            least_dist = mean_dist
            least_point_index = i
    
    long_means = long_samples[least_point_index]
    high_unc = np.percentile(long_samples, (100 + 68.27) / 2, axis=0) - long_means
    low_unc = long_means - np.percentile(long_samples, (100 - 68.27) / 2, axis=0)

    return long_means, (high_unc + low_unc) / 2


def min_func_ls_sq(params, result):
    if result.is_set():
        raise Exception()
    thetas = [params[f"d{i}"].value for i in range(N_FREE_DIM)]
    ll = -log_probability(thetas)
    return ll

def get_single_density(seed, info, result):
    local_rng = random.Random()
    local_rng.seed(seed)
    while not result.is_set():
        val = [local_rng.random() * 4 * info.mean_density for _ in range(N_FREE_DIM)]
        params = lmfit.Parameters()
        params.add("d0", value=val[0], min=MIN_DENSITY * info.mean_density, max=MAX_DENSITY*info.mean_density)
        params.add("d1", value=val[1], min=MIN_DENSITY * info.mean_density, max=MAX_DENSITY*info.mean_density)
        params.add("d2", value=val[2], min=MIN_DENSITY * info.mean_density, max=MAX_DENSITY*info.mean_density)
        params.add("d3", value=val[3], min=MIN_DENSITY * info.mean_density, max=MAX_DENSITY*info.mean_density)
        params.add("d4", value=val[4], min=MIN_DENSITY * info.mean_density, max=MAX_DENSITY*info.mean_density)
        params.add("d5", value=val[5], min=MIN_DENSITY * info.mean_density, max=MAX_DENSITY*info.mean_density)
        params.add("d6", value=val[6], min=MIN_DENSITY * info.mean_density, max=MAX_DENSITY*info.mean_density)
        params.add("d7", value=val[7], min=MIN_DENSITY * info.mean_density, max=MAX_DENSITY*info.mean_density)
        params.add("d8", value=val[8], min=MIN_DENSITY * info.mean_density, max=MAX_DENSITY*info.mean_density)
        try:
            min_result = lmfit.minimize(min_func_ls_sq, params, args=(result,),method='lbfgsb')
        except:
            print("Thread bailed")
            return
        if min_result.redchi < 1e8:
            print(lmfit.fit_report(min_result))
            print("Thread successfully completed")
            result.set(min_result)
            return
        else:
            print(f"Attempt failed with redchi {min_result.redchi}")

def get_densities_ls_sq(output_name, info):
    threads = []
    result = MinResult()
    print(f"Starting {NUM_THREADS} threads")
    for i in range(NUM_THREADS):
        seed = random.randint(0, 0xffff_ffff_ffff_ffff)
        t = Thread(target=get_single_density, args=(seed, info, result))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    best_result = result.get()
    print("Best result:", best_result)

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
        masks = get_grids_centroid(N_ALL_DIM, asteroid.grid_line, asteroid.indicator_map, asteroid.indicator)[1]
        with open(name + "-grids.npy", 'wb') as f:
            np.save(f, masks)
    else:
        with open(name + "-grids.npy", 'rb') as f:
            masks = np.load(f)
    
    rlms = asteroid.moment_field()

    rlm_mat_complex = np.einsum("iabc,jabc->ji", masks, rlms) * division**3

    # Set order
    rlm_mat = np.zeros((N_ALL_DIM, N_ALL_DIM))
    rlm_mat[0, :] = rlm_mat_complex[8,:].real # K22
    rlm_mat[1, :] = rlm_mat_complex[6,:].real # K20
    rlm_mat[2, :] = rlm_mat_complex[15,:].real / surface_am # R K33
    rlm_mat[3, :] = rlm_mat_complex[15,:].imag / surface_am # I K33
    rlm_mat[4, :] = rlm_mat_complex[14,:].real / surface_am # R K32
    rlm_mat[5, :] = rlm_mat_complex[14,:].imag / surface_am # I K32
    rlm_mat[6, :] = rlm_mat_complex[13,:].real / surface_am # R K31
    rlm_mat[7, :] = rlm_mat_complex[13,:].imag / surface_am # I K31
    rlm_mat[8, :] = rlm_mat_complex[12,:].real / surface_am # K30
    rlm_mat[9, :] = rlm_mat_complex[ 3,:].real * surface_am # R K11
    rlm_mat[10, :] = rlm_mat_complex[3,:].imag * surface_am # I K11
    rlm_mat[11, :] = rlm_mat_complex[2,:].real * surface_am # K10
    rlm_mat[12, :] = rlm_mat_complex[7,:].real # R K21
    rlm_mat[13, :] = rlm_mat_complex[7,:].imag # I K21
    rlm_mat[14, :] = rlm_mat_complex[8,:].imag # I K22
    rlm_mat[15, :] = rlm_mat_complex[0,:].real # K00
    radius_vec = rlm_mat_complex[-1, :].real # radius

    rlm_fixed_inv = inv(rlm_mat[np.arange(9, 16),:][:,np.arange(9, 16)])
    rlm_cross = rlm_mat[np.arange(9, 16),:][:,np.arange(0, 9)]
    rlm_prod = rlm_fixed_inv @ rlm_cross
    return AsteroidInfo(mean_density, data, data_inv_covs, rlm_mat, radius_vec, rlm_fixed_inv, rlm_prod, masks)

def display(densities, true_densities, uncertainty_ratios,
    grid_line, error, asteroid_name, duration=5):

    if not os.path.isdir(f"../figs/{asteroid_name}"):
        os.mkdir(f"../figs/{asteroid_name}")

    warnings.filterwarnings("ignore")

    densities /= np.nanmean(densities)

    if true_densities is not None:
        true_densities /= np.nanmean(true_densities)
        ratios = (densities - true_densities) / (densities * uncertainty_ratios)
        make_slices(ratios, grid_line, "$\\Delta\\sigma$", 'coolwarm', f"../figs/{asteroid_name}/fe-r", error, percentile=95, balance=True)
        make_gif(ratios, grid_line, "$\\Delta\\sigma$", 'coolwarm', f"../figs/{asteroid_name}/fe-r.gif", duration=duration, percentile=95, balance=True)
        difference = (true_densities - densities)

    print("Plotting density")
    make_slices(densities, grid_line, "$\\rho$", 'plasma', f"../figs/{asteroid_name}/fe-d", error)
    make_gif(densities, grid_line, "$\\rho$", 'plasma', f"../figs/{asteroid_name}/fe-d.gif", duration)
    
    print("Plotting uncertainty")
    make_slices(uncertainty_ratios, grid_line, "$\\sigma_\\rho / \\rho$", 'Greys_r', f"../figs/{asteroid_name}/fe-u", error, 95)
    make_gif(uncertainty_ratios, grid_line, "$\\sigma_\\rho / \\rho$", 'Greys_r', f"../figs/{asteroid_name}/fe-u.gif", duration, 95)

    if true_densities is not None:
        print("Plotting differences")
        make_slices(difference, grid_line, "$\\Delta\\rho$", 'PuOr_r', f"../figs/{asteroid_name}/fe-s", error, 95, balance=True)
        make_gif(difference, grid_line, "$\\Delta\\rho$", 'PuOr_r', f"../figs/{asteroid_name}/fe-s.gif", duration, 95, balance=True)

    warnings.filterwarnings("default")

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
    print("Grid line", asteroid.grid_line)
    asteroid_info = load(name, asteroid, surface_am, sample_path, division, generate)

    if asteroid_info is None:
        return np.nan
    
    means, unc = get_densities_mcmc(name+"-fe", asteroid_info, generate)
    if np.any(np.isnan(unc)):
        return np.nan

    if map:
        densities, uncertainty_ratios = get_map(asteroid_info, means, unc, asteroid)
        true_densities = asteroid.get_true_densities()
        true_densities[~asteroid.indicator_map] = np.nan
        error = log_probability(means[:N_FREE_DIM], asteroid_info) / N_FREE_DIM * -2
        display(densities, true_densities, uncertainty_ratios,
        asteroid.grid_line, error, name)

    return unc
    
if __name__ == "__main__":
    k22, k20, surface_am = -0.05200629, -0.2021978, 1000 # For the shape
    for i in range(20):
        pipeline(f"den-core-sph-{i}", "../samples/den-core-sph-0-samples.npy", Indicator.ell(surface_am, k22, k20),
            surface_am, DIVISION, MAX_RADIUS, True, used_bulk_am=922.9234884822591)