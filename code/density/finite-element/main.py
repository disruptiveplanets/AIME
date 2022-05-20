import emcee, corner, sys, lmfit, random, os, time
import numpy as np
from scipy.linalg import pinvh, inv
from matplotlib import pyplot as plt
from multiprocessing import Pool, Lock
from grids import get_grids_centroid
sys.path.append("..")
from core import Asteroid, Indicator
from scipy.optimize import minimize
from threading import Thread
from numdifftools import Hessian

N_WALKERS = 32
MAX_L = 3
N_ALL_DIM = (MAX_L + 1)**2
N_FREE_DIM = N_ALL_DIM - 7
MAX_N_STEPS = 100_000
DIVISION = 49
RLM_EPSILON = 1e-20
MAX_RADIUS = 2000
SMALL_SIGMA = 1e-5
CERTAIN_INDICES = [0, 1, 2, 3, 5, 7, -1]
MAX_DENSITY = 10
MIN_DENSITY = 0
UNCERTAINTY_RATIO = 0.005
AM = 1000
VERY_LARGE_SLOPE = 1e30
NUM_THREADS = os.cpu_count()

def get_cov(path):
    with open(path, 'rb') as f:
        flat_samples = np.load(f).reshape(-1, N_FREE_DIM + 1)[:, 1:]
    cov = np.cov(flat_samples.transpose())
    data = np.mean(flat_samples, axis=0)
    return data, cov

def get_theta_long(theta_short):
    # Get the densities consistent with making the mass 1 and com 0 and rotation
    return np.append(theta_short, rlm_fixed_inv[:,-1] - rlm_prod @ theta_short)

def get_klms(theta_long):
    unscaled_klms = rlm_mat @ theta_long
    radius_sqr = radius_vec @ theta_long
    scaled_klms = np.array(unscaled_klms) / radius_sqr
    scaled_klms[-1] *= radius_sqr # Do not scale mass term
    return scaled_klms


def log_prior(theta_long):
    if np.any(theta_long < mean_density * MIN_DENSITY):
        return VERY_LARGE_SLOPE * np.min(theta_long) / mean_density
    if np.any(theta_long > mean_density * MAX_DENSITY):
        return -VERY_LARGE_SLOPE * (np.max(theta_long) / mean_density - MAX_DENSITY)
    return 0.0
    
def log_like(theta_long):
    diff_klms = get_klms(theta_long)[:N_FREE_DIM] - data
    return -0.5 * diff_klms.transpose() @ data_inv_covs @ diff_klms

def log_probability(theta_short):
    theta_long = get_theta_long(theta_short)
    lp = log_prior(theta_long)
    ll = log_like(theta_long)
    return ll + lp

def mcmc_fit(theta_start, output_name, generate=True):
    if generate:
        backend = emcee.backends.HDFBackend(output_name+".h5")
        backend.reset(N_WALKERS, N_FREE_DIM)
        old_tau = np.inf

        with Pool() as pool:
            sampler = emcee.EnsembleSampler(N_WALKERS, N_FREE_DIM, log_probability, backend=backend, pool=pool)

            pos = np.zeros((N_WALKERS, N_FREE_DIM))
            for i in range(N_FREE_DIM):
                pos[:,i] = np.random.randn(N_WALKERS) * UNCERTAINTY_RATIO * mean_density + theta_start[i]
            for p in pos:
                print("Log prob of pos:", log_probability(p))


            for sample in sampler.sample(pos, iterations=MAX_N_STEPS, progress=True):
                if sampler.iteration % 200 == 0:
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
                if sampler.iteration % 1000 == 0:
                    print("Last log prob", np.mean(sampler.get_last_sample().log_prob))

            sampler._previous_state = sample
    else:
        backend = emcee.backends.HDFBackend(output_name+".h5", read_only=True)
        sampler = emcee.EnsembleSampler(N_WALKERS, N_FREE_DIM, log_probability, backend=backend)

    return sampler

class MinResult:
    def __init__(self):
        self.lock = Lock()
        self.val = None
    
    def set(self, v):
        self.lock.acquire()
        self.val = v
        self.lock.release()

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

def minimize_func(theta, result):
    if result.is_set():
        raise Exception()
    return -log_probability(theta)

def min_func_mcmc(seed, result):
    local_rng = random.Random()
    local_rng.seed(seed)
    while not result.is_set():
        val = [local_rng.random() * 4 * mean_density for _ in range(N_FREE_DIM)]
        try:
            min_result = minimize(minimize_func, x0=val, method="Nelder-Mead", args=(result,), options = {"maxiter": 500 * len(val)})
        except:
            print("Thread bailed")
            return
        if min_result.success and min_result.fun < 1e4:
            print(f"Thread successfully completed with redchi {min_result.fun}")
            result.set(min_result.x)
            return
        else:
            print(f"Attempt failed with redchi {min_result.fun}")

def get_theta_start_mcmc():
    threads = []
    result = MinResult()
    print(f"Starting {NUM_THREADS} threads")
    for i in range(NUM_THREADS):
        seed = random.randint(0, 0xffff_ffff_ffff_ffff)
        t = Thread(target=min_func_mcmc, args=(seed, result))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    result = result.get()
    print(Hessian(lambda theta: -log_probability(theta))(result))
    return result
    
def get_densities_mcmc(output_name, generate=True):
    if generate:
        theta_start = get_theta_start_mcmc()
        print("Theta start:", theta_start)
        sampler = mcmc_fit(theta_start, output_name, True)
        try:
            tau = sampler.get_autocorr_time()
            print("Autocorrelation time", tau)
            flat_samples = sampler.get_chain(discard=2 * np.max(tau), thin=np.max(tau) / 2, flat=True)
        except Exception:
            flat_samples = sampler.get_chain(discard=1000, thin=1, flat=True)

        with open(output_name + ".npy", 'wb') as f:
            np.save(f, flat_samples)
    else:
        with open(output_name + ".npy", 'rb') as f:
            flat_samples = np.load(f)

    fig = corner.corner(flat_samples)
    fig.savefig(output_name + ".png")

    means = np.mean(flat_samples, axis=0)
    unc = np.std(flat_samples, axis=0)
    return means, unc


def min_func_ls_sq(params, result):
    if result.is_set():
        raise Exception()
    thetas = [params[f"d{i}"].value for i in range(N_FREE_DIM)]
    ll = -log_probability(thetas)
    return ll

def get_single_density(seed, result):
    local_rng = random.Random()
    local_rng.seed(seed)
    while not result.is_set():
        val = [local_rng.random() * 4 * mean_density for _ in range(N_FREE_DIM)]
        params = lmfit.Parameters()
        params.add("d0", value=val[0], min=MIN_DENSITY * mean_density, max=MAX_DENSITY*mean_density)
        params.add("d1", value=val[1], min=MIN_DENSITY * mean_density, max=MAX_DENSITY*mean_density)
        params.add("d2", value=val[2], min=MIN_DENSITY * mean_density, max=MAX_DENSITY*mean_density)
        params.add("d3", value=val[3], min=MIN_DENSITY * mean_density, max=MAX_DENSITY*mean_density)
        params.add("d4", value=val[4], min=MIN_DENSITY * mean_density, max=MAX_DENSITY*mean_density)
        params.add("d5", value=val[5], min=MIN_DENSITY * mean_density, max=MAX_DENSITY*mean_density)
        params.add("d6", value=val[6], min=MIN_DENSITY * mean_density, max=MAX_DENSITY*mean_density)
        params.add("d7", value=val[7], min=MIN_DENSITY * mean_density, max=MAX_DENSITY*mean_density)
        params.add("d8", value=val[8], min=MIN_DENSITY * mean_density, max=MAX_DENSITY*mean_density)
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

def get_densities_ls_sq(output_name):
    threads = []
    result = MinResult()
    print(f"Starting {NUM_THREADS} threads")
    for i in range(NUM_THREADS):
        seed = random.randint(0, 0xffff_ffff_ffff_ffff)
        t = Thread(target=get_single_density, args=(seed, result))
        t.start()
        threads.append(t)
    for t in threads:
        t.join()
    best_result = result.get()
    print("Best result:", best_result)


if __name__ == "__main__":
    k22a, k20a = -0.05200629, -0.2021978
    path = "../samples/den-asym-0-samples.npy"

    asteroid = Asteroid("fe", path, AM, DIVISION, MAX_RADIUS, Indicator.ell(AM, k22a, k20a), None)

    mean_density = 1 / (np.sum(asteroid.indicator_map) * DIVISION**3)

    data, cov = get_cov(path)
    data_inv_covs = pinvh(cov)
    masks = get_grids_centroid(N_ALL_DIM, asteroid.grid_line, asteroid.indicator_map, asteroid.indicator)[1]
    rlms = asteroid.moment_field()

    rlm_mat_complex = np.einsum("iabc,jabc->ji", masks, rlms) * DIVISION**3

    # Set order
    rlm_mat = np.zeros((N_ALL_DIM, N_ALL_DIM))
    rlm_mat[0, :] = rlm_mat_complex[8,:].real # K22
    rlm_mat[1, :] = rlm_mat_complex[6,:].real # K20
    rlm_mat[2, :] = rlm_mat_complex[15,:].real / AM # R K33
    rlm_mat[3, :] = rlm_mat_complex[15,:].imag / AM # I K33
    rlm_mat[4, :] = rlm_mat_complex[14,:].real / AM # R K32
    rlm_mat[5, :] = rlm_mat_complex[14,:].imag / AM # I K32
    rlm_mat[6, :] = rlm_mat_complex[13,:].real / AM # R K31
    rlm_mat[7, :] = rlm_mat_complex[13,:].imag / AM # I K31
    rlm_mat[8, :] = rlm_mat_complex[12,:].real / AM # K30
    rlm_mat[9, :] = rlm_mat_complex[3,:].real * AM # R K11
    rlm_mat[10, :] = rlm_mat_complex[3,:].imag * AM # I K11
    rlm_mat[11, :] = rlm_mat_complex[2,:].real * AM # K10
    rlm_mat[12, :] = rlm_mat_complex[7,:].real # R K21
    rlm_mat[13, :] = rlm_mat_complex[7,:].imag # I K21
    rlm_mat[14, :] = rlm_mat_complex[8,:].imag # I K22
    rlm_mat[15, :] = rlm_mat_complex[0,:].real # K00
    radius_vec = rlm_mat_complex[-1, :].real # radius

    rlm_fixed_inv = inv(rlm_mat[np.arange(9, 16),:][:,np.arange(9, 16)])
    rlm_cross = rlm_mat[np.arange(9, 16),:][:,np.arange(0, 9)]
    rlm_prod = rlm_fixed_inv @ rlm_cross

    print("Uniform probability:", log_probability([mean_density] * 9))

    means, unc = get_densities_mcmc("fe", False)
    #means, unc = get_densities_ls_sq("fe")

    print("Means", means / np.mean(means))
    print("Standard deviations", unc / np.mean(means))
    print("Ratios", unc / means)
    print("Klm error", get_klms(means)[:N_FREE_DIM] - data)

    ## Still need to propagate uncertainties to the other densities