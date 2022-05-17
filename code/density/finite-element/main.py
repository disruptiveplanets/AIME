import emcee, corner, sys
import numpy as np
from scipy.linalg import pinvh
from matplotlib import pyplot as plt
from multiprocessing import Pool
from grids import get_grids_centroid
sys.path.append("..")
from core import Asteroid, Indicator
from scipy.optimize import minimize

N_WALKERS = 64
MAX_L = 3
N_DIM = (MAX_L + 1)**2 + 1
MAX_N_STEPS = 100_000
DIVISION = 99
RLM_EPSILON = 1e-20
MAX_RADIUS = 2000
SMALL_SIGMA = 1e-5
CERTAIN_INDICES = [0, 1, 2, 3, 5, 7, -1]

def widen(mat):
    for i in CERTAIN_INDICES:
        mat[i][i] = SMALL_SIGMA**2
    return mat

def get_klms(theta):
    i = 0
    weighted_masks = np.einsum("ijkl,i->jkl", masks, theta)
    klms = np.einsum("ijkl,jkl->i", rlms, weighted_masks) * DIVISION**3

    #klms /= klms[0] # Assume mass
    radius = np.sqrt(klms[-1])
    radius = np.sqrt(data[-1])
    i = 0
    for l in range(0, MAX_L + 1):
        for m in range(-MAX_L, MAX_L + 1):
            klms[i] /= radius**l
        i += 1
    return klms

def log_prior(theta):
    if np.any(theta < 0):
        return -np.inf
    if np.any(theta > mean_density * 10):
        return -np.inf
    return 0.0
    
def log_like(theta):
    diff_klms = get_klms(theta) - data
    return -0.5 * np.abs(diff_klms.transpose() @ data_inv_covs @ diff_klms)

def log_probability(theta):
    lp = log_prior(theta)
    ll = log_like(theta)
    if not np.isfinite(lp):
        return -np.inf
    return ll + lp

def mcmc_fit(theta_start, output_name, generate=True):
    if generate:
        backend = emcee.backends.HDFBackend(output_name+".h5")
        backend.reset(N_WALKERS, N_DIM)
        old_tau = np.inf

        with Pool() as pool:
            sampler = emcee.EnsembleSampler(N_WALKERS, N_DIM, log_probability, backend=backend, pool=pool)

            pos = np.zeros((N_WALKERS, N_DIM))
            for i in range(N_DIM):
                pos[:,i] = np.random.randn(N_WALKERS) * 0.5 * mean_density + theta_start[i]

            for sample in sampler.sample(pos, iterations=MAX_N_STEPS, progress=True):
                if sampler.iteration % 100 == 0:
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
            sampler._previous_state = sample
    else:
        backend = emcee.backends.HDFBackend(output_name+".h5", read_only=True)
        sampler = emcee.EnsembleSampler(N_WALKERS, N_DIM, log_probability, backend=backend)

    return sampler

def min_func(theta):
    return -log_probability(theta)
def pool_func(seed):
    min_result = minimize(min_func, x0=seed, method="Nelder-Mead", options = {"maxiter": 500 * len(seed)})
    if min_result.success:
        return min_result.x, min_result.fun
    else:
        return None, np.nan

def get_theta_start():
    seeds = np.random.uniform(low=0, high=4, size=(48, N_DIM)) * mean_density
    with Pool() as pool:
        mins = pool.map(pool_func, seeds)
    min_f = None
    min_x = None
    num_nan = 0
    for x, f in mins:
        if np.isnan(f):
            num_nan += 1
            continue
        if min_f is None or f < min_f:
            min_f = f
            min_x = x
    print(f"{num_nan} / {seeds.shape[0]} minimizations did not succeed.")
    print(min_x)
    return min_x
    

def get_densities(output_name, generate=True):
    if generate:
        theta_start = get_theta_start()
        sampler = mcmc_fit(theta_start, output_name, True)
        try:
            tau = sampler.get_autocorr_time()
            print("Autocorrelation time", tau)
            flat_samples = sampler.get_chain(discard=2 * np.max(tau), thin=np.max(tau) / 2, flat=True)
        except Exception:
            flat_samples = sampler.get_chain(discard=1000, thin=1, flat=True)

        fig = corner.corner(flat_samples)
        fig.savefig(output_name + ".png")

        with open(output_name + ".npy", 'wb') as f:
            np.save(f, flat_samples)
    else:
        with open(output_name + ".npy", 'rb') as f:
            flat_samples = np.load(f)

    means = np.mean(flat_samples, axis=0)
    print(get_klms(means) - data)

if __name__ == "__main__":
    am = 1000
    k22a, k20a = -0.05200629, -0.2021978

    asteroid = Asteroid("fe", "../samples/den-asym-0-samples.npy", am, DIVISION, MAX_RADIUS, Indicator.ell(am, k22a, k20a), None)

    mean_density = 1 / (np.sum(asteroid.indicator_map) * DIVISION**3)

    data = asteroid.data
    data_inv_covs = pinvh(widen(asteroid.sigma_data))
    masks = get_grids_centroid(N_DIM, asteroid.grid_line, asteroid.indicator_map, asteroid.indicator)[1]
    rlms = asteroid.moment_field()

    get_densities("fe", True)
    #plt.show()