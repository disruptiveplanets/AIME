import numpy as np
import matplotlib.pyplot as plt
import emcee, asteroids, time, sys
from multiprocessing import Pool

EARTH_RADIUS = 6370000
EARTH_MASS =5.972e24
CADENCE = 60.0
REGENERATE_DATA = True
N_WALKERS = 32
N_STEPS = 5000
REL_SIGMA = 0.1

ASTEROIDS_MAX_K = 2 # Remember to change the counterpart in backend.hpp

MULTIPROCESSING = True

reload = False
if len(sys.argv) == 2:
    if sys.argv[1] == "reload":
        reload = True
        REGENERATE_DATA = False

#np.random.seed(123)

spin = [0.00012, 0.00022, 0.00032]
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
jlms = [5.972e24, 0, 0, 4.972e22]
klms = [
    1e6, 1e5, 5e5,
    0, 0, 0, 0, 0, 0, 0 #m3
    ]

theta_true = (
    0, 1.2e6, 1.1e5, -4.9e5,
)
theta_start = (
    0.1, 1.0e6, 1.0e5, -5.0e5,
)
theta_range = (
    (-np.pi, np.pi), (0.5e6, 2e6), (0.5e5, 2.0e5), (-2.5e5, -1.0e6),
)

def fit_function(theta):
    #start = time.time()
    resolved_data = asteroids.simulate(CADENCE, jlms, theta[1:],
        spin[0], spin[1], spin[2], theta[0], impact_parameter, speed)
    #print(time.time() - start)
    return np.asarray(resolved_data)

def get_list_from_file(filename):
    f = open(filename, 'r')
    l = [float(line) for line in f.readlines()]
    f.close()
    return np.asarray(l)
def save_to_file(filename, l):
    f = open(filename, 'w')
    for entry in l:
        f.write(str(entry) + '\n')
    f.close()

def randomize(y):
    yerr = y * REL_SIGMA
    return y + np.random.randn(len(y)) * yerr, yerr

# Generate some synthetic data from the model.
if REGENERATE_DATA:
    start = time.time()
    y = fit_function(theta_true)
    print("Took {} s".format(time.time() - start))
    save_to_file("simulated-data.dat", y)
else:
    y = get_list_from_file("simulated-data.dat")

x = np.arange(len(y))
y, yerr = randomize(y)

plt.figure(figsize=(12, 4))
x_display = np.arange(len(y) / 3)
plt.errorbar(x_display, y[::3], yerr=yerr[::3], label = 'x', fmt='.')
plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], label = 'y', fmt='.')
plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], label = 'z', fmt='.')
plt.xlabel("Time (Cadences)")
plt.ylabel("Spin (rad/s)")
plt.legend()
plt.show()

# Confirm that theta start works.
asteroids.simulate(CADENCE, jlms, theta_start[1:],
    spin[0], spin[1], spin[2], theta_start[0], impact_parameter, speed)

def log_likelihood(theta, y, yerr):
    # Normal likelihood
    try:
        model = fit_function(theta)
    except RuntimeError:
        return -np.inf # Zero likelihood
    sigma2 = yerr ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

def log_prior(theta):
    for i, param in enumerate(theta):
        if param > max(theta_range[i]) or param < min(theta_range[i]):
            return -np.inf
    return 0.0

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, y, yerr)

pos = theta_start + 1e-4 * np.random.randn(N_WALKERS, len(theta_start))
nwalkers, ndim = pos.shape

save_filename = "asteroids.h5"
backend = emcee.backends.HDFBackend(save_filename)
if not reload:
    backend.reset(nwalkers, ndim)
else:
    print("Initial size: {}".format(backend.iteration))

if MULTIPROCESSING:
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,
            args=(x, y, yerr), backend=backend, pool=pool)

        if not reload:
            sampler.run_mcmc(pos, N_STEPS, progress=True)
        else:
            sampler.run_mcmc(None, N_STEPS, progress=True)
else:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,
        args=(x, y, yerr), backend=backend)
    if not reload:
        sampler.run_mcmc(pos, N_STEPS, progress=True)
    else:
        sampler.run_mcmc(None, N_STEPS, progress=True)

if reload:
    print("New size: {}".format(backend.iteration))
