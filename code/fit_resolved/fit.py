import numpy as np
import matplotlib.pyplot as plt
import emcee, asteroids, time, sys
from multiprocessing import Pool

EARTH_RADIUS = 6370000
EARTH_MASS =5.972e24
CADENCE = 3600.0
REGENERATE_DATA = False
N_WALKERS = 32
N_STEPS = 5000

MULTIPROCESSING = True

reload = False
if len(sys.argv) == 2:
    if sys.argv[1] == "reload":
        reload = True
        REGENERATE_DATA = False

np.random.seed(123)

spin = 0.000050189
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
jlms = [5.972e24, 5.972e22, -5.972e22, 4.972e22]

theta_true = (
    1.0, 2.0, 3.0, 4.0, 5.0,#m2
    0.142, 1.512, -0.513, 0.0, 3.512, 2.261, -1.523, #m3
)
theta_start = (
    1.0, 1.0, 1.0, 1.0, 1.0,#m2
    0.0, 0.0, 0.0, 0, 0.0, 0.0, 0.0, #m3
)
L = int(np.sqrt(len(theta_true)+4))
theta_range = ((-5, 5),) * len(theta_start)
theta_labels = []
for l in range(2, L):
    for m in range(-l, l+1):
        theta_labels.append("M" + str(l) + ',' + str(m))

def fit_function(theta):
    #start = time.time()
    resolved_data = asteroids.simulate(CADENCE, jlms, theta, spin,
        impact_parameter, speed)
    #print(time.time() - start)
    return resolved_data

def get_list_from_file(filename):
    f = open(filename, 'r')
    l = [float(line) for line in f.readlines()]
    f.close()
    return l
def save_to_file(filename, l):
    f = open(filename, 'w')
    for entry in l:
        f.write(str(entry) + '\n')
    f.close()

# Generate some synthetic data from the model.
if REGENERATE_DATA:
    start = time.time()
    y = fit_function(theta_true)
    print("Took {} s".format(time.time() - start))
    save_to_file("simulated-data.dat", y)
else:
    y = get_list_from_file("simulated-data.dat")
x = np.arange(len(y))
yerr = np.abs(0.1 * np.random.rand(len(y)) * y)
y += yerr * np.random.randn(len(y))

plt.figure(figsize=(12, 4))
x_display = np.arange(len(y) / 3)
plt.errorbar(x_display, y[::3], yerr=yerr[::3], label = 'x', fmt='.')
plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], label = 'y', fmt='.')
plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], label = 'z', fmt='.')
plt.xlabel("Time (Cadences)")
plt.ylabel("Spin (rad/s)")
plt.legend()
plt.show()



def log_likelihood(theta, y, yerr):
    # Normal likelihood
    model = fit_function(theta)
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
