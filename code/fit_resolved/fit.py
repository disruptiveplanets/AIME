import numpy as np
import matplotlib.pyplot as plt
import emcee, asteroids, time, sys
from multiprocessing import Pool
from random_vector import *


ASTEROIDS_MAX_K = 2 # Remember to change the counterpart in backend.hpp
EARTH_RADIUS = 6370000
N_WALKERS = 32
MAX_N_STEPS = 10000

if len(sys.argv) not in [2, 3]:
    raise Exception("Please pass a file to describe the fit")
output_name = sys.argv[1]
f = open("../../staged/" + output_name+".dat", 'r')
cadence = int(f.readline())
impact_parameter = EARTH_RADIUS * int(f.readline())
speed = float(f.readline())
spin = [float(x) for x in f.readline().split(',')]
jlms = [float(x) for x in f.readline().split(',')]
theta_true = [float(x) for x in f.readline().split(',')]
theta_start = [float(x) for x in f.readline().split(',')]
theta_spread = [float(x) for x in f.readline().split(',')]
theta_high = np.asarray([float(x) for x in f.readline().split(',')])
theta_low = np.asarray([float(x) for x in f.readline().split(',')])
sigma = float(f.readline()) * np.sqrt(spin[0]**2 + spin[1]**2 + spin[2]**2)
while output_name[-1] == '\n':
    output_name = output_name[:-1]
f.close()
assert(len(theta_true) == len(theta_start) == len(theta_spread) == len(theta_high) == len(theta_low))
assert(len(theta_true) == (ASTEROIDS_MAX_K + 1)**2 - 5)
assert(np.all(theta_high > theta_low))

print("Cadence {}, impact parameter {}, speed {}".format(cadence, impact_parameter, speed))
print("Spin", spin)
print("Jlms", jlms)
print("Theta true", theta_true)
print("theta start", theta_start)
print("theta spread", theta_spread)
print("Theta high", theta_high)
print("Theta low", theta_low)
print("Sigma", sigma)
print("Name", output_name)
N_DIM = len(theta_start)

reload = False
if len(sys.argv) == 3 and sys.argv[2] == "reload":
    reload = True
    REGENERATE_DATA = False

def fit_function(theta):
    resolved_data = asteroids.simulate(cadence, jlms, theta[1:],
        spin[0], spin[1], spin[2], theta[0], impact_parameter, speed)
    return np.asarray(resolved_data)

def log_likelihood(theta, y, yerr):
    # Normal likelihood
    try:
        model = fit_function(theta)
    except RuntimeError:
        return -np.inf # Zero likelihood
    return -np.sum((y - model) ** 2 /  yerr ** 2)

def log_prior(theta):
    for i, param in enumerate(theta):
        if param > theta_high[i] or param < theta_low[i]:
            return -np.inf
    return 0.0

def log_probability(theta, x, y, yerr):
    prior = log_prior(theta)
    like = log_likelihood(theta, y, yerr)
    if not np.isfinite(prior) or not np.isfinite(like):
        return -np.inf
    return prior + like

# Generate some synthetic data from the model.
start = time.time()
y = fit_function(theta_true)
print("Data generation took {} s".format(time.time() - start))
x = np.arange(len(y))
y, yerr = randomize_rotate(y, sigma)

print("DOF:", len(y))

plt.figure(figsize=(12, 4))
x_display = np.arange(len(y) / 3)
plt.errorbar(x_display, y[::3], yerr=yerr[::3], label = 'x', fmt='.')
plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], label = 'y', fmt='.')
plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], label = 'z', fmt='.')
plt.xlabel("Time (Cadences)")
plt.ylabel("Spin (rad/s)")
plt.legend()
plt.show()

backend = emcee.backends.HDFBackend(output_name+".h5")

if not reload:
    pos = theta_start + theta_spread * np.random.randn(N_WALKERS, N_DIM)
    backend.reset(N_WALKERS, N_DIM)
else:
    pos=None
    print("Initial size: {}".format(backend.iteration))

old_tau = np.inf
with Pool() as pool:
    sampler = emcee.EnsembleSampler(N_WALKERS, N_DIM, log_probability,
        args=(x, y, yerr), backend=backend, pool=pool)

    if reload:
        pos = sampler._previous_state

    for sample in sampler.sample(pos, iterations=MAX_N_STEPS, progress=True):
        if sampler.iteration % 100 != 0:
            continue

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

if reload:
    print("New size: {}".format(backend.iteration))
else:
    print("Done")
