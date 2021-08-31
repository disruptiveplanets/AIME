import numpy as np
import matplotlib.pyplot as plt
import emcee, asteroids, time, sys
from multiprocessing import Pool
from random_vector import *
import scipy


ASTEROIDS_MAX_K = 2 # Remember to change the counterpart in backend.hpp
ASTEROIDS_MAX_J = 0 # Remember to change the counterpart in backend.hpp
EARTH_RADIUS = 6370000
N_WALKERS = 32
MAX_N_STEPS = 50_000
NUM_MINIMIZE_POINTS = 8

if len(sys.argv) not in [2, 3]:
    raise Exception("Please pass a file to describe the fit")
output_name = sys.argv[1]
f = open("../../staged/" + output_name+".dat", 'r')
cadence = int(f.readline())
impact_parameter = EARTH_RADIUS * int(f.readline())
radius = float(f.readline())
speed = float(f.readline())
spin = [float(x) for x in f.readline().split(',')]
jlms = [float(x) for x in f.readline().split(',')]
theta_true = [float(x) for x in f.readline().split(',')]
theta_high = np.asarray([float(x) for x in f.readline().split(',')])
theta_low = np.asarray([float(x) for x in f.readline().split(',')])

sigma = float(f.readline()) * np.sqrt(spin[0]**2 + spin[1]**2 + spin[2]**2)
while output_name[-1] == '\n':
    output_name = output_name[:-1]
f.close()
assert(len(theta_true) == len(theta_high) == len(theta_low))
assert(len(theta_true) == (ASTEROIDS_MAX_K + 1)**2 - 6)
assert(len(jlms) == (ASTEROIDS_MAX_J + 1)**2)
assert(np.all(theta_high > theta_low))

print("Cadence {}, impact parameter {}, speed {}".format(cadence, impact_parameter, speed))
print("Spin", spin)
print("Jlms", jlms)
print("Radius", radius)
print("Theta true", theta_true)
print("Theta high", theta_high)
print("Theta low", theta_low)
print("Sigma", sigma)
print("Name", output_name)
N_DIM = len(theta_true)

reload = False
if len(sys.argv) == 3 and sys.argv[2] == "reload":
    reload = True
    REGENERATE_DATA = False

def fit_function(theta):
    resolved_data = asteroids.simulate(cadence, jlms, theta[1:], radius,
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

def log_probability(theta, y, yerr):
    prior = log_prior(theta)
    like = log_likelihood(theta, y, yerr)
    if not np.isfinite(prior) or not np.isfinite(like):
        return -np.inf
    return prior + like

####################################################################
# Generate synthetic data
####################################################################
start = time.time()
y = fit_function(theta_true)
print("Data generation took {} s".format(time.time() - start))
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

####################################################################
# Minimize likelihood
####################################################################
print()
def minus_log_like(theta, y, yerr):
    # Normal likelihood
    try:
        model = fit_function(theta)
    except RuntimeError:
        return 1e10 # Zero likelihood
    return np.sum((y - model) ** 2 /  yerr ** 2)

bounds = np.asarray([(theta_low[i], theta_high[i]) for i in range(len(theta_high))])
bound_widths = np.asarray([theta_high[i] - theta_low[i] for i in range(len(theta_high))])

parameter_points = []

while len(parameter_points) < NUM_MINIMIZE_POINTS:
    point = np.asarray(theta_low) + bound_widths * np.asarray([random.random() for i in range(len(theta_high))])
    try:
        model = fit_function(point)
    except RuntimeError:
        continue
    parameter_points.append(point)

def get_minimum(point):
    bfgs_min = scipy.optimize.minimize(minus_log_like, point, args=(y, yerr),
        method='L-BFGS-B', bounds=bounds, options={"eps": 1e-10})
    if not bfgs_min.success:
        print("One of the minimum finding points failed.")
        return None, None, None
    else:
        try:
            scipy.linalg.inv(bfgs_min.hess_inv.todense())
        except:
            print("Matrix could not be inverted")
            return None, None, None
        return bfgs_min.fun, bfgs_min.x, scipy.linalg.inv(bfgs_min.hess_inv.todense())


with Pool() as pool:
    results = pool.map(get_minimum, parameter_points)

min_like = None
min_theta = None
min_hessian = None
for like, theta, hess in results:
    if like is None: continue
    if min_like is None or like < min_like:
        min_like = like
        min_theta = theta
        min_hessian = hess

print("Maximum log likelihood was {} (redchi: {}) at parameters {}".format(
    -min_like, min_like / len(y), min_theta))

theta_start = min_theta

def populate(sigmas, count):
    evals, diagonalizer = scipy.linalg.eigh(sigmas)
    print(1 / evals)
    diagonal_points = 1/np.sqrt(np.abs(evals)) * (np.random.randn(count * N_DIM).reshape(count, N_DIM))
    global_points = np.asarray([np.matmul(diagonalizer, d) for d in diagonal_points])
    return global_points

def populate_ball(sigmas, count):
    evals, diagonalizer = scipy.linalg.eigh(sigmas)
    weights = [np.sum([diagonalizer[i][j] / evals[j] for j in range(N_DIM)]) for i in range(N_DIM)]
    print(weights)
    global_points = np.sqrt(np.abs(weights)) * (np.random.randn(count * N_DIM).reshape(count, N_DIM))
    return global_points


####################################################################
# Run MCMC
####################################################################
print()
backend = emcee.backends.HDFBackend(output_name+".h5")

if not reload:
    pos = populate(-min_hessian, N_WALKERS)
    backend.reset(N_WALKERS, N_DIM)
else:
    pos=None
    print("Initial size: {}".format(backend.iteration))

old_tau = np.inf
with Pool() as pool:
    sampler = emcee.EnsembleSampler(N_WALKERS, N_DIM, log_probability,
        args=(y, yerr), backend=backend, pool=pool)

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
    sampler._previous_state = sample

if reload:
    print("New size: {}".format(backend.iteration))
else:
    print("Done")
