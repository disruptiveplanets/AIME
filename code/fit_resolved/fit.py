import numpy as np
import matplotlib.pyplot as plt
import emcee, asteroids, time, sys
from multiprocessing import Pool
from random_vector import *
from scipy import optimize, linalg


ASTEROIDS_MAX_K = 3 # Remember to change the counterpart in backend.hpp
ASTEROIDS_MAX_J = 0 # Remember to change the counterpart in backend.hpp
EARTH_RADIUS = 6_370_000
N_WALKERS = 32
MAX_N_STEPS = 50_000
NUM_MINIMIZE_POINTS = 48
NUM_FITS = 3
EPSILON = 1e-10 # If ABNORMAL_TERMINATION_IN_LNSRCH occurs, EPSILON may be too large.

DISTANCE_RATIO_CUT = 2
MIN_THETA_DIST = 0.01

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
        spin[0], spin[1], spin[2], theta[0], impact_parameter, speed, -1)
    return np.asarray(resolved_data)

def minimize_function(theta):
    resolved_data = asteroids.simulate(cadence, jlms, theta[1:], radius,
        spin[0], spin[1], spin[2], theta[0], impact_parameter, speed, DISTANCE_RATIO_CUT)
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

cadence_cut = len(minimize_function(theta_true))

print("DOF:", len(y))

plt.figure(figsize=(12, 4))
x_display = np.arange(len(y) / 3)
plt.errorbar(x_display, y[::3], yerr=yerr[::3], label = 'x', fmt='.')
plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], label = 'y', fmt='.')
plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], label = 'z', fmt='.')
plt.xlabel("Time (Cadences)")
plt.ylabel("Spin (rad/s)")
plt.axvline(x=cadence_cut // 3, color='k')
plt.legend()
plt.show()

####################################################################
# Minimize likelihood
####################################################################
y_min = y[:cadence_cut]
yerr_min = yerr[:cadence_cut]

print()
def redchi(theta):
    # Normal likelihood
    try:
        model = minimize_function(theta)
    except RuntimeError:
        return 1e10 # Zero likelihood
    return np.sum((y_min - model) ** 2 /  yerr_min ** 2) / len(y_min)

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
    bfgs_min = optimize.minimize(redchi, point, method='L-BFGS-B', options={"eps": EPSILON}, bounds=bounds)
    if not bfgs_min.success:
        print("One of the minimum finding points failed.")
        print(bfgs_min)
    try:
        return bfgs_min.fun, bfgs_min.x, linalg.inv(bfgs_min.hess_inv.todense())
    except:
        print("Something broke (variables not defined or matrix inversion failed)")
        return None, None, None


with Pool() as pool:
    results = pool.map(get_minimum, parameter_points)

print("All red chis:")
for like, theta, _ in results:
    if like is None: continue
    print(like, "at parameters", theta)

def populate(sigmas, count):
    evals, diagonalizer = linalg.eigh(sigmas)
    print(1 / evals)
    diagonal_points = 1/np.sqrt(np.abs(evals)) * (np.random.randn(count * N_DIM).reshape(count, N_DIM))
    global_points = np.asarray([np.matmul(diagonalizer, d) for d in diagonal_points])
    return global_points

def populate_ball(sigmas, count):
    evals, diagonalizer = linalg.eigh(sigmas)
    weights = [np.sum([diagonalizer[i][j] / evals[j] for j in range(N_DIM)]) for i in range(N_DIM)]
    print(weights)
    global_points = np.sqrt(np.abs(weights)) * (np.random.randn(count * N_DIM).reshape(count, N_DIM))
    return global_points


####################################################################
# Run MCMC
####################################################################

def mcmc_fit(theta_start, hess, index):
    print()
    backend = emcee.backends.HDFBackend(output_name+"-{}.h5".format(index))

    if not reload:
        pos = populate(-hess, N_WALKERS) + theta_start
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

sorted_results = sorted(results, key=lambda x: x[0])
distinct_results = []

for redchi, theta, hess in sorted_results:
    if len(distinct_results) >= NUM_FITS:
        break
    for _, accepted_theta, accepted_hess in distinct_results:
        if np.sum([(theta[i] - accepted_theta[i])**2 for i in range(len(theta))]) < MIN_THETA_DIST * MIN_THETA_DIST:
            break
    print("Chose redchi", redchi, "at", theta)
    distinct_results.append((redchi, theta, hess))

for i, (_, theta, hess) in enumerate(distinct_results):
    mcmc_fit(theta, hess, i)
