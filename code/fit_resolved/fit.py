TEST = False

import numpy as np
import matplotlib.pyplot as plt
import emcee, time, sys
if not TEST:
    import asteroids_0_2, asteroids_0_3
if TEST:
    import test_0_2 as asteroids_0_2
    import test_0_2 as asteroids_0_3
from multiprocessing import Pool
from random_vector import *
from scipy import optimize, linalg
from mpmath import mp
import plotille
import display, collect

import numdifftools as nd


REPORT_INITIAL = True
PLOT_POSES = True
RECALC_HESS = True

EARTH_RADIUS = 6_370_000
N_WALKERS = 32
MAX_N_STEPS = 100_000
NUM_MINIMIZE_POINTS = 48
NUM_FITS = [3, 1]
EPSILON = 1e-10 # If ABNORMAL_TERMINATION_IN_LNSRCH occurs, EPSILON may be too large.
MIN_SPREAD = 1e-4 ** 2

DISTANCE_RATIO_CUT = 10
MIN_THETA_DIST = 0.01

if len(sys.argv) not in [2, 3]:
    raise Exception("Please pass a file to describe the fit")
output_name = sys.argv[1]
f = open("../../staged/" + output_name+".dat", 'r')
ASTEROIDS_MAX_J, ASTEROIDS_MAX_K = f.readline().split(', ')
ASTEROIDS_MAX_J = int(ASTEROIDS_MAX_J)
ASTEROIDS_MAX_K = int(ASTEROIDS_MAX_K)
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

NUM_FITS = NUM_FITS[:ASTEROIDS_MAX_K-1]

if REPORT_INITIAL:
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
    if ASTEROIDS_MAX_K == 3:
        resolved_data = asteroids_0_3.simulate(cadence, jlms, theta[1:], radius,
            spin[0], spin[1], spin[2], theta[0], impact_parameter, speed, -1)
    elif ASTEROIDS_MAX_K == 2:
        resolved_data = asteroids_0_2.simulate(cadence, jlms, theta[1:], radius,
            spin[0], spin[1], spin[2], theta[0], impact_parameter, speed, -1)
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

cadence_cut = len(asteroids_0_2.simulate(cadence, jlms, theta_true[1:], radius,
    spin[0], spin[1], spin[2], theta_true[0], impact_parameter, speed, DISTANCE_RATIO_CUT))

plt.figure(figsize=(12, 4))
x_display = np.arange(len(y) / 3)
plt.errorbar(x_display, y[::3], yerr=yerr[::3], label = 'x', fmt='.')
plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], label = 'y', fmt='.')
plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], label = 'z', fmt='.')
plt.xlabel("Time (Cadences)")
plt.ylabel("Spin (rad/s)")
plt.axvline(x=cadence_cut // 3, color='k')
plt.legend()
#plt.show()

def is_identity(h):
    for i in range(len(h)):
        if abs(h[i][i] - 1) > EPSILON:
            return False
    for i in range(len(h)-1):
        for j in range(i+1, len(h)):
            if abs(h[i][j]) > EPSILON:
                return False
    return True


####################################################################
# Minimize likelihood
####################################################################
print()
y_min = y[:cadence_cut]
yerr_min = yerr[:cadence_cut]

def minimize_function(theta, simulate_func):
    resolved_data = simulate_func(cadence, jlms, theta[1:], radius,
        spin[0], spin[1], spin[2], theta[0], impact_parameter, speed, DISTANCE_RATIO_CUT)
    return np.asarray(resolved_data)

def minimize_log_prob(float_theta, fix_theta, simulate_func):
    # Normal likelihood
    theta = list(fix_theta) + list(float_theta)
    try:
        model = minimize_function(theta, simulate_func)
    except RuntimeError:
        return 1e10 # Zero likelihood

    # Flat priors
    for i, t in enumerate(theta):
        if t > theta_high[i] or t < theta_low[i]:
            return 1e10 # Zero likelihood

    # Return redchi
    return np.sum((y_min - model) ** 2 /  yerr_min ** 2)

def get_minimum(arg):
    point, fix_theta, l, bounds = arg
    if l == 2:
        simulate_func = asteroids_0_2.simulate
    elif l == 3:
        simulate_func = asteroids_0_3.simulate
    bfgs_min = optimize.minimize(minimize_log_prob, point, args=(fix_theta, simulate_func),
        method='L-BFGS-B', options={"eps": EPSILON}, bounds=bounds)
    if not bfgs_min.success:
        print("One of the minimum finding points failed.")
        #print(bfgs_min)
        return None, None, None, None

    if RECALC_HESS:
        def hess_func(min_coords):
            r = minimize_log_prob(min_coords, fix_theta, simulate_func)
            return r

        grad = nd.Gradient(hess_func, step=EPSILON)(bfgs_min.x)
        hess = nd.Hessian(hess_func, step=EPSILON)(bfgs_min.x)

        try:
            hess_inv = linalg.inv(hess)
            #print("Min hess:", bfgs_min.hess_inv.todense())
        except:
            print("Singular matrix: {}".format(hess))
            return None, None, None, None
    else:
        if is_identity(bfgs_min.hess_inv.todense()):
            print("One of the Hessians was the identity.")
            return None, None, None, None
        hess_inv = bfgs_min.hess_inv.todense()

    new_evals = []
    new_evecs = []
    evals, evecs = mp.eigsy(mp.matrix(hess_inv))
    for e in evals:
        ### Correct for non positive definite hessians
        #new_evals.append(float(e))
        new_evals.append(abs(float(e)))
    if np.any(np.asarray(new_evals) < 0):
        print("One of the Hessians was not positive definite")
        return None, None, None, None

    for k in range(int(len(evecs))):
        new_evecs.append(np.array([evecs[j, k] for j in range(int(len(evecs)))],
        dtype=np.float64))

    # Return redchi, minimizing params, hessian
    return bfgs_min.fun / len(y_min), bfgs_min.x, new_evals, new_evecs

def minimize(l, num_return, fix_theta):
    assert l <= ASTEROIDS_MAX_K
    assert len(fix_theta) == max((l)**2 - 6, 0)

    if l == 2:
        simulate = asteroids_0_2.simulate
    elif l == 3:
        simulate = asteroids_0_3.simulate
    else:
        raise Exception("l={} is not supported".format(l))

    theta_indices = range(max(0, (l)**2 - 6), (l+1)**2 - 6)
    bounds = np.asarray([(theta_low[i], theta_high[i]) for i in theta_indices])
    bound_widths = np.asarray([theta_high[i] - theta_low[i] for i in theta_indices])

    # Choose the seed parameters
    parameter_points = []

    while len(parameter_points) < NUM_MINIMIZE_POINTS:
        randoms = np.asarray([random.random() for i in theta_indices])
        point = np.asarray(theta_low)[theta_indices] + bound_widths * randoms
        try:
            model = minimize_function(point, simulate)
        except RuntimeError:
            continue
        parameter_points.append((point, fix_theta, l, bounds))

    # Perform the minimization
    with Pool() as pool:
        results = pool.map(get_minimum, parameter_points)

    # Extract the lowest redchis
    sorted_results = sorted(results, key=
        lambda x: x[0] if (x[0] is not None and not np.isnan(x[0])) else 1e10)

    distinct_results = []
    for redchi, theta, evals, evecs in sorted_results:
        if redchi is None:
            continue
        if len(distinct_results) >= num_return:
            break
        choose = True
        for distinct_theta, _, _, _ in distinct_results:
            if np.sum([(theta[i] - distinct_theta[i])**2 for i in range(len(theta))]) < MIN_THETA_DIST**2:
                choose = False
                break
        if choose:
            distinct_results.append((theta, evals, evecs, redchi))
    return distinct_results

# Stepped minimization
queue = [([], [], [], 0)]
for i, num_fits in enumerate(NUM_FITS):
    new_queue = []
    for fix_theta, evals, evecs, _ in queue:
        tier_results = minimize(i+2, num_fits, fix_theta)
        for result_theta, result_evals, result_evecs, result_redchi in tier_results:
            print("Deg {} redchi: {}".format(i, result_redchi))
            new_queue.append((fix_theta + list(result_theta), evals + list(result_evals),
                evecs + list(result_evecs), result_redchi))
    queue = new_queue

# Stitch together the hessians
kernel = []
for theta, evals, evecs, redchi in queue:
    resized_evecs = []
    max_len = (1 + ASTEROIDS_MAX_K)**2 - 6
    for e in evecs:
        l_value = (len(e) - 1) // 2
        start_index = l_value**2 - 6
        if l_value < 2:
            l_value = 2
            start_index = 0
        evec = [0] * start_index + list(e) + [0] * (max_len - start_index - len(e))
        resized_evecs.append(np.array(evec))
    print("The kernel includes a point with theta {}, redchi {}".format(theta, redchi))
    kernel.append((theta, evals, np.array(resized_evecs)))

print("There are {} MCMC starting points, and there should be {}".format(len(kernel), np.prod(NUM_FITS)))

####################################################################
# Run MCMC
####################################################################
print()

def populate(evals, diagonalizer, count):
    diagonal_points = np.sqrt(np.maximum(MIN_SPREAD, evals)) * (np.random.randn(count * N_DIM).reshape(count, N_DIM))
    global_points = np.asarray([np.matmul(diagonalizer, d) for d in diagonal_points])
    return global_points

def mcmc_fit(theta_start, evals, evecs, index):
    backend = emcee.backends.HDFBackend(output_name+"-{}.h5".format(index))

    if not reload:
        pos = populate(evals, evecs, N_WALKERS) + theta_start
        if PLOT_POSES:
            for i in range(len(pos[0])):
                print(plotille.histogram(
                    pos[:,i],
                    bins=8,
                    width=80,
                    height=15,
                    X_label='theta {}'.format(i),
                    Y_label='Counts',
                    linesep='\n',
                ))
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

for i, (theta, evals, evecs) in enumerate(kernel):
    pass
    #mcmc_fit(theta, evals, evecs, i)


####################################################################
# Process data
####################################################################
print()

i = 0
while True:
    print("Showing i = {}".format(i))
    try:
        disp = display.Display("{0}".format(output_name), "{0}-{1}".format(output_name, i))
    except Exception as e:
        if i == 0:
            raise e
        break
    disp.show_redchi()
    disp.show_params()
    disp.show_corner()
    disp.show_compare()
    disp.show_results()
    plt.show()
    if not collect.collect(output_name + "-" + str(i), output_name):
        break
    del disp
    i += 1
