TEST = False

import numpy as np
import matplotlib.pyplot as plt
import emcee, time, sys, logging
if not TEST:
    import asteroids_0_2, asteroids_0_3, asteroids_2_2, asteroids_2_3, asteroids_3_2, asteroids_3_3
if TEST:
    import test_loglike as asteroids_0_2
    import test_loglike as asteroids_0_3
    import test_loglike as asteroids_2_2
    import test_loglike as asteroids_2_3
    import test_loglike as asteroids_3_2
    import test_loglike as asteroids_3_3
from multiprocessing import Pool
import random_vector, random
from scipy import optimize
from scipy.linalg import pinvh
from mpmath import mp
import plotille

import numdifftools as nd

try:
    import coloredlogs
    coloredlogs.install()
except:
    pass

logging.basicConfig(level=logging.INFO)

PLOT_POSES = True
NUM_MINIMIZE_POINTS_PER = 4 # Used to be 8
NUM_L3_MINIMIZATIONS = 1# Changed to 1 to decrease time.
DISTANCE_RATIO_CUTS = [
    [2.0, None],
    [-2.0, None]
]
THRESHOLD_LIKE_RAT = 2




#GM = 1.26686534e17 # Jupiter
GM = 3.986004418e14 # Earth





EARTH_RADIUS = 6_370_000
UNCERTAINTY_MODEL = random_vector.TILT_UNIFORM_TRUE
INTEGRAL_LIMIT_FRAC = 1.0e-3

N_WALKERS = 32
MAX_N_STEPS = 100_000
MIN_THETA_DIST = 1e-4
LARGE_NUMBER = 1e100
MIN_SPACING = np.array([1.0e-6, 1.0e-6, 1.0e-6,
    1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1, 1.0e-1])# Comes in blocks corresponding to the block diagonals
MAX_SPACING = np.array([0.0003, 0.0003, 0.0003,
    LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER, LARGE_NUMBER])
DEFAULT_THIN = 10

if len(sys.argv) not in [2, 3]:
    raise Exception("Please pass a file to describe the fit")
output_name = sys.argv[1]
f = open("../../staged/" + output_name+".txt", 'r')
ASTEROIDS_MAX_J, ASTEROIDS_MAX_K = f.readline().split(', ')
ASTEROIDS_MAX_J = int(ASTEROIDS_MAX_J)
ASTEROIDS_MAX_K = int(ASTEROIDS_MAX_K)
cadence = int(float(f.readline()))
perigee = EARTH_RADIUS * float(f.readline())
radius = float(f.readline())
speed = float(f.readline())
spin = [float(x) for x in f.readline().split(',')]
jlms = [float(x) for x in f.readline().split(',')]
theta_true = [float(x) for x in f.readline().split(',')]
theta_high = np.asarray([float(x) for x in f.readline().split(',')])
theta_low = np.asarray([float(x) for x in f.readline().split(',')])
sigma = [float(d) for d in f.readline().split(", ")]# theta, ratio
last_line = f.readline()
VELOCITY_MUL = 1 if last_line == '' else float(last_line)
while output_name[-1] == '\n':
    output_name = output_name[:-1]
f.close()

SIGMA_FACTOR = (sigma[0] / 0.01)**2#theta

if ASTEROIDS_MAX_K == 2:
    cut_index = 3
elif ASTEROIDS_MAX_K == 3:
    cut_index = 10
else:
    raise Exception("Have not implemented max k >= 4")
MIN_SPACING = MIN_SPACING[:cut_index]
MAX_SPACING = MAX_SPACING[:cut_index]
N_DIM = len(theta_true)
reload = False
if len(sys.argv) == 3 and sys.argv[2] == "reload":
    reload = True
    REGENERATE_DATA = False

def fit_function(theta, target_length=None):
    if ASTEROIDS_MAX_K == 3:
        if ASTEROIDS_MAX_J == 0:
            resolved_data = asteroids_0_3.simulate(cadence, jlms, theta, radius,
                spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS, 0, False, VELOCITY_MUL)
        elif ASTEROIDS_MAX_J == 2:
            resolved_data = asteroids_2_3.simulate(cadence, jlms, theta, radius,
                spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS, 0, False, VELOCITY_MUL)
        elif ASTEROIDS_MAX_J == 3:
            resolved_data = asteroids_3_3.simulate(cadence, jlms, theta, radius,
                spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS, 0, False, VELOCITY_MUL)
    elif ASTEROIDS_MAX_K == 2:
        if ASTEROIDS_MAX_J == 0:
            resolved_data = asteroids_0_2.simulate(cadence, jlms, theta, radius,
                spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS, 0, False, VELOCITY_MUL)
        elif ASTEROIDS_MAX_J == 2:
            resolved_data = asteroids_2_2.simulate(cadence, jlms, theta, radius,
                spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS, 0, False, VELOCITY_MUL)
        elif ASTEROIDS_MAX_J == 3:
            resolved_data = asteroids_3_2.simulate(cadence, jlms, theta, radius,
                spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS, 0, False, VELOCITY_MUL)
    if target_length is not None:
        while len(resolved_data)//3 < target_length:
            resolved_data.append(resolved_data[-3])
            resolved_data.append(resolved_data[-3])
            resolved_data.append(resolved_data[-3])
    return np.asarray(resolved_data).reshape(-1, 3)

def log_likelihood(theta, y, sigma):
    # Normal likelihood
    try:
        model = fit_function(theta, len(y))
    except RuntimeError:
        return -np.inf # Zero likelihood
   
    return random_vector.log_likelihood(UNCERTAINTY_MODEL, y, model, len(y), sigma)

def log_prior(theta):
    for i, param in enumerate(theta):
        if param > theta_high[i] or param < theta_low[i]:
            return -np.inf
    return 0.0

def log_probability(theta, y, sigma):
    prior = log_prior(theta)
    like = log_likelihood(theta, y, sigma)
    if not np.isfinite(prior) or not np.isfinite(like):
        return -np.inf
    return prior + like

def get_data_cut(drc):
    enforce_drc = False if drc is None else True
    drc = 0 if drc is None else drc
    if ASTEROIDS_MAX_J == 0:
        return len(asteroids_0_2.simulate(cadence, jlms, theta_true, radius,
            spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS,
            drc, enforce_drc, VELOCITY_MUL))//3
    elif ASTEROIDS_MAX_J == 2:
        return len(asteroids_2_2.simulate(cadence, jlms, theta_true, radius,
            spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS,
            drc, enforce_drc, VELOCITY_MUL))//3
    elif ASTEROIDS_MAX_J == 3:
        return len(asteroids_3_2.simulate(cadence, jlms, theta_true, radius,
            spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS,
            drc, enforce_drc, VELOCITY_MUL))//3
    raise Exception("Cannot handle this J.")

data_cuts = [[get_data_cut(drc) for drc in tier] for tier in DISTANCE_RATIO_CUTS]

def minimize_function(theta, simulate_func, l_index, cut_index):
    drc = DISTANCE_RATIO_CUTS[l_index][cut_index]
    enforce_drc = False if drc is None else True

    drc = 0 if drc is None else drc
    resolved_data = simulate_func(cadence, jlms, theta, radius,
        spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS, drc, enforce_drc, VELOCITY_MUL)
    
    want_length = data_cuts[l_index][cut_index]
    while len(resolved_data)//3 < want_length:
        resolved_data.append(resolved_data[-3])
        resolved_data.append(resolved_data[-3])
        resolved_data.append(resolved_data[-3])

    while len(resolved_data)//3 > want_length:
        del resolved_data[-1]
        del resolved_data[-1]
        del resolved_data[-1]

    return np.asarray(resolved_data).reshape(-1, 3)

def minimize_log_prob(y, float_theta, fix_theta, simulate_func, l_index, cut_index):
    # Normal likelihood
    theta = list(fix_theta) + list(float_theta)
    try:
        model = minimize_function(theta, simulate_func, l_index, cut_index)
    except RuntimeError as e:
        return LARGE_NUMBER # Zero likelihood

    # Flat priors
    for i, t in enumerate(theta):
        if t > theta_high[i] or t < theta_low[i]:
            return LARGE_NUMBER # Zero likelihood

    # Return chisq
    chisq = 0
    want_length = data_cuts[l_index][cut_index]
    return -random_vector.log_likelihood(UNCERTAINTY_MODEL, y, model, want_length, sigma)

def get_minimum(arg):
    y, point, fix_theta, l, bounds = arg

    redchi_record = []

    if l == 2:
        if ASTEROIDS_MAX_J == 0:
            simulate_func = asteroids_0_2.simulate
        elif ASTEROIDS_MAX_J == 2:
            simulate_func = asteroids_2_2.simulate
        elif ASTEROIDS_MAX_J == 3:
            simulate_func = asteroids_3_2.simulate
    elif l == 3:
        if ASTEROIDS_MAX_J == 0:
            simulate_func = asteroids_0_3.simulate
        elif ASTEROIDS_MAX_J == 2:
            simulate_func = asteroids_2_3.simulate
        elif ASTEROIDS_MAX_J == 3:
            simulate_func = asteroids_3_3.simulate

    # Do an initial fit with cut data
    for cut_index in range(len(DISTANCE_RATIO_CUTS)):
        def minimize_log_prob_cut(x):
            return minimize_log_prob(y, x, fix_theta, simulate_func, l - 2, cut_index)

        start_redchi = 2 *  minimize_log_prob_cut(point) / data_cuts[l-2][cut_index] / 3

        '''bfgs_min = optimize.minimize(minimize_log_prob_cut, point,
            method='L-BFGS-B', bounds=bounds, 
            options={'disp': None,
                'maxls': 20,
                'iprint': -1,
                'gtol': 1e-05,
                'eps': 1e-8,
                'maxiter': 15000,
                'ftol': 2.220446049250313e-09,
                'maxcor': 10,
                'maxfun': 15000})'''

        bfgs_min = optimize.minimize(minimize_log_prob_cut, point, method='Nelder-Mead', options={
            'maxiter': 200 * len(point),
        })

        point = bfgs_min.x

        if not bfgs_min.success:
            logging.warning(f"One of the minimum finding points failed at l {l} and cut index {cut_index} and redchi {bfgs_min.fun / data_cuts[l-2][cut_index] / 3}.")
            logging.error(bfgs_min)
            return None, None, None, None

        end_redchi = 2 * bfgs_min.fun / data_cuts[l-2][cut_index] / 3
        if abs(start_redchi - end_redchi) < 1e-3 and start_redchi / true_redchi > THRESHOLD_LIKE_RAT:
            logging.warning(f"At cut index {cut_index}, the redchi {start_redchi} was unchanged by the fit.")
            logging.debug(bfgs_min)
        else:
            logging.debug(f"Minimization index {cut_index} passed successfully (redchi {start_redchi} -> {-end_redchi})")

        redchi_record.append((start_redchi, end_redchi))

    # Now do fit on full data

    def minimize_log_prob_uncut(x):
        return minimize_log_prob(y, x, fix_theta, simulate_func, l-2, -1)

    result = bfgs_min.x
    
    # Now that an answer has been achieved, find the hessian

    grad = nd.Gradient(minimize_log_prob_uncut)(result)
    hess = nd.Hessian(minimize_log_prob_uncut)(result)
    minimizing_likelihood = bfgs_min.fun
    #grad = derivatives.Gradient(hess_func, step=1e-10)(result)
    #hess = derivatives.Hessian(hess_func, step=1e-10)(result)

    if l == 2:
        logging.debug(minimize_log_prob_uncut(result+np.array([1.0e-8, 0, 0]))-minimizing_likelihood)
        logging.debug(minimize_log_prob_uncut(result-np.array([1.0e-8, 0, 0]))-minimizing_likelihood)
        logging.debug(minimize_log_prob_uncut(result+np.array([0, 1.0e-8, 0]))-minimizing_likelihood)
        logging.debug(minimize_log_prob_uncut(result-np.array([0, 1.0e-8, 0]))-minimizing_likelihood)
        logging.debug(minimize_log_prob_uncut(result+np.array([0, 0, 1.0e-8]))-minimizing_likelihood)
        logging.debug(minimize_log_prob_uncut(result-np.array([0, 0, 1.0e-8]))-minimizing_likelihood)

    logging.debug("Hessian: {}".format(hess))

    logging.info(f"REDCHI RECORD: {redchi_record}, THETA: {result}")
    logging.debug(f"GRADIENT: {grad}")
    logging.debug(f"THETA: {result}")

    new_evals = []
    new_evecs = []
    try:
        evals, evecs = mp.eigsy(mp.matrix(hess))
    except Exception:
        logging.error(f"Could not find eigenvectors of the Hessian {hess}")
        return None, None, None, None
    for e in evals:
        ### Correct for non positive definite hessians
        new_evals.append(float(e))
    new_evals = np.array(new_evals)
    
    if np.any(np.asarray(new_evals) < 0.0):
        logging.warning(f"The Hessian was not positive definite. \nHessian {hess}\nRedchi {2 * minimizing_likelihood / len(y) / 3}, \nTheta {result}, \nEigenvalues {new_evals}, \nGradient {grad}")
        #logging.debug(bfgs_min)
        if l == 2:
            #return None, None, None, None
            new_evals[0] = -new_evals[0]
            if np.any(np.asarray(new_evals) < 0.0):
                return None, None, None, None
            else:
                logging.warning("Allowing the Hessian to pass with reversed first eigenvalue.")

        if l == 3:
            new_evals[3:] = 1 / LARGE_NUMBER
            logging.warning("Allowing the Hessian to pass with maxed eigenvalues for l=3 components.")
            
    logging.debug(f"Eigenvalues: {new_evals}")

    for k in range(int(len(evecs))):
        new_evecs.append(np.array([evecs[j, k] for j in range(int(len(evecs)))],
        dtype=np.float64))
    new_evecs = np.array(new_evecs)

    #logging.info("Reconstructed Hessian:", np.matmul(new_evecs.transpose(), np.matmul(np.diag(new_evals), new_evecs)))

    # Return redchi, minimizing params, hessian
    return 2 * minimizing_likelihood / len(y) / 3, result, new_evals, new_evecs


def minimize(y, l, fix_theta):
    assert l <= ASTEROIDS_MAX_K
    assert len(fix_theta) == max((l)**2 - 6, 0)

    if l == 2:
        if ASTEROIDS_MAX_J == 0:
            simulate = asteroids_0_2.simulate
        elif ASTEROIDS_MAX_J == 2:
            simulate = asteroids_2_2.simulate
        elif ASTEROIDS_MAX_J == 3:
            simulate = asteroids_3_2.simulate
    elif l == 3:
        if ASTEROIDS_MAX_J == 0:
            simulate = asteroids_0_3.simulate
        elif ASTEROIDS_MAX_J == 2:
            simulate = asteroids_2_3.simulate
        elif ASTEROIDS_MAX_J == 3:
            simulate = asteroids_3_3.simulate
    else:
        raise Exception("l={} is not supported".format(l))

    theta_indices = range(max(0, (l)**2 - 6), (l+1)**2 - 6)
    bounds = np.asarray([(theta_low[i], theta_high[i]) for i in theta_indices])
    bound_widths = np.asarray([theta_high[i] - theta_low[i] for i in theta_indices])

    # Choose the seed parameters
    parameter_points = []

    while len(parameter_points) < NUM_MINIMIZE_POINTS_PER * len(theta_true):
        randoms = np.asarray([random.random() for i in theta_indices])
        point = np.asarray(theta_low)[theta_indices] + bound_widths * randoms
        try:
            model = minimize_function(point, simulate, 0, 0) # Terminate the run early with drc=2 because I just want the runtime error
        except RuntimeError:
            continue
        parameter_points.append([y, point, fix_theta, l, bounds])

    #parameter_points[0][0] = np.array(theta_true) + np.random.randn(3)*0.001

    # Perform the minimization
    with Pool() as pool:
        results = pool.map(get_minimum, parameter_points)

    num_succeeded = 0
    for a, _, _, _ in results:
        if a is not None:
            num_succeeded += 1

    logging.info(f">>>>>>>>> Of {NUM_MINIMIZE_POINTS_PER * len(theta_true)} minimizations run, {num_succeeded} succeeded. <<<<<<<<<<")

    # Extract the lowest redchis
    sorted_results = sorted(results, key=
        lambda x: x[0] if (x[0] is not None and not np.isnan(x[0])) else LARGE_NUMBER)

    distinct_results = []
    for redchi, theta, evals, evecs in sorted_results:
        if redchi is None:
            continue
        if redchi / true_redchi > THRESHOLD_LIKE_RAT: # divide by 2 because redchi is actually 0.5 of true redchi.
            continue
        choose = True
        for distinct_theta, _, _, _ in distinct_results:
            if np.all(np.abs(distinct_theta - theta)[1:3] < MIN_THETA_DIST):
                # Only scan for theta 1 and 2 (K2m)
                choose = False
                break
            if len(theta) > 3 and len (distinct_results) >= NUM_L3_MINIMIZATIONS:
                # Do not choose more than NUM_L3_MINIMIZATIONS l=3 tiers.
                choose = False
                break

        if choose:
            distinct_results.append((theta, evals, evecs, redchi))
    return distinct_results

def populate(evals, diagonalizer, count, start):
    spacing = 1 / np.sqrt(evals)
    if np.any(spacing < MIN_SPACING):
        logging.error(f"Had to widen some sigmas. Sigmas were {spacing}")
        spacing = np.maximum(MIN_SPACING * SIGMA_FACTOR, spacing)
    if not np.all(np.isfinite(spacing)):
        logging.error(f"Some sigmas were inf. Sigmas were {spacing}")
        spacing[~np.isfinite(spacing)] = MAX_SPACING[~np.isfinite(spacing)] * SIGMA_FACTOR
    if np.any(spacing > 1):
        logging.error(f"Some sigmas were greater than one. Sigmas were {spacing}")
        spacing[spacing > 1] = 1.0
    logging.info(f"Sigmas: {spacing}")
    logging.info(f"Evecs: {diagonalizer}")

    diagonal_points = spacing * (np.random.randn(count * N_DIM).reshape(count, N_DIM))
    global_points = np.asarray([np.matmul(diagonalizer.transpose(), d) for d in diagonal_points]) + start

    for i, point in enumerate(global_points):
        while True:
            logprob = log_probability(point, y, sigma)
            if np.isfinite(logprob):
                break
            point = np.matmul(diagonalizer.transpose(), spacing * np.random.randn(N_DIM)) + start
            global_points[i] = point
        logging.info("{}, {}".format(point, logprob / len(y) / 3))

    return global_points

def mcmc_fit(theta_start, evals, evecs, index):
    backend = emcee.backends.HDFBackend(output_name+"-{}.h5".format(index))

    if not reload:
        pos = populate(evals, evecs, N_WALKERS, theta_start)


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
        logging.info("Initial size: {}".format(backend.iteration))

    old_tau = np.inf
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(N_WALKERS, N_DIM, log_probability,
            args=(y, sigma), backend=backend, pool=pool)

        if reload:
            pos = sampler._previous_state

        if not emcee.walkers_independent(pos):
            f = open("errors.dat", 'a')
            f.write(output_name + f": walkers weren't independent (eigenvalues {evals})\n")
            f.close()
            logging.warning(f"Walkers weren't independent. Eigenvalues {evals}")
            sys.exit()

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
                    logging.info("Converged")
                    break
                old_tau = tau

            if sampler.iteration % (MAX_N_STEPS//10) == 0:
                redchis = -2 * sample.log_prob / len(y) / 3
                logging.info(f"Minimum redchi at {sampler.iteration/MAX_N_STEPS * 100}\%: {np.min(redchis)}")
        sampler._previous_state = sample

    if reload:
        logging.info("New size: {}".format(backend.iteration))
    else:
        logging.info("Done")

if __name__ == "__main__":
    assert(len(theta_true) == len(theta_high) == len(theta_low))
    assert(len(theta_true) == (ASTEROIDS_MAX_K + 1)**2 - 6)
    assert(len(jlms) == (ASTEROIDS_MAX_J + 1)**2)
    assert(np.all(theta_high > theta_low))

    logging.info("Cadence {}, perigee {}, speed {}".format(cadence, perigee, speed))
    logging.info("Spin {}".format(spin))
    logging.info("Jlms {}".format(jlms))
    logging.info("Radius {}".format(radius))
    logging.info("Theta true {}".format(theta_true))
    logging.info("Theta high {}".format(theta_high))
    logging.info("Theta low {}".format(theta_low))
    logging.info("Sigma {}".format(sigma))
    logging.info("Name {}".format(output_name))
    logging.info("Velocity multiplier {}".format(VELOCITY_MUL))
    ####################################################################
    # Generate synthetic data
    ####################################################################
    start = time.time()
    y = fit_function(theta_true)
    logging.info("Data generation took {} s".format(time.time() - start))
    y, y_inv_covs = random_vector.randomize(UNCERTAINTY_MODEL, y, sigma)
    if UNCERTAINTY_MODEL == random_vector.TILT_UNIFORM_TRUE:
        UNCERTAINTY_ARGUMENT = sigma
    else:
        UNCERTAINTY_ARGUMENT = y_inv_covs

    np.save(f"{output_name}-data.npy", y)
    np.save(f"{output_name}-unc.npy", y_inv_covs)

    logging.info(f"DOF: {len(y)}")

    plt.figure(figsize=(12, 4))
    x_display = np.arange(len(y))
    uncs = np.array([2 * np.sqrt(np.diagonal(pinvh(a))) for a in y_inv_covs])
    plt.fill_between(x_display, y[:,0]+uncs[:,0], y[:,0]-uncs[:,0], alpha=0.5, color="C0")
    plt.fill_between(x_display, y[:,1]+uncs[:,1], y[:,1]-uncs[:,1],alpha=0.5, color="C1")
    plt.fill_between(x_display, y[:,2]+uncs[:,2], y[:,2]-uncs[:,2],  alpha=0.5, color="C2")
    plt.scatter(x_display, y[:,0], label='x', s=1, color="C0")
    plt.scatter(x_display, y[:,1], label='y', s=1, color="C1")
    plt.scatter(x_display, y[:,2], label='z', s=1, color="C2")
    plt.xlabel("Time (Cadences)")
    plt.ylabel("Spin (rad/s)")
    for tier in data_cuts:
        for drc in tier:
            plt.axvline(x=drc, color='k')
    plt.legend()
    # plt.show()



    ####################################################################
    # Minimize likelihood
    ####################################################################

    if ASTEROIDS_MAX_K == 3:
        if ASTEROIDS_MAX_J == 0:
            real_sim_func = asteroids_0_3.simulate
        elif ASTEROIDS_MAX_J == 2:
            real_sim_func = asteroids_2_3.simulate
        elif ASTEROIDS_MAX_J == 3:
            real_sim_func = asteroids_3_3.simulate
    elif ASTEROIDS_MAX_K == 2:
        if ASTEROIDS_MAX_J == 0:
            real_sim_func = asteroids_0_2.simulate
        elif ASTEROIDS_MAX_J == 2:
            real_sim_func = asteroids_2_2.simulate
        elif ASTEROIDS_MAX_J == 3:
            real_sim_func = asteroids_3_2.simulate


    true_redchi = 2 * minimize_log_prob(y, theta_true, [], real_sim_func, 0, -1) / len(y) / 3
    logging.info("TRUE REDCHI: {}".format(true_redchi))

    # Stepped minimization
    queue = [([], [], [], 0)]
    for i in range(2, ASTEROIDS_MAX_K + 1):
        new_queue = []
        for fix_theta, evals, evecs, _ in queue:
            tier_results = minimize(y, i, fix_theta)
            for result_theta, result_evals, result_evecs, result_redchi in tier_results:
                logging.info("Deg {} redchi: {}. Theta: {}".format(i, result_redchi, result_theta))
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
        if redchi / true_redchi > THRESHOLD_LIKE_RAT:
            logging.warning(f"The point with theta {theta} and redchi {redchi} was excluded from the kernel")
            continue
        logging.info("The kernel includes a point with theta {}, redchi {}".format(theta, redchi))
        kernel.append((theta, evals, np.array(resized_evecs)))

    logging.info("There are {} MCMC starting points".format(len(kernel)))

    if len(kernel) == 0:
        f = open("errors.dat", 'a')
        f.write(output_name+": kernel is empty\n")
        f.close()
        logging.critical("Kernel is empty")
        sys.exit()

    ####################################################################
    # Run MCMC
    ####################################################################


    for i, (theta, evals, evecs) in enumerate(kernel):
        mcmc_fit(theta, evals, evecs, i)