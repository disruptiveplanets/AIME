TEST = False

import numpy as np
import matplotlib.pyplot as plt
import emcee, time, sys
if not TEST:
    import asteroids_0_2, asteroids_0_3, asteroids_2_2, asteroids_2_3, asteroids_3_2, asteroids_3_3
if TEST:
    import test_0_2 as asteroids_0_2
    import test_0_2 as asteroids_0_3
    import test_0_2 as asteroids_2_2
    import test_0_2 as asteroids_2_3
from multiprocessing import Pool
import random_vector, random
from scipy import optimize, linalg
from mpmath import mp
import plotille
import display, collect

import numdifftools as nd

TEST_THETA = [
    [0.54744608, 0.04032506, -0.16942554],
    [0.54740554, 0.04017766, -0.16945756],
    [0.54746540, 0.04038743, -0.16948097],
    [0.54758347, 0.04047878, -0.16934459],
    [0.54743905, 0.04029860, -0.16938477],
    [0.54743763, 0.04029555, -0.16946423],
    [0.54766856, 0.04038121, -0.1694405 ],
    [0.54748953, 0.04004225, -0.16952856],
    [0.54738075, 0.04031553, -0.1694349 ],
    [0.54748118, 0.04027882, -0.16931596],
    [0.54761604, 0.04013331, -0.169445  ],
    [0.54749511, 0.04020166, -0.16946253],
    [0.54733365, 0.04035511, -0.16951817],
    [0.54746701, 0.04038966, -0.16952725],
    [0.54760029, 0.04019186, -0.16932237],
    [0.54742583, 0.04016524, -0.16932933],
    [0.54743737, 0.04040893, -0.16957725],
    [0.54739316, 0.04019740, -0.16936259],
    [0.54734638, 0.04034705, -0.16942803],
    [0.54745556, 0.04025321, -0.16932022],
    [0.54753959, 0.04023222, -0.16943905],
    [0.54730535, 0.04027055, -0.16967787],
    [0.54736582, 0.04020499, -0.16946332],
    [0.54734158, 0.04023380, -0.16957096],
    [0.54746875, 0.04029686, -0.16945413],
    [0.54743498, 0.04041504, -0.16944435],
    [0.54741520, 0.04025687, -0.16943495],
    [0.54750947, 0.04056622, -0.16948461],
    [0.54717167, 0.04038458, -0.16931208],
    [0.54750349, 0.04044016, -0.16922549],
    [0.54728596, 0.04013732, -0.16942784],
    [0.54722685, 0.04026163, -0.16944982],]


GM = 3.986004418e14
UNCERTAINTY_MODEL = 1

INTEGRAL_LIMIT_FRAC = 1.0e-3

EARTH_RADIUS = 6_370_000
N_WALKERS = 32
MAX_N_STEPS = 10_000
NUM_MINIMIZE_POINTS = 48
EPSILON = 1e-10
MIN_SPREAD = 1e-4 ** 2

DISTANCE_RATIO_CUT = 10
MIN_THETA_DIST = 0.01
LARGE_NUMBER = 1e100

if len(sys.argv) not in [2, 3]:
    raise Exception("Please pass a file to describe the fit")
output_name = sys.argv[1]
f = open("../../staged/" + output_name+".txt", 'r')
ASTEROIDS_MAX_J, ASTEROIDS_MAX_K = f.readline().split(', ')
ASTEROIDS_MAX_J = int(ASTEROIDS_MAX_J)
ASTEROIDS_MAX_K = int(ASTEROIDS_MAX_K)

NUM_FITS = f.readline().split(', ')
NUM_FITS = [int(i) for i in NUM_FITS]

cadence = int(f.readline())
perigee = EARTH_RADIUS * float(f.readline())
radius = float(f.readline())
speed = float(f.readline())
spin = [float(x) for x in f.readline().split(',')]
jlms = [float(x) for x in f.readline().split(',')]
theta_true = [float(x) for x in f.readline().split(',')]
theta_high = np.asarray([float(x) for x in f.readline().split(',')])
theta_low = np.asarray([float(x) for x in f.readline().split(',')])

sigma = float(f.readline())
while output_name[-1] == '\n':
    output_name = output_name[:-1]
f.close()
assert(len(theta_true) == len(theta_high) == len(theta_low))
assert(len(theta_true) == (ASTEROIDS_MAX_K + 1)**2 - 6)
assert(len(jlms) == (ASTEROIDS_MAX_J + 1)**2)
assert(np.all(theta_high > theta_low))

NUM_FITS = NUM_FITS[:ASTEROIDS_MAX_K-1]

N_DIM = len(theta_true)

def fit_function(theta, target_length=None):
    if ASTEROIDS_MAX_K == 3:
        if ASTEROIDS_MAX_J == 0:
            resolved_data = asteroids_0_3.simulate(cadence, jlms, theta[1:], radius,
                spin[0], spin[1], spin[2], theta[0], perigee, speed, GM, EARTH_RADIUS, -1)
        elif ASTEROIDS_MAX_J == 2:
            resolved_data = asteroids_2_3.simulate(cadence, jlms, theta[1:], radius,
                spin[0], spin[1], spin[2], theta[0], perigee, speed, GM, EARTH_RADIUS, -1)
        elif ASTEROIDS_MAX_J == 3:
            resolved_data = asteroids_3_3.simulate(cadence, jlms, theta[1:], radius,
                spin[0], spin[1], spin[2], theta[0], perigee, speed, GM, EARTH_RADIUS, -1)
    elif ASTEROIDS_MAX_K == 2:
        if ASTEROIDS_MAX_J == 0:
            resolved_data = asteroids_0_2.simulate(cadence, jlms, theta[1:], radius,
                spin[0], spin[1], spin[2], theta[0], perigee, speed, GM, EARTH_RADIUS, -1)
        elif ASTEROIDS_MAX_J == 2:
            resolved_data = asteroids_2_2.simulate(cadence, jlms, theta[1:], radius,
                spin[0], spin[1], spin[2], theta[0], perigee, speed, GM, EARTH_RADIUS, -1)
        elif ASTEROIDS_MAX_J == 3:
            resolved_data = asteroids_3_2.simulate(cadence, jlms, theta[1:], radius,
                spin[0], spin[1], spin[2], theta[0], perigee, speed, GM, EARTH_RADIUS, -1)
    if target_length is not None:
        while len(resolved_data)//3 < target_length:
            resolved_data.append(resolved_data[-3])
            resolved_data.append(resolved_data[-3])
            resolved_data.append(resolved_data[-3])
    return np.asarray(resolved_data).reshape(-1, 3)


def log_likelihood(theta, y, y_inv_covs):
    # Normal likelihood
    try:
        model = fit_function(theta, len(y))
    except RuntimeError:
        return -np.inf # Zero likelihood
    # Return chisq
    chisq = 0
    for i in range(len(y)):
        chisq += np.matmul(y[i] - model[i], np.matmul(y_inv_covs[i], y[i] - model[i]))
    return -chisq

def log_prior(theta):
    for i, param in enumerate(theta):
        if param > theta_high[i] or param < theta_low[i]:
            return -np.inf
    return 0.0

def log_probability(theta, y, y_inv_covs):
    prior = log_prior(theta)
    like = log_likelihood(theta, y, y_inv_covs)
    if not np.isfinite(prior) or not np.isfinite(like):
        return -np.inf
    return prior + like

####################################################################
# Generate synthetic data
####################################################################
start = time.time()
y = fit_function(theta_true)
print("Data generation took {} s".format(time.time() - start))
if UNCERTAINTY_MODEL == 1:
    y, y_inv_covs = random_vector.randomize_rotate_uniform(y, sigma)
if UNCERTAINTY_MODEL == 2:
    raise Exception("Unimplemented")

print("DOF:", len(y))


for theta in TEST_THETA:
    print(log_probability(theta, y, y_inv_covs))