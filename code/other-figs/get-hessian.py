
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool
from scipy import optimize, linalg
from mpmath import mp
import numdifftools as nd
from shutil import copyfile
copyfile("../fit_resolved/asteroids_0_2.cpython-38-x86_64-linux-gnu.so", "asteroids_0_2.cpython-38-x86_64-linux-gnu.so")
import asteroids_0_2, sys
sys.path.insert(0, "../fit_resolved")
from random_vector import randomize_rotate


THETA = [0.39269908169, 0.05200629, -0.2021978]
SPIN = [0.0001, 0.0002, 0.0003]
IMPACT_PARAMETER = 5 * 6_378_000
SPEED = 6000
RADIUS = 1000
JLMS = [1]
CADENCE = 120
MIN_SPREAD = 1e-4 ** 2
AMPLIFY = 200
EPSILON = 1e-12

SIGMA = 0.1

def fit_function(theta):
    resolved_data = asteroids_0_2.simulate(CADENCE, JLMS, theta[1:], RADIUS,
        SPIN[0], SPIN[1], SPIN[2], theta[0], IMPACT_PARAMETER, SPEED, -1)
    return np.asarray(resolved_data)

####################################################################
# Minimize likelihood
####################################################################
y_min = fit_function(THETA)
_, yerr_min = randomize_rotate(y_min, SIGMA)

def minimize_log_prob(theta):
    # Normal likelihood
    try:
        model = fit_function(theta)
    except RuntimeError:
        return 1e10 # Zero likelihood

    redchi = np.sum((y_min - model) ** 2 /  yerr_min ** 2) / len(y_min)
    return redchi


grad = nd.Gradient(minimize_log_prob, step=EPSILON)(THETA)
print(grad)

hess = nd.Hessian(minimize_log_prob, step=EPSILON)(THETA)
print(hess)
