from scipy.optimize import minimize
import time, sys, os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir + "/fit_resolved/")
import numpy as np
import matplotlib.pyplot as plt
import asteroids
from random_vector import *


ASTEROIDS_MAX_K = 2 # Remember to change the counterpart in backend.hpp
ASTEROIDS_MAX_J = 0 # Remember to change the counterpart in backend.hpp
EARTH_RADIUS = 6370000

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
theta_start = [float(x) for x in f.readline().split(',')]
theta_spread = [float(x) for x in f.readline().split(',')]
theta_high = np.asarray([float(x) for x in f.readline().split(',')])
theta_low = np.asarray([float(x) for x in f.readline().split(',')])

sigma = float(f.readline()) * np.sqrt(spin[0]**2 + spin[1]**2 + spin[2]**2)
while output_name[-1] == '\n':
    output_name = output_name[:-1]
f.close()
assert(len(theta_true) == len(theta_start) == len(theta_spread) == len(theta_high) == len(theta_low))
assert(len(theta_true) == (ASTEROIDS_MAX_K + 1)**2 - 6)
assert(len(jlms) == (ASTEROIDS_MAX_J + 1)**2)
assert(np.all(theta_high > theta_low))
bounds = [(theta_low[i], theta_high[i]) for i in range(len(theta_high))]

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

def fit_function(theta):
    resolved_data = asteroids.simulate(cadence, jlms, theta[1:], radius,
        spin[0], spin[1], spin[2], theta[0], impact_parameter, speed)
    return np.asarray(resolved_data)

def red_chi(theta, y, yerr):
    # Normal likelihood
    try:
        model = fit_function(theta)
    except RuntimeError:
        return 1e10 / len(y) # Zero likelihood
    return np.sum((y - model) ** 2 /  yerr ** 2) / len(y)

start = time.time()
y = fit_function(theta_true)
print("Data generation took {} s".format(time.time() - start))
x = np.arange(len(y))
y, yerr = randomize_rotate(y, sigma)


bfgs_min = minimize(red_chi, theta_start, args=(y, yerr), method='L-BFGS-B', bounds=bounds, options={"eps": 1e-10})
print(bfgs_min)
print(bfgs_min.hess_inv)
print(bfgs_min.hess_inv.H)

#bfgs_real_min = minimize(minus_log_likelihood, theta_true, args=(y, yerr), method='L-BFGS-B', bounds=bounds, options={"eps": 1e-10})
#print(bfgs_min)
