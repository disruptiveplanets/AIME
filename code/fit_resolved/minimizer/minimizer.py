import numpy as np
import matplotlib.pyplot as plt
import sys, os, inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

import asteroids, time
from multiprocessing import Pool
from random_vector import *
from scipy import optimize, linalg

np.random.seed(123)


ASTEROIDS_MAX_K = 2 # Remember to change the counterpart in backend.hpp
ASTEROIDS_MAX_J = 0 # Remember to change the counterpart in backend.hpp
EARTH_RADIUS = 6370000
NUM_MINIMIZE_POINTS = 8
EPSILON = 1e-11

CADENCE_CUT = 650


if len(sys.argv) not in [2, 3]:
    raise Exception("Please pass a file to describe the fit")
output_name = sys.argv[1]
f = open("../../../staged/" + output_name+".dat", 'r')
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

N_DIM = len(theta_true)

reload = False
if len(sys.argv) == 3 and sys.argv[2] == "reload":
    reload = True
    REGENERATE_DATA = False

def fit_function(theta):
    resolved_data = asteroids.simulate(cadence, jlms, theta[1:], radius,
        spin[0], spin[1], spin[2], theta[0], impact_parameter, speed, CADENCE_CUT)
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

plt.figure(figsize=(12, 4))
x_display = np.arange(len(y) / 3)
plt.errorbar(x_display, y[::3], yerr=yerr[::3], label = 'x', fmt='.')
plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], label = 'y', fmt='.')
plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], label = 'z', fmt='.')
plt.xlabel("Time (Cadences)")
plt.ylabel("Spin (rad/s)")
plt.axvline(x=CADENCE_CUT, color='k')
plt.legend()
plt.show()

####################################################################
# Minimize likelihood
####################################################################
print()
def redchi(theta, y, yerr):
    # Normal likelihood
    try:
        model = fit_function(theta)
    except RuntimeError:
        return 1e10 # Zero likelihood
    return np.sum((y - model) ** 2 /  yerr ** 2) / len(y)

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
    bfgs_min = optimize.minimize(redchi, point, args=(y, yerr), method='L-BFGS-B', options={"eps": EPSILON}, bounds=bounds)
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

min_like = None
min_theta = None
min_hessian = None
for like, theta, hess in results:
    if like is None: continue
    if min_like is None or like < min_like:
        min_like = like
        min_theta = theta
        min_hessian = hess

print("Redchi was {} at parameters {}".format(min_like / len(y), min_theta))

print("All red chis:")
for like, theta, _ in results:
    if like is None: continue
    print(like, "at parameters", theta)

####################################################################
# Plot minimizing locations:
####################################################################

for file in os.listdir():
    # Generate these dat files with visualize-likelihood/populator.py
    if file[-4:] != ".dat":
        continue
    f = open(file, 'r')
    lines = f.readlines()
    f.close()
    index_x, index_y = lines[0].split(', ')
    index_x = int(index_x)
    index_y = int(index_y)

    theta_true = []
    for t in lines[1].split(', '):
        theta_true.append(float(t))

    xs = []
    for x in lines[2].split(', '):
        xs.append(float(x))

    ys = []
    for y in lines[3].split(', '):
        ys.append(float(y))

    red_chis = []
    for line in lines[4:]:
        red_line = []
        for rc in line.split(', '):
            red_line.append(float(rc))
        red_chis.append(red_line)

    finites = np.asarray(red_chis).reshape(len(red_chis) * len(red_chis[1]))
    finites = finites[np.isfinite(finites)]

    plt.figure()
    c = plt.pcolor(xs, ys, red_chis, vmax=np.nanpercentile(finites, 80), vmin=np.nanpercentile(finites, 1))
    plt.colorbar(c)
    plt.xlabel("Theta {}".format(index_x))
    plt.ylabel("Theta {}".format(index_y))
    thetas_x = []
    thetas_y = []
    for i, (_, theta, _) in enumerate(results):
        if theta is None: continue
        plt.plot([parameter_points[i][index_x], theta[index_x]], [parameter_points[i][index_y], theta[index_y]], linestyle='dashed', color="C1")
        thetas_x.append(theta[index_x])
        thetas_y.append(theta[index_y])
    plt.scatter(thetas_x, thetas_y,  color='C1', s=12, label="end points")
    plt.scatter(min_theta[index_x], min_theta[index_y], color='red', s=12, label="Minimum chisq")
    plt.plot(theta_true[index_x], theta_true[index_y], marker='*', markerfacecolor='red', markersize=8, markeredgecolor='black')
    plt.legend()
    plt.savefig("theta-{}-{}.png".format(index_x, index_y))
plt.show()
