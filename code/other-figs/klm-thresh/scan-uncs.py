# PROBLEM: this data doesn't have K3m.

import sys, os
import numpy as np
sys.path.append("../../density")
from likelihood import Likelihood as Model
#from likelihood import Harmonic as Model
from core import Asteroid, Indicator, TrueShape

DIVISION = 99
NUM_SAMPLES = 153

sigmas = np.zeros((NUM_SAMPLES, 10))

percentiles = {}
with open("shape-data/percentiles.dat", 'r') as f:
    for line in f.readlines():
        if line == '': continue
        elements = line.split(':')
        name = elements[0]
        perc_array = []
        for percs in elements[1:]:
            perc_array.append([float(x) for x in percs.split(',')])
        percentiles[name] = np.array(perc_array)

for i in range(NUM_SAMPLES):
    d = "param-{:03}".format(i)
    if not os.path.exists(f"shape-data/{d}/{d}-0-samples.npy"):
        continue
    with open(f"shape-data/{d}/{d}.txt", 'r') as f:
        f.readline()
        cadence = int(float(f.readline()))
        perigee = float(f.readline()) # In Earth radii
        radius = float(f.readline())
        speed = float(f.readline())
        spin = [float(x) for x in f.readline().split(',')]
        jlms = [float(x) for x in f.readline().split(',')]
        theta_true = [float(x) for x in f.readline().split(',')]
        theta_high = np.asarray([float(x) for x in f.readline().split(',')])
        theta_low = np.asarray([float(x) for x in f.readline().split(',')])
        sigma = [float(d) for d in f.readline().split(", ")]# theta, ratio
        last_line = f.readline()

    print("ASTEROID", i)

    am = radius
    k20 = theta_true[2]
    k22 = theta_true[1]

    a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22)
    b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22)
    c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20)

    max_radius = int(max(a, b, c) + 4 * DIVISION)
    
    asteroid = Asteroid("shape-{i}", f"shape-data/{d}/{d}-0-samples.npy", am, DIVISION, max_radius, Indicator.ell(radius, k22, k20), None)
    method = Model(asteroid)
    method.solve()
    uncs = method.map_unc()

    density_uncertainty = np.nanmedian(uncs)
    
    ## Divide original sigmas by density_uncertainty
    ptiles = percentiles[f"{d}-0-samples.npy"]
    original_sigma = ((ptiles[:,2] - ptiles[:,3]) - (ptiles[:,-2] - ptiles[:,3]))/2

    print(density_uncertainty)

    sigmas[i] = original_sigma / density_uncertainty

with open("sigmas.npy", 'wb') as f:
    np.save(f, sigmas)

print("Done")