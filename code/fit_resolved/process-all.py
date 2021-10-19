import numpy as np
import matplotlib.pyplot as plt
import emcee, time, sys, os
from multiprocessing import Pool
from random_vector import *
from scipy import optimize, linalg
from mpmath import mp
import plotille
import display, collect

import numdifftools as nd

EARTH_RADIUS = 6_370_000

def process(output_name):
    f = open("../../staged/" + output_name+".txt", 'r')
    ASTEROIDS_MAX_J, ASTEROIDS_MAX_K = f.readline().split(', ')
    ASTEROIDS_MAX_J = int(ASTEROIDS_MAX_J)
    ASTEROIDS_MAX_K = int(ASTEROIDS_MAX_K)

    NUM_FITS = f.readline().split(', ')
    NUM_FITS = [int(i) for i in NUM_FITS]

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

    ####################################################################
    # Process data
    ####################################################################
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

    ####################################################################
    # Save samples
    ####################################################################

    index = 0
    while True:
        try:
            reader = emcee.backends.HDFBackend(output_name+"-{}.h5".format(index), read_only=True)
        except:
            print("Max index was {} (no h5)".format(index-1))
            break

        try:
            tau = reader.get_autocorr_time()
            burnin = int(2 * np.max(tau))
            thin = int(0.5 * np.min(tau))
        except:
            print("Could not find autocorrelation time because the chain is too short.")
            thin = 1
            burnin = 1000

        try:
            samples = reader.get_chain(discard=burnin, thin=thin)
        except:
            print("Max index was {} (broken h5)".format(index-1))
            break
        #log_prob_samples = reader.get_log_prob(discard=burnin, thin=thin)
        #log_prior_samples = reader.get_blobs(discard=burnin, thin=thin)


        flat_samples = np.array([samples[:,:,i].flatten() for i in range(3)])
        np.savetxt(output_name+"-{}-samples.dat".format(index), flat_samples)
        index += 1


if len(sys.argv) not in [2, 3]:
    for item in os.listdir("../../staged/"):
        try:
            process(item[:-4])
        except:
            print("{} failed".format(item))
else:
    process(sys.argv[1])
