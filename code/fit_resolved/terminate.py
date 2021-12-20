TEST = False

import numpy as np
import matplotlib.pyplot as plt
import emcee, time, sys, os
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


GM = 3.986004418e14
EARTH_RADIUS = 6_370_000

def terminal(output_name):

    f = open("../../staged/" + output_name+".txt", 'r')
    f.readline()
    num_trials = np.prod([int(i) for i in f.readline().split(', ')])
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

    N_DIM = len(theta_true)

    ####################################################################
    # Process data
    ####################################################################
    print()
    print(output_name.upper())

    i = 0
    while True:
        try:
            disp = display.Display("{0}".format(output_name), "{0}-{1}".format(output_name, i))
        except Exception as e:
            if i == 0:
                print(f"Run {output_name} failed.")
                return
                #raise e
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

    for index in range(num_trials):
        reader = emcee.backends.HDFBackend(output_name+"-{}.h5".format(index), read_only=True)

        try:
            tau = reader.get_autocorr_time()
            burnin = int(2 * np.max(tau))
            thin = int(0.5 * np.min(tau))
        except:
            print("Could not find autocorrelation time because the chain is too short.")
            thin = 1
            burnin = 1000

        samples = reader.get_chain(discard=burnin, thin=thin)
        #log_prob_samples = reader.get_log_prob(discard=burnin, thin=thin)
        #log_prior_samples = reader.get_blobs(discard=burnin, thin=thin)

        flat_samples = np.array([samples[:,:,i].flatten() for i in range(len(theta_true))])

        np.savetxt(output_name+"-{}-samples.dat".format(index), flat_samples)

if __name__ == "__main__":
    if len(sys.argv) == 2:
        output_name = sys.argv[1]
        terminal(output_name)

    else:
        for name in os.listdir('../../staged'):
            terminal(name[:-4])
            plt.close()
            plt.cla() # Close the figure in multiple ways to prevent memory overflow
            plt.clf()
