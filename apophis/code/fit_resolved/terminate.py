import numpy as np
import matplotlib.pyplot as plt
import emcee, sys, os
from multiprocessing import Pool
import display, collect

import numdifftools as nd


THRESHOLD_REDCHI = 2
DEFAULT_THIN = 10

def terminal(output_name, do_not_duplicate=True):

    f = open("../../staged/" + output_name+".txt", 'r')
    f.readline()
    cadence = int(float(f.readline()))
    radius = float(f.readline())
    theta_true = [float(x) for x in f.readline().split(',')]
    theta_high = np.asarray([float(x) for x in f.readline().split(',')])
    theta_low = np.asarray([float(x) for x in f.readline().split(',')])

    sigma = [float(d) for d in f.readline().split(',')]
    while output_name[-1] == '\n':
        output_name = output_name[:-1]
    f.close()

    N_DIM = len(theta_true)

    # Make sure this file hasn't been done before
    if do_not_duplicate:
        already_done = True
        for index in range(10000):
            if not os.path.exists(f'{output_name}-{index}-samples.npy'):
                if os.path.exists(f'{output_name}-{index}.h5'):
                    # H5 exists, but samples does not exist
                    already_done = False
                    break
                else:
                    # Neither exists, so this index was not reached.
                    break
        if already_done:
            plt.cla()
            plt.clf()
            plt.close('all')
            return
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
                plt.cla()
                plt.clf()
                plt.close('all')
                return
                #raise e
            break
        disp.show_redchi()
        disp.show_params()
        n_data = disp.show_corner()
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

    for index in range(1000):
        reader = emcee.backends.HDFBackend(output_name+"-{}.h5".format(index), read_only=True)


        try:
            tau = reader.get_autocorr_time()
            burnin = int(2 * np.max(tau))
            thin = DEFAULT_THIN#int(0.5 * np.min(tau))
        except:
            print("Could not find autocorrelation time because the chain is too short.")
            thin = DEFAULT_THIN
            burnin = 1000
        try:
            samples = reader.get_chain(discard=burnin, thin=thin)
        except Exception:
            break

        log_prob_samples = reader.get_log_prob(discard=burnin, thin=thin) / 3 / n_data 
        #log_prior_samples = reader.get_blobs(discard=burnin, thin=thin)
        redchis = -2 * np.nanmin(log_prob_samples, axis=0)
        
        mask = np.where(redchis < THRESHOLD_REDCHI)[0]

        print(f"Saving samples from {len(mask)}/{len(redchis)} walkers")

        with open(f"{output_name}-{index}-samples.npy", 'wb') as f:
            np.save(f, samples[:,mask,:], allow_pickle=False)

    plt.cla()
    plt.clf()
    plt.close('all')

def wrap_terminal(name, b=True): # False if you want to run over all files
    try:
        return terminal(name, b)
    except Exception as e:
        print("Terminal {} failed: {}".format(name, e))
        return None

if __name__ == "__main__":
    if len(sys.argv) == 2:
        output_name = sys.argv[1]
        terminal(output_name, False)

    else:
        run_names = []
        for name in os.listdir('../../staged'):
            if "ob" in name:
                continue
            run_names.append(name[:-4])

        with Pool(16) as pool:
            pool.map(wrap_terminal, run_names)
