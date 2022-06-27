# PROBLEM: this data doesn't have K3m.

import sys, os
import numpy as np
sys.path.append("../../density")
from core import Indicator, TrueShape, UncertaintyTracker
sys.path.append("../../density/mcmc")
from mcmc_core import MCMCAsteroid
from fe import FiniteElement
from contextlib import contextmanager

DIVISION = 99
FILE_PATH = "../../../data/"
NUM_TRIALS = 5
DEGREES_OF_FREEDOM = 5
MAX_REPEAT_COUNT = 100
N_AVERAGE_SAMPLES = 1000

NAMES = {
    "scan-perigee": (2,), 
    "probe-s-theta": (2,), 
    "probe-s-rho": (2,), 
    "scan-cadence": (2,),
    "scan-period": (2,), 
    "scan-am": (2,),
    "scan-vex": (2,), 
    "observation-gap": (2,),
    "scan-spin-pole": (2,),
    "cad-period-sync-contour": (2,2),
    "cad-speed-contour": (2,2)
}


@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

def get_indices(name, index_lengths):
    index_pos = 0
    indices = []
    for index, index_length in enumerate(index_lengths):
        if index_pos == 0:
            indices.append(int(name[index_pos-index_length:]))
        else:
            indices.append(int(name[index_pos-index_length:index_pos]))
        index_pos -= index_length + 1
    return np.array(indices, dtype=int)
    
def scan_directory(directory, index_lengths):
    if not os.path.exists(directory):
        raise Exception(f"Directory {directory} does not exist")

    max_indices = np.zeros_like(index_lengths)
    for dname in os.listdir(directory):
        if os.path.isdir(directory+'/'+dname):
            indices = get_indices(dname, index_lengths)
            max_indices = np.maximum(max_indices, indices)
    max_indices += 1
    uncs = np.ones(np.append(max_indices, [5])) * np.nan

    for dname in os.listdir(directory):
        run_name = directory+'/'+dname
        if not os.path.isdir(run_name):
            continue
        for fname in os.listdir(run_name):
            if not fname.endswith("-0-samples.npy"):
                continue
            indices = tuple(get_indices(dname, index_lengths))
            if not np.all(np.isnan(uncs[indices])):
                raise Exception(f"Index {indices} was already processed")
            unc = get_unc_for_file(run_name, run_name+'/'+fname)
            uncs[indices] = unc
    return uncs


def get_unc_for_file(dname, fname):
    # fname ends with -0-samples.npy
    with open(fname, 'rb') as f:
        array = np.load(f)
        if len(array) == 0: return np.nan
    
        # flat_samples = array.reshape(-1, array.shape[-1])

        # median = np.percentile(flat_samples, 50, axis=0)
        # up_sigma = np.percentile(flat_samples, 50 + 68.27 / 2, axis=0)
        # down_sigma = np.percentile(flat_samples, 50 - 68.27 / 2, axis=0)
        # original_sigma = ((up_sigma - median) - (down_sigma - median))/2

    with open(f"{fname[:-14]}.txt", 'r') as f:
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

    # Make asteroid
    am = radius

    division = DIVISION * radius / 1000

    if am < 2 * division:
        return np.nan
    k20 = theta_true[2]
    k22 = theta_true[1]

    a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22)
    b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22)
    c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20)

    max_radius = int(max(a, b, c) + 4 * division)
    
    short_name = fname[:-14]
    short_name = short_name[short_name.rfind('/')+1:]
    print(short_name)

    # Do not regenerate if the file was already done.
    unc_tracker = UncertaintyTracker()
    for run_index in range(NUM_TRIALS):
        repeat_num = 0
        repeat = True
        while repeat:
            repeat = False
            generate = not os.path.exists(f"ast-{short_name}-{run_index}-fe.npy")
            print(f"Trial {run_index}. Generating: {generate}")
            with suppress_stdout():
                asteroid = MCMCAsteroid(f"ast-{short_name}-{run_index}", fname, Indicator.ell(radius, k22, k20), TrueShape.uniform(),
                    am, division, max_radius, DEGREES_OF_FREEDOM, am)
                this_unc_tracker = asteroid.pipeline(FiniteElement, False, generate=generate, n_samples=N_AVERAGE_SAMPLES)
            if this_unc_tracker is None:
                print("Failed. Had to repeat")
                repeat = True
                repeat_num += 1
            if repeat_num == MAX_REPEAT_COUNT:
                # Give up. Keep repeat = True
                break # Don't add anything to the uncs list
        if not repeat:
            # Success
            unc_tracker += this_unc_tracker

    density_map, uncertainty_map = unc_tracker.generate()
    if density_map is None:
        return np.nan
    uncertainty_ratio = uncertainty_map / density_map
    return np.array([
        np.nanpercentile(uncertainty_ratio, 100 - (100 - 95.45) / 2),
        np.nanpercentile(uncertainty_ratio, 100 - (100 - 68.27) / 2),
        np.nanpercentile(uncertainty_ratio, 50),
        np.nanpercentile(uncertainty_ratio, (100 - 68.27) / 2),
        np.nanpercentile(uncertainty_ratio, (100 - 95.45) / 2),
    ])

def scan_specific(directory, threshold):
    index_lengths = NAMES[directory]
    uncs = scan_directory(FILE_PATH + directory, index_lengths)
    print(uncs)
    with open(directory + ".npy", 'wb') as f:
        np.save(f, uncs)
    if len(uncs.shape) == 1:
        where_more_than_one = np.where(uncs>threshold)[0]
        if len(where_more_than_one) > 0:
            threshold_index_left = where_more_than_one[0]
            threshold_index_right = where_more_than_one[-1]

        else:
            threshold_index_left = None
            threshold_index_right = None
        print("Increasing threshold", threshold_index_left, "\t", "Decreasing threshold", threshold_index_right)
        
def scan_all(thresholds):
    print("Thresholds", thresholds)
    for name, index_lengths in NAMES.items():
        uncs = scan_directory(FILE_PATH + name, index_lengths)
        print(uncs)
        with open(name+".npy", 'wb') as f:
            np.save(f, uncs)

        if len(uncs.shape) == 1:
            print(name, end='')
            for threshold in thresholds:
                where_more_than_one = np.where(uncs>threshold)[0]
                if len(where_more_than_one) > 0:
                    threshold_index_left = where_more_than_one[0]
                    threshold_index_right = where_more_than_one[-1]

                else:
                    threshold_index_left = None
                    threshold_index_right = None
                    
                print(f"\tinc.\t{threshold_index_left}\tdec.\t{threshold_index_right}")

if __name__ == "__main__":
    scan_specific(sys.argv[1], 1)
    #scan_specific(sys.argv[1], 0.2)
    #scan_all((1,))
