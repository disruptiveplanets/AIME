# PROBLEM: this data doesn't have K3m.

import sys, os
import numpy as np
sys.path.append("../density")
from likelihood import Likelihood as Model
#from harmonic import Harmonic as Model
#from finite_element import FiniteElement as Model
from core import Asteroid, Indicator

DIVISION = 29
FILE_PATH = "../../data/"

NAMES = [
    #"scan-perigee", 
    #"probe-s-theta", 
    #"probe-s-rho", 
    #"scan-cadence",
    #"scan-period", 
    #"scan-am", 
    #"scan-vex", 
    #"observation-gap"
    "scan-spin-pole",
    ]

from contextlib import contextmanager
import sys, os

@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:  
            yield
        finally:
            sys.stdout = old_stdout

def scan_directory(directory):
    if not os.path.exists(directory):
        raise Exception(f"Directory {directory} does not exist")

    max_index = 0
    for dname in os.listdir(directory):
        if os.path.isdir(directory+'/'+dname):
            max_index = max(max_index, int(dname[-2:]))
            
    uncs = np.ones(max_index + 1) * np.nan

    for dname in os.listdir(directory):
        if not os.path.isdir(directory+'/'+dname):
            continue
        for fname in os.listdir(directory+'/'+dname):
            if not fname.endswith("-0-samples.npy"):
                continue
            unc = get_unc_for_file(directory+'/'+dname, directory+'/'+dname+'/'+fname)
            index = int(dname[-2:])
            uncs[index] = unc
    return uncs


def get_unc_for_file(dname, fname):
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
    if am < 2 * DIVISION:
        return np.nan
    k20 = theta_true[2]
    k22 = theta_true[1]

    a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22)
    b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22)
    c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20)

    max_radius = int(max(a, b, c) + 4 * DIVISION)
    
    asteroid = Asteroid(f"ast-{fname}", fname, am, DIVISION, max_radius, Indicator.ell(radius, k22, k20), None)

    with suppress_stdout():
        method = Model(asteroid)
    method.solve()
    uncs = method.map_unc()

    density_uncertainty = np.nanmedian(uncs)
    return density_uncertainty

def scan_specific(directory, threshold):
    from matplotlib import pyplot as plt
    uncs = scan_directory(FILE_PATH + sys.argv[1])
    print(uncs)
    with open(sys.argv[1]+".npy", 'wb') as f:
        np.save(f, uncs)
    where_more_than_one = np.where(uncs>threshold)[0]
    if len(where_more_than_one) > 0:
        threshold_index_left = where_more_than_one[0]
        threshold_index_right = where_more_than_one[-1]

    else:
        threshold_index_left = None
        threshold_index_right = None
    print("Increasing threshold", threshold_index_left, "\t", "Decreasing threshold", threshold_index_right)
    
    plt.plot(uncs)
    plt.xlabel("index")
    plt.ylabel("Distribution uncertainty")
    plt.axhline(y=threshold)
    plt.show()

def scan_all(thresholds):
    print("Thresholds", thresholds)
    for name in NAMES:
        uncs = scan_directory(FILE_PATH + name)
        print(uncs)
        with open(name+".npy", 'wb') as f:
            np.save(f, uncs)
        
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
    #scan_specific(sys.argv[1], 1)
    scan_all((1, 0.2))
