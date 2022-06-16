# PROBLEM: this data doesn't have K3m.

import sys, os
import numpy as np
from code.density.core import TrueShape
sys.path.append("../../../density")
from core import Indicator
sys.path.append("../../../density/mcmc")
from mcmc_core import MCMCAsteroid
from fe import FiniteElement

DOFS = [9, 7, 5, 3, 2]
NUM_TRIALS = 5
index = int(sys.argv[2])
assert(0 <= index < len(DOFS) * NUM_TRIALS)
NUM_DOF = DOFS[index % NUM_TRIALS]
TRIAL_INDEX = index // NUM_TRIALS
DIVISION = 99
FILE_PATH = "../../../../data/"
MAX_REPEAT_COUNT = 100

print(f"{NUM_DOF} degrees of freedom")
print(f"Trial {TRIAL_INDEX}")
print(f"Division {DIVISION}")

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

    uncs = np.ones(max_indices) * np.nan
    deviations = np.ones(max_indices) * np.nan

    for dname in os.listdir(directory):
        run_name = directory+'/'+dname
        if not os.path.isdir(run_name):
            continue
        for fname in os.listdir(run_name):
            if not fname.endswith("-0-samples.npy"):
                continue
            indices = tuple(get_indices(dname, index_lengths))
            if not np.isnan(uncs[indices]):
                raise Exception(f"Index {indices} was already processed")
            deviation, unc = get_unc_for_file(run_name, run_name+'/'+fname)
            uncs[indices] = unc
            deviations[indices] = deviation

    return deviations, uncs


def get_unc_for_file(dname, fname):
    # fname ends with -0-samples.npy
    with open(fname, 'rb') as f:
        array = np.load(f)
        if len(array) == 0: return np.nan, np.nan, np.nan
    
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
    generate = not os.path.exists(f"cast-{NUM_DOF}-{TRIAL_INDEX}-{short_name}-fe.npy")
    repeat_num = 0
    repeat = True
    while repeat:
        repeat = False
        with suppress_stdout():
            asteroid = MCMCAsteroid(f"cast-{NUM_DOF}-{TRIAL_INDEX}-{short_name}", fname, Indicator.ell(radius, k22, k20), TrueShape.usniform(),
                am, division, max_radius, NUM_DOF, am, generate=generate)
            deviation, uncs = asteroid.pipeline(FiniteElement, False, generate=generate)

        if np.any(np.isnan(uncs)):
            print("Failed. Had to repeat")
            repeat = True
            repeat_num += 1
        if repeat_num == MAX_REPEAT_COUNT:
            # Give up
            return np.nan, np.nan

    return deviation, uncs

def scan_specific(directory):
    index_lengths = NAMES[directory]
    devs, uncs = scan_directory(FILE_PATH + directory, index_lengths)
    print(uncs)
    with open(f"c{directory}-{NUM_DOF}-{TRIAL_INDEX}.npy", 'wb') as f:
        np.save(f, (devs, uncs))

if __name__ == "__main__":
    scan_specific(sys.argv[1])
