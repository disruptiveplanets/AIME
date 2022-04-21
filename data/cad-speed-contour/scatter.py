# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

LOW_ORDER_INDEX = 1
BASE_NAME = f"cont-s"
NUM_DIVISIONS = 48 # Run with 12 cores per process.
LOW_ORDER = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)][LOW_ORDER_INDEX]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]

def get_text(ratio, cadence):
    return """0, 3
{}
5
1000
6
0.00006464182, 0.00012928364, -0.00012928364
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
0.78539816339, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
-0.78539816339, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5
{}""".format(cadence, LOW_ORDER[0], LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6], ratio)

cadences = 60 * np.linspace(2, 60, 20)
factors = 10**np.linspace(np.log10(0.5), np.log10(2), 10)

for icad, cad in enumerate(cadences):
    for ifactor, factor in enumerate(factors):
        with open("../../staged/{}-{:02}-{:02}.txt".format(BASE_NAME, icad, ifactor), 'w') as f:
            f.write(get_text(factor, cad))

