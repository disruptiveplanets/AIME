# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

LOW_ORDER_INDEX = 1
BASE_NAME = f"eq-spin-wide"
LOW_V_EX= 1000
HIGH_V_EX = 10000
NUM_DIVISIONS = 36 # Run with 12 cores per process.
LOW_ORDER = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)][LOW_ORDER_INDEX]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]

SPIN_MAG = 0.00019392547

def get_text(phi, gamma0):
    return """0, 3
120
5
1000
10000
{}, {}, 0
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
0.78539816339, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
-0.78539816339, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5""".format(SPIN_MAG * np.cos(phi), SPIN_MAG * np.sin(phi), gamma0, LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6])

for ratio_index, phi in enumerate(np.linspace(0, 2 * np.pi, NUM_DIVISIONS)):
    f = open("../../staged/{}-{:02}.txt".format(BASE_NAME, ratio_index), 'w')
    f.write(get_text(phi, 0.39269908169))
    f.close()