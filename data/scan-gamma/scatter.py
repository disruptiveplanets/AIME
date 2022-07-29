# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

LOW_ORDER_INDEX = 1
BASE_NAME = f"gamma"
LOW_GAMMA_0= 0
HIGH_GAMMA_0 = 2 * np.pi
NUM_DIVISIONS = 48 # Run with 12 cores per process.
LOW_ORDER = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)][LOW_ORDER_INDEX]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]
PERIGEE_GAMMA = -2.58971954411
MU = 3.986004418e14

def get_text(gamma_0):
    quartile = int((gamma_0 + np.pi / 4) / (np.pi / 2))
    high_lim = np.pi / 4 + quartile * np.pi / 2
    low_lim = - np.pi / 4 + quartile * np.pi / 2
    return """0, 3
120
5
1000
6000
0.00006464182, 0.00012928364, -0.00012928364
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
{}, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
{}, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5""".format(gamma_0, LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6], high_lim, low_lim)

for ratio_index, perigee in enumerate(np.linspace(LOW_GAMMA_0, HIGH_GAMMA_0, NUM_DIVISIONS)):
    f = open("../../staged/{}-{:02}.txt".format(BASE_NAME, ratio_index), 'w')
    f.write(get_text(perigee))
    f.close()
