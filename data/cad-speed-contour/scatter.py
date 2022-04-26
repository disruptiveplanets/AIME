# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

LOW_ORDER_INDEX = 1
BASE_NAME = f"cont-s"
NUM_DIVISIONS = 48 # Run with 12 cores per process.
LOW_ORDER = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)][LOW_ORDER_INDEX]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]
PERIGEE_GAMMA = -2.58971954411

def get_text(ratio, cadence):
    gamma_0 = (PERIGEE_GAMMA - 13.728 / 9 / ratio * 2 * np.pi) % (2 * np.pi)
    quartile = int((gamma_0 + np.pi / 4) / (np.pi / 2))
    high_lim = np.pi / 4 + quartile * np.pi / 2
    low_lim = - np.pi / 4 + quartile * np.pi / 2
    return """0, 3
{}
5
1000
6000
0.00006464182, 0.00012928364, -0.00012928364
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
{}, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
{}, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5
{}""".format(cadence, gamma_0, LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6], high_lim, low_lim, ratio)

cadences = 60 * np.linspace(2, 60, 20)
factors = 10**np.linspace(np.log10(0.5), np.log10(2), 10)

for icad, cad in enumerate(cadences):
    for ifactor, factor in enumerate(factors):
        with open("../../staged/{}-{:02}-{:02}.txt".format(BASE_NAME, icad, ifactor), 'w') as f:
            f.write(get_text(factor, cad))

