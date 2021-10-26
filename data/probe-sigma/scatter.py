# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

BASE_NAME = "sigma-1"
LOW_LOG_SIGMA = -5
HIGH_LOG_SIGMA = 0
NUM_SIGMA_DIVISIONS = 48 # Run with 12 cores per process.
SPECIALS = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]

def get_text(i, sigma):
    return """0, 3
3, 1
120
5
1000
4000
0.00012, 0.00022, 0.00032
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
0.78539816339, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
-0.78539816339, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
{}""".format(SPECIALS[i][0], SPECIALS[i][1], SPECIALS[i][2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6], sigma)

for i in range(len(SPECIALS)):
    for sigma_index, log_sigma in enumerate(np.linspace(LOW_LOG_SIGMA, HIGH_LOG_SIGMA, NUM_SIGMA_DIVISIONS)):
        f = open("../../staged/{}-{}-{:03}.txt".format(BASE_NAME, i, sigma_index), 'w')
        f.write(get_text(i, 10**log_sigma))
        f.close()
