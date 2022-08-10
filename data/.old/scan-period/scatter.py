# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

LOW_ORDER_INDEX = 1
BASE_NAME = f"per"
LOW_LOG_PERIOD = np.log10(2.4)
HIGH_LOG_PERIOD = np.log10(24)
NUM_DIVISIONS = 48 # Run with 8 cores per process.
LOW_ORDER = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)][LOW_ORDER_INDEX]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]

# Fix product, vary ratio

def get_text(period):
    omegax, omegay, omegaz = np.array([0.00006464182, 0.00012928364, -0.00012928364]) * (9 / period)
    return """0, 3
120
5
1000
6000
{}, {}, {}
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
0.78539816339, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
-0.78539816339, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5""".format(omegax, omegay, omegaz, LOW_ORDER[0], LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6])

for index, period in enumerate(10**np.linspace(LOW_LOG_PERIOD, HIGH_LOG_PERIOD, NUM_DIVISIONS)):
    f = open("../../staged/{}-{:02}.txt".format(BASE_NAME, index), 'w')
    f.write(get_text(period))
    f.close()
