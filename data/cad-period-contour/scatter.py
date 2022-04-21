# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

LOW_ORDER_INDEX = 1
BASE_NAME = f"cont-p"
NUM_DIVISIONS = 48 # Run with 12 cores per process.
LOW_ORDER = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)][LOW_ORDER_INDEX]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]

def get_text(period, cadence):
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
{}""".format(cadence, 0.00006464182 * 9 / period, 0.00012928364 * 9 / period, -0.00012928364 * 9 / period,
    LOW_ORDER[0], LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6])

cadences = 60 * np.linspace(2, 60, 20)
periods = 10**np.linspace(np.log10(4.5), np.log10(18), 10)

for icad, cad in enumerate(cadences):
    for iperiod, period in enumerate(periods):
        with open("../../staged/{}-{:02}-{:02}.txt".format(BASE_NAME, icad, iperiod), 'w') as f:
            f.write(get_text(period, cad))

