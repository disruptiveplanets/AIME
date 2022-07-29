# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

LOW_ORDER_INDEX = 1
BASE_NAME = f"cont-ps"
NUM_DIVISIONS = 48 # Run with 12 cores per process.
LOW_ORDER = [(None, 0, -0.09766608), (None, 0.05200629, -0.2021978)][LOW_ORDER_INDEX]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]
PERIGEE_GAMMA = -2.58971954411

def get_text(period, cadence):
    gamma_0 = (PERIGEE_GAMMA - 13.728 / period * 2 * np.pi) % (2 * np.pi)
    quartile = int((gamma_0 + np.pi / 4) / (np.pi / 2))
    high_lim = np.pi / 4 + quartile * np.pi / 2
    low_lim = - np.pi / 4 + quartile * np.pi / 2
    return """0, 3
{}
5
1000
6000
{}, {}, {}
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
{}, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
{}, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5""".format(cadence, 0.00006464182 * 9 / period, 0.00012928364 * 9 / period, -0.00012928364 * 9 / period,
    gamma_0, LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6], high_lim, low_lim)

cadences = 60 * np.linspace(2, 60, 20)
periods = 10**np.linspace(np.log10(4.5), np.log10(18), 10)

for icad, cad in enumerate(cadences):
    for iperiod, period in enumerate(periods):
        with open("../../staged/{}-{:02}-{:02}.txt".format(BASE_NAME, icad, iperiod), 'w') as f:
            f.write(get_text(period, cad))

