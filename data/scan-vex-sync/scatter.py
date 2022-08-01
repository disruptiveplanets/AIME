# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

LOW_ORDER_INDEX = 1
BASE_NAME = f"vex-sync"
LOW_V_EX= 1000
HIGH_V_EX = 10000
NUM_DIVISIONS = 48 # Run with 12 cores per process.
LOW_ORDER = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)][LOW_ORDER_INDEX]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]
PERIGEE_GAMMA = -2.58971954411
MU = 3.986004418e14

def get_text(vex):
    perigee = 5 * 6_370_000
    semi_major_ax = MU / vex**2
    eccentricity = 1 + perigee / semi_major_ax
    semi_latus_rectum = -semi_major_ax * (1 - eccentricity**2)
    cos_true_anomaly = (semi_latus_rectum / (10 * perigee) - 1) / eccentricity
    eccentric_anomaly = np.arccosh((eccentricity + cos_true_anomaly) / (1 + eccentricity * cos_true_anomaly))
    time_to_perigee = np.sqrt(semi_major_ax**3 / MU) * (eccentricity * np.sinh(eccentric_anomaly) - eccentric_anomaly)
    print(vex, '\t', time_to_perigee / 3600)
    radians = 2 * np.pi / (9 * 3600) * time_to_perigee
    gamma_0 = (PERIGEE_GAMMA - radians) % (2 * np.pi)
    quartile = int((gamma_0 + np.pi / 4) / (np.pi / 2))
    high_lim = np.pi / 4 + quartile * np.pi / 2
    low_lim = - np.pi / 4 + quartile * np.pi / 2
    return """0, 3
120
5
1000
{}
0.00006464182, 0.00012928364, -0.00012928364
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
{}, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
{}, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5""".format(vex, gamma_0, LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6], high_lim, low_lim)

for ratio_index, perigee in enumerate(np.linspace(LOW_V_EX, HIGH_V_EX, NUM_DIVISIONS)):
    f = open("../../staged/{}-{:02}.txt".format(BASE_NAME, ratio_index), 'w')
    f.write(get_text(perigee))
    f.close()
