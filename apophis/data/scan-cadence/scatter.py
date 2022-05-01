# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

BASE_NAME = f"cad"
LOW_CADENCE= 2 * 60
HIGH_CADENCE = 60 * 60
NUM_DIVISIONS = 48 # Run with 8 cores per process.
LOW_ORDER = (0.39269908169, 0.0316092, -0.0747126)
PRECESS = 1.39269908169
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]

# Fix product, vary ratio

def get_text(cadence):
    return """0, 3
{}
263
{}
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
1.57079632679, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
0.0000001, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5""".format(cadence, PRECESS, LOW_ORDER[0], LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6])

for ratio_index, cadence in enumerate(np.linspace(LOW_CADENCE, HIGH_CADENCE, NUM_DIVISIONS)):
    f = open("../../staged/{}-{:02}.txt".format(BASE_NAME, ratio_index), 'w')
    f.write(get_text(cadence))
    f.close()
