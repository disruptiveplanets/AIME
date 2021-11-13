# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

BASE_NAME = "ff"
SPECIALS = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]
SIGMA = 0.01

def get_text(i, flyby):
    return """3, 3
3, 1
120
{}
1000
4000
0.00012, 0.00022, 0.00032
1.0, -0.29703, 0, 0, 4.45545, 0, 0, 0, -8.91089, -44.5545, 0, 0, 0, 133.663, 0, 0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
0.78539816339, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
-0.78539816339, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
0.001""".format(flyby * 60, SPECIALS[i][0], SPECIALS[i][1], SPECIALS[i][2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6])

for i in range(len(SPECIALS)):
    for flyby_index, flyby in enumerate(np.arange(1, 3, 0.2)):
        name ="../../staged/{}-{}-{:02}.txt".format(BASE_NAME, i, flyby_index)
        f = open(name, 'w')
        f.write(get_text(i, flyby))
        f.close()
