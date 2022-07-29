# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

BASE_NAME = "ob"
SPECIALS = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]
NUM_GENERATE = 12

def get_text(i, j20):
    return """2, 2
120
5
1000
4000
0.00024682682, 0, -0.00024682682
1.0, 0, 0, 0, 0, 0, 0, 0, {}
{}, {}, {}
0.78539816339, 0.12499, -0.0001
-0.78539816339, -0.12499, -0.24999
0.01, 1e-5""".format(j20, SPECIALS[i][0], SPECIALS[i][1], SPECIALS[i][2])

for i in range(len(SPECIALS)):
    if i == 0: continue
    for j2_index, oblateness in enumerate(10**np.linspace(-2, -1, NUM_GENERATE+1)[1:]):
        name ="../../staged/{}-{}-{:02}.txt".format(BASE_NAME, i, j2_index+24)
        f = open(name, 'w')
        f.write(get_text(i, -oblateness * 2))
        f.close()
