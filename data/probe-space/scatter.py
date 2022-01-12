# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

BASE_NAME = "param-"
SIGMA = 0.01
NUM_SUBDIVISIONS = 16 # run with 24?
TRIANGLE_BASE = 0.24448
TRIANGLE_HEIGHT = 0.24448
INITIAL_ANGLE = 0.39269908169
SPECIALS = [(0, -0.09766608), (0.05200629 , -0.2021978), (-0.05200629 , -0.2021978)]

def get_text(params):
    return """0, 2
120
5
1000
6000
0.00006464182, 0.00012928364, -0.00012928364
1.0
{}, {}, {}
0.78539816339, 0.12499, -0.0001
-0.78539816339, -0.12499, -0.24999
{}, 0""".format(INITIAL_ANGLE, params[0], params[1], SIGMA)

tot_index = 0
l2zs = []
l2xs = []
for i in range(NUM_SUBDIVISIONS+1):
    l2z = -TRIANGLE_HEIGHT / NUM_SUBDIVISIONS * i
    half_width = TRIANGLE_BASE / NUM_SUBDIVISIONS / 2 * i
    for j in range(0, i+1):
        l2x = j * TRIANGLE_BASE / NUM_SUBDIVISIONS  - half_width
        print("{:03} Generating theta2={}, theta3={}".format(tot_index, l2x, l2z))

        f = open("../../staged/{}{:03}.txt".format(BASE_NAME, tot_index), 'w')
        f.write(get_text((l2x, l2z)))
        f.close()

        tot_index += 1
        l2xs.append(l2x)
        l2zs.append(l2z)

plt.scatter(l2xs, l2zs)
plt.scatter(np.array(SPECIALS)[:,0], np.array(SPECIALS)[:,1])
plt.xlim(-0.125, 0.125)
plt.ylim(-0.25, 0)
plt.show()
