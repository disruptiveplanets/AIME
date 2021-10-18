## Find coordinates that minimize the distance to all points in the triangle
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import minimize

TRIANGLE_BASE = 0.25
TRIANGLE_HEIGHT = 0.25
NUM_SUBDIVISIONS = 1000

tot_index = 0
points = []
for i in range(NUM_SUBDIVISIONS + 1):
    l2z = -TRIANGLE_HEIGHT / NUM_SUBDIVISIONS * i
    half_width = TRIANGLE_BASE / NUM_SUBDIVISIONS / 2 * i
    for j in range(1 + i):
        l2x = j * TRIANGLE_BASE / NUM_SUBDIVISIONS  - half_width
        points.append((l2x, l2z))

def dist(a, b):
    return np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)

def get_val(p):
    p1 = [0, p[0]]
    p2 = [p[1], p[2]]
    p3 = [-p[1], p[2]]
    val=0
    for p in points:
        d1 = dist(p1, p)
        d2 = dist(p2, p)
        d3 = dist(p3, p)
        val += min(d1, d2, d3)**2
    return val

out = minimize(get_val, (-5/12 * TRIANGLE_HEIGHT, 3/16 * TRIANGLE_HEIGHT, -19/24 * TRIANGLE_HEIGHT))
print(out.x)

plt.scatter(np.array(points)[:,0], np.array(points)[:,1], alpha=0.2)
plt.scatter([0, out.x[1], -out.x[1]], [out.x[0], out.x[2], out.x[2]])
plt.axvline(x=0)
plt.axhline(y=-19/24 * TRIANGLE_HEIGHT)
plt.axhline(y=-5/12 * TRIANGLE_HEIGHT)
plt.axvline(x=3 / 16 * TRIANGLE_HEIGHT)
plt.axvline(x=-3/16 * TRIANGLE_HEIGHT)
plt.show()
