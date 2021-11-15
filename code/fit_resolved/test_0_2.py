import numpy as np
from math import factorial

DATA_WIDTH = 500
EARTH_RADIUS = 6_370_000
x = np.linspace(0, 12 * np.pi, DATA_WIDTH)

def simulate(cadence, jlms, theta, radius, spinx, spiny, spinz, initial_roll, impact_parameter,
    speed, mu, central_radius, cut):

    params = [initial_roll]
    for j in theta:
        params.append(j)

    y1 = np.zeros_like(x)
    y2 = np.zeros_like(x)
    y3 = np.zeros_like(x)
    for n, p in enumerate(params):
        mul = 1
        if n >= 3:
            mul = 0.000001
        y1 += mul * p**2 * np.sin(x * n)
        y2 += mul * p**2 * np.cos(x * n)
        y3 += mul * p**2 * np.cos(x * (n+1))

    y = []
    for i in range(DATA_WIDTH):
        y.append(y1[i])
        y.append(y2[i])
        y.append(y3[i])

    return y
