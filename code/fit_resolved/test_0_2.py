import numpy as np
from math import factorial

DATA_WIDTH = 500
EARTH_RADIUS = 6_370_000
x = np.linspace(0, 12 * np.pi, DATA_WIDTH)

def simulate(cadence, jlms, theta, radius, spinx, spiny, spinz, initial_roll, impact_parameter,
    speed, mu, central_radius, cut=None):

    consider_x = x[:int(DATA_WIDTH / 2 * (1 + cut / 10)) - 1] if cut >= 0 else x

    params = [initial_roll]
    for j in theta:
        params.append(j)

    y1 = np.zeros_like(consider_x)
    y2 = np.zeros_like(consider_x)
    y3 = np.zeros_like(consider_x)
    for n, p in enumerate(params):
        mul = 1
        if n >= 3:
            mul = 0.000001
        y1 += mul * p**2 * np.sin(consider_x * n)
        y2 += mul * p**2 * np.cos(consider_x * n)
        y3 += mul * p**2 * np.cos(consider_x * (n+1))

    y = []
    for i in range(consider_x.shape[0]):
        y.append(y1[i])
        y.append(y2[i])
        y.append(y3[i])

    return y
