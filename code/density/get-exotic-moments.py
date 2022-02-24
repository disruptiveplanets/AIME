import numpy as np
import sys
sys.path.append("../../code/density")
from core import Asteroid, Indicator
from multiprocessing import Pool

division = 9
max_radius = 1500
am = 1000
k22a, k20a = 0.05200629, -0.2021978
k22s, k20s = 0, -0.09766608
b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20s - 12 * k22s)
a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20s + 12 * k22s)
c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20s)


asteroids = [
    ("sph", Indicator.sph(am), lambda x,y,z: 1),
    ("ells", Indicator.ell(am, k22s, k20s), lambda x,y,z: 1),
    ("ella", Indicator.ell(am, k22a, k20a), lambda x,y,z: 1),
    ("tet", Indicator.tet(am), lambda x,y,z: 1),
    ("db", Indicator.dumbbell(am), lambda x,y,z: 1),
    ("in", Indicator.ell(am, k22a, k20a), lambda x,y,z: np.exp(-0.5 * x*x/(a*a) + y*y/(b*b) + z*z/(c*c))),
    ("out", Indicator.ell(am, k22a, k20a), lambda x,y,z: np.exp(0.5 * x*x/(a*a) + y*y/(b*b) + z*z/(c*c))),
]

def get_klms(index):
    name, indicator, generator = asteroids[index]
    asteroid = Asteroid("", am, division, max_radius, indicator)
    rlms = asteroid.moment_field(max_l=3)

    density = asteroid.map_np(generator)

    i = 0
    klms = []
    for l in range(0, 4):
        for m in range(-l, l+1):
            klms.append(np.sum(rlms[i] * density) / am**l * division**3)
            i += 1
    klms.append(np.sqrt(np.sum(rlms[i]) / division**3))

    return name, klms / klms[0]

with Pool() as pool:
    results = pool.map(get_klms, range(len(asteroids)))

for name, klms in results:
    print(name)
    for k in klms:
        print(k)
    print()
    