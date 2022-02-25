import numpy as np
from core import Asteroid, Indicator
from multiprocessing import Pool

division = 39
max_radius = 2000
am = 1000
k22a, k20a = -0.05200629, -0.2021978
k22s, k20s = 0, -0.09766608
b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20a - 12 * k22a)
a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20a + 12 * k22a)
c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20a)


blob_length = 500
blob_rad = 500
blob_vol = np.pi * 4 / 3 * blob_rad**3
ellipsoid_vol = np.pi * 4 / 3 * a * b * c
density_factor = 4
lump_shift = blob_rad * density_factor * blob_vol / (ellipsoid_vol + density_factor * blob_vol)
print(lump_shift)
def lump(x, y, z):
    o = np.ones_like(x)
    o[(x-blob_length+lump_shift)**2 + y**2 + z**2 < blob_rad**2] = density_factor + 1
    return o


asteroids = [
    ("sph", Indicator.sph(am), lambda x,y,z: 1),
    ("ells", Indicator.ell(am, k22s, k20s), lambda x,y,z: 1),
    ("ella", Indicator.ell(am, k22a, k20a), lambda x,y,z: 1),
    ("tet", Indicator.tet(am), lambda x,y,z: 1),
    ("db", Indicator.dumbbell(am), lambda x,y,z: 1),
    ("in", Indicator.ell(am, k22a, k20a), lambda x,y,z: np.exp(-0.5 * x*x/(a*a) + y*y/(b*b) + z*z/(c*c))),
    ("out", Indicator.ell(am, k22a, k20a), lambda x,y,z: np.exp(0.5 * x*x/(a*a) + y*y/(b*b) + z*z/(c*c))),
    ("blob", Indicator.ell_x_shift(am, k22a, k20a, -lump_shift), lump),
]

def get_klms(index):
    name, indicator, generator = asteroids[index]
    asteroid = Asteroid("", am, division, max_radius, indicator)
    density = asteroid.map_np(generator)

    rlms = asteroid.moment_field(max_l=3)

    i = 0
    klms = []
    for l in range(0, 4):
        for m in range(-l, l+1):
            klms.append(np.sum(rlms[i] * density) / am**l * division**3)
            i += 1
    am_this = np.sqrt(np.sum(rlms[i] * density) * division**3 / klms[0])
    klms = np.array(klms) / klms[0]

    return name, np.append(klms, am_this)

with Pool(4) as pool:
    results = pool.map(get_klms, range(len(asteroids)))

for name, klms in results:
    print(name)
    for k in klms:
        print(k)
    print()
    