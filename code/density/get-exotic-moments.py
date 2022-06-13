import numpy as np
from core import Asteroid, Indicator, TrueShape
from multiprocessing import Pool
import matplotlib.pyplot as plt

division = 9
max_radius = 2000
ELLIPSOID_AM = 1000
k22a, k20a = -0.05200629, -0.2021978
k22s, k20s = 0, -0.09766608
b = np.sqrt(5/3) * ELLIPSOID_AM * np.sqrt(1 - 2 * k20a - 12 * k22a)
a = np.sqrt(5/3) * ELLIPSOID_AM * np.sqrt(1 - 2 * k20a + 12 * k22a)
c = np.sqrt(5/3) * ELLIPSOID_AM * np.sqrt(1 + 4 * k20a)



blob_displacement = 500
blob_rad = 300
blob_vol = np.pi * 4 / 3 * blob_rad**3
ellipsoid_vol = np.pi * 4 / 3 * a * b * c
density_factor = 5
lump_shift = blob_displacement * (blob_vol * density_factor) / ellipsoid_vol
print("Blob Mass fraction:", (blob_vol * density_factor) / (blob_vol * density_factor+ellipsoid_vol))
print("Blob Lump shift:", lump_shift)

core_displacement = 300
core_rad = 500
core_vol = np.pi * 4 / 3 * core_rad**3
ellipsoid_vol = np.pi * 4 / 3 * a * b * c
density_factor = 1.5
core_shift = core_displacement * (core_vol * density_factor) / ellipsoid_vol
print("Core mass fraction:", (core_vol * density_factor) / (core_vol * density_factor+ellipsoid_vol))
print("Core shift:", core_shift)


asteroids = [
    #("sph", Indicator.sph(ELLIPSOID_AM), lambda x,y,z: 1),
    #("ells", Indicator.ell(ELLIPSOID_AM, k22s, k20s), lambda x,y,z: 1),
    #("ella", Indicator.ell(ELLIPSOID_AM, k22a, k20a), lambda x,y,z: 1),
    #("tet", Indicator.tet(ELLIPSOID_AM), lambda x,y,z: 1),
    #("db", Indicator.dumbbell(ELLIPSOID_AM), lambda x,y,z: 1),
    #("in", Indicator.ell(ELLIPSOID_AM, k22a, k20a), TrueShape.in_(ELLIPSOID_AM)),
    #("out", Indicator.ell(ELLIPSOID_AM, k22a, k20a), TrueShape.out(ELLIPSOID_AM)),
    #("in-sph", Indicator.sph(ELLIPSOID_AM), TrueShape.in_sph(ELLIPSOID_AM)),
    #("out-sph", Indicator.sph(ELLIPSOID_AM), TrueShape.out_sph(ELLIPSOID_AM)),
    #("blob", Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -lump_shift), TrueShape.blob(ELLIPSOID_AM, k22a, k20a)),
    #("rot-blob", Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -lump_shift), TrueShape.rot_blob(ELLIPSOID_AM, k22a, k20a)),
    #("core-sph", Indicator.ell(ELLIPSOID_AM, k22a, k20a), TrueShape.core_sph(ELLIPSOID_AM, 1.5, 500), False),
    ("core-move", Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift), TrueShape.core_shift(1.5, 500, core_displacement), True),
]

def get_klms(index):
    name, indicator, generator, recalc_am = asteroids[index]
    if recalc_am: 
        asteroid = Asteroid(name, "", ELLIPSOID_AM, division, max_radius, indicator, None, used_bulk_am=ELLIPSOID_AM)
        rlms = asteroid.moment_field(max_l=3, surface_am=ELLIPSOID_AM)
        surface_am_sqr = np.sum(rlms[-1] * asteroid.indicator_map) / np.sum(asteroid.indicator_map)
        am = np.sqrt(surface_am_sqr)
        print("New am", am)
    else: 
        am = ELLIPSOID_AM
    asteroid = Asteroid(name, "", am, division, max_radius, indicator, None, used_bulk_am=am)
    density = asteroid.map_np(generator)

    # plt.imshow(density[len(density)//2, :, :])
    # plt.figure()
    # plt.imshow(density[:, len(density)//2, :])
    # plt.figure()
    # plt.imshow(density[:, :, len(density)//2])
    # plt.show()

    rlms = asteroid.moment_field(max_l=3, surface_am=am)

    i = 0
    klms = []
    for l in range(0, 4):
        for m in range(-l, l+1):
            klms.append(np.sum(rlms[i] * density) * division**3)
            i += 1
    am_this = np.sqrt(np.sum(rlms[i] * density) * division**3 / klms[0])
    klms = np.array(klms) / klms[0]

    # Rescale by am
    i = 0
    for l in range(0, 4):
        for m in range(-l, l+1):
            klms[i] /= am_this**l
            i += 1

    return name, np.append(klms, am_this)

with Pool() as pool:
    results = pool.map(get_klms, range(len(asteroids)))

for name, klms in results:
    print(name)
    for k in klms:
        print(k)
    print()
    
