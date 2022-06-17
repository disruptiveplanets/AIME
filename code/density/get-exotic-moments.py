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
density_factor_low = 0.5
density_factor_high = 2
core_shift_low = core_displacement * (core_vol * density_factor_low) / ellipsoid_vol
core_shift_high = core_displacement * (core_vol * density_factor_high) / ellipsoid_vol
print("Core mass fraction:", (core_vol * density_factor) / (core_vol * density_factor+ellipsoid_vol))
print(f"Core shift low: {core_shift_low}\tCore shift high: {core_shift_high}")


asteroids = [
    #("sph", Indicator.sph(ELLIPSOID_AM), lambda x,y,z: 1),
    #("ells", Indicator.ell(ELLIPSOID_AM, k22s, k20s), lambda x,y,z: 1),
    #("ella", Indicator.ell(ELLIPSOID_AM, k22a, k20a), lambda x,y,z: 1, False),
    #("tet", Indicator.tet(ELLIPSOID_AM), lambda x,y,z: 1),
    #("db", Indicator.dumbbell(ELLIPSOID_AM), lambda x,y,z: 1),
    #("in", Indicator.ell(ELLIPSOID_AM, k22a, k20a), TrueShape.in_(ELLIPSOID_AM)),
    #("out", Indicator.ell(ELLIPSOID_AM, k22a, k20a), TrueShape.out(ELLIPSOID_AM)),
    #("in-sph", Indicator.sph(ELLIPSOID_AM), TrueShape.in_sph(ELLIPSOID_AM)),
    #("out-sph", Indicator.sph(ELLIPSOID_AM), TrueShape.out_sph(ELLIPSOID_AM)),
    #("blob", Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -lump_shift), TrueShape.blob(ELLIPSOID_AM, k22a, k20a), False),
    #("rot-blob", Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -lump_shift), TrueShape.rot_blob(ELLIPSOID_AM, k22a, k20a)),
    #("core-ell", Indicator.ell(ELLIPSOID_AM, k22a, k20a), TrueShape.core(ELLIPSOID_AM, k22a, k20a, 3, 0.65), False),
    #("core-sph-3", Indicator.ell(ELLIPSOID_AM, k22a, k20a), TrueShape.core_sph(3, 500), False),
    #("core-sph-1.5", Indicator.ell(ELLIPSOID_AM, k22a, k20a), TrueShape.core_sph(1.5, 500), False),
    ("core-move-3", Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_high), TrueShape.core_shift(3, 500, core_displacement), True),
    #("core-move-1.5", Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_low), TrueShape.core_shift(1.5, 500, core_displacement), True),
]

def get_klms(index):
    name, indicator, true_shape, recalc_am = asteroids[index]
    if recalc_am: 
        asteroid = Asteroid(name, ELLIPSOID_AM, division, max_radius, indicator, None)
        rlms = asteroid.moment_field(surface_am=ELLIPSOID_AM)
        surface_am_sqr = np.sum(rlms[-1] * asteroid.indicator_map) / np.sum(asteroid.indicator_map)
        surface_am = np.sqrt(surface_am_sqr)
        print("New surface_am", surface_am)
    else: 
        surface_am = ELLIPSOID_AM

    asteroid = Asteroid(name, surface_am, division, max_radius, indicator, None)
    density = asteroid.map_np(true_shape)

    # plt.imshow(density[len(density)//2, :, :])
    # plt.figure()
    # plt.imshow(density[:, len(density)//2, :])
    # plt.figure()
    # plt.imshow(density[:, :, len(density)//2])
    # plt.show()

    rlms = asteroid.moment_field(surface_am=surface_am)

    i = 0
    klms = []
    for l in range(0, 4):
        for m in range(-l, l+1):
            klms.append(np.sum(rlms[i] * density) * division**3)
            i += 1
    bulk_am = np.sqrt(np.sum(rlms[i] * density) * division**3 / klms[0])
    klms = np.array(klms)
    klms[1:] /= bulk_am**2 # Add in nonlinear /a^2 term, but keep the first term as the mass.
    klms /= klms[0] # Normalize by mass

    # Convert to old am by getting rid of surface_am
    i = 1
    for l in range(1, 4):
        for m in range(-l, l+1):
            klms[i] *= (surface_am / bulk_am)**(l - 2)
            i += 1

    return name, np.append(klms, bulk_am)

with Pool() as pool:
    results = pool.map(get_klms, range(len(asteroids)))

for name, klms in results:
    print(name)
    for k in klms:
        print(k)
    print()
    
