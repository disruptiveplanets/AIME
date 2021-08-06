import asteroids, time, sys
import matplotlib.pyplot as plt
from random_vector import *
import numpy as np

CADENCE = 120

EARTH_RADIUS = 6370000
EARTH_MASS =5.972e24

spin = [0.00012, 0.00022, 0.00032]
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
jlms = [EARTH_MASS, 0, 6e22, 7e22]
klmss = [
    [1.2e6, 1.1e5, -4.9e5],
    [1.2e6, 1.1e5, -4.9e5],
]
SIGMA = 0.02 * np.sqrt(spin[0]**2 + spin[1]**2 + spin[2]**2)
initial_rolls = [
    0,
    0,
]

SHOW_ERROR = False


plt.figure(figsize=(12, 4))

for data_iter in range(len(initial_rolls)):
    initial_roll = initial_rolls[data_iter]
    klms = klmss[data_iter]

    start = time.time()
    try:
        resolved_data = asteroids.simulate(CADENCE, jlms, klms,
            spin[0], spin[1], spin[2], initial_roll, impact_parameter, speed)
    except RuntimeError as err:
        print(err)
        sys.exit()
    end = time.time()
    print("Time taken: {} s".format(end - start))

    f = open("resolved.dat", 'w')
    for i, dat in enumerate(resolved_data):
        f.write(str(dat) + ' ')
        if i %3 == 2:
            f.write('\n')
    f.close()

    if SHOW_ERROR:
        y, yerr = randomize_flat(resolved_data, SIGMA)

        x_display = np.arange(len(y) / 3)
        if data_iter == 0:
            plt.plot(x_display, y[::3], yerr=yerr[::3], label='x', fmt='.', color='C0', alpha=0.5)
            #plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], label = 'y', fmt='.')
            #plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], label = 'z', fmt='.')
        else:
            plt.errorbar(x_display, y[::3], yerr=yerr[::3], fmt='.', color='C0', alpha=0.5)
            #plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], fmt='.')
            #plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], fmt='.')
    else:
        x_display = np.arange(len(resolved_data) / 3)
        if data_iter == 0:
            plt.plot(x_display, resolved_data[::3], label='x', color='C0', alpha=0.5)
            plt.plot(x_display, resolved_data[1::3], label = 'y', color='C1', alpha=0.5)
            plt.plot(x_display, resolved_data[2::3], label = 'z', color='C2', alpha=0.5)
        else:
            plt.plot(x_display, resolved_data[::3], color='C0', alpha=0.5, linestyle=':')
            plt.plot(x_display, resolved_data[1::3], color='C1', alpha=0.5, linestyle=':')
            plt.plot(x_display, resolved_data[2::3], color='C2', alpha=0.5, linestyle=':')

plt.xlabel("Time (Cadences)")
plt.ylabel("Spin (rad/s)")
plt.legend()
plt.savefig("resolved_data.png")
plt.show()
