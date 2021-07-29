import asteroids, time, sys
import matplotlib.pyplot as plt
from random_vector import randomize
import numpy as np

CADENCE = 600.0

EARTH_RADIUS = 6370000
EARTH_MASS =5.972e24

spin = [0.00012, 0.00022, 0.00032]
impact_parameter = 5 * EARTH_RADIUS
speed = 4000
jlms = [EARTH_MASS, 0, 6e22, 7e22]
klms = [
    1.2e6, 1.1e5, -4.9e5,
]
SIGMA = 0.2 * np.sqrt(spin[0]**2 + spin[1]**2 + spin[2]**2)
initial_roll = 0

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


y, yerr = randomize(resolved_data, SIGMA)

plt.figure(figsize=(12, 4))
x_display = np.arange(len(y) / 3)
plt.errorbar(x_display, y[::3], yerr=yerr[::3], label = 'x', fmt='.')
plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], label = 'y', fmt='.')
plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], label = 'z', fmt='.')
plt.xlabel("Time (Cadences)")
plt.ylabel("Spin (rad/s)")
plt.legend()
plt.savefig("resolved_data.png")
plt.show()
