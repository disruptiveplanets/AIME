import numpy as np
import matplotlib.pyplot as plt

f = open("torque-mags.txt", 'r')
start = False
times = []
spins = []
for line in f.readlines():
    if line == '\n':
        start = True
        continue
    if line == '':
        continue
    if line[0] in ['S', 'M']:# Last line
        continue
    if start:
        spin_mag, time = line.split(' ')
        times.append(float(time))
        spins.append(float(spin_mag))
    
f.close()

fig, ax1 = plt.subplots()

ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Spin mag in non-inertial frame (Hz)")
ax1.plot(times, spins, color='k')
ax1.scatter(times, spins, marker='.', color='k')

fig.tight_layout()
plt.show()
plt.savefig("img-torque-mags.png")
