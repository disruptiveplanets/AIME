import numpy as np
import matplotlib.pyplot as plt

f = open("spin-mags.txt", 'r')
start = False
t = 0
times = []
spins = []
dists = []
for line in f.readlines():
    if line == '\n':
        start = True
        continue
    if line == '':
        continue
    if line[0] == 'S':# Last line
        continue
    if start:
        dt, spin_mag, dist = line.split(' ')
        t += float(dt)
        times.append(t)
        spins.append(float(spin_mag))
        dists.append(float(dist))
    
f.close()

fig, ax1 = plt.subplots()

ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Spin mag in non-inertial frame (Hz)")
ax1.plot(times, spins, color='k')
ax1.scatter(times, spins, marker='.', color='k')

ax2 = ax1.twinx()

ax2.set_ylabel("Position (m)")
ax2.plot(times, dists, color='C1')
ax2.scatter(times, dists, color='C1', marker='.')

print("Max spin:", np.max(spins))
print("Min pos:", np.min(dists))

fig.tight_layout()
plt.show()
plt.savefig("img-spin-mags.png")
