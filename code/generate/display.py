import matplotlib.pyplot as plt
import numpy as np

FILE_NAME = "2-params-resolved.dat"
f = open(FILE_NAME, 'r')

xs = []
ys = []
zs = []
for line in f.readlines():
    if line == "": continue
    x, y, z = line.split(" ")
    xs.append(float(x))
    ys.append(float(y))
    zs.append(float(z))

time = np.arange(0, len(zs), 1) * 120/3600

plt.figure(figsize=(8, 4))

plt.plot(time, xs, label='x')
plt.plot(time, ys, label='y')
plt.plot(time, zs, label='z')
plt.xlim(time[0], time[-1])

plt.legend()
plt.tight_layout()
plt.show()
