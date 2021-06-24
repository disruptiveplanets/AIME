import os
from matplotlib import pyplot as plt
import numpy as np

raw_files = os.listdir()
spinxs = []
spinys = []
spinzs = []
files = []
maxlen = 0
for filename in raw_files:
    if filename[-4:] == '.txt' and filename[0] != '_':
        files.append(filename)
for filename in files:
    f = open(filename, 'r')
    xs = []
    ys = []
    zs = []
    for line in f.readlines():
        x, y, z = line[1:-2].split(' ')
        xs.append(float(x))
        ys.append(float(y))
        zs.append(float(z))

    maxlen = max(maxlen, len(xs))
    spinxs.append(xs)
    spinys.append(ys)
    spinzs.append(zs)

for i in range(len(spinxs)):
    buffer = (maxlen - len(spinxs[i])) / 2
    xs = np.arange(buffer, maxlen - buffer)
    plt.plot(xs, spinxs[i], c='C'+str(i), linestyle='solid', label=files[i])
    plt.plot(xs, spinys[i], c='C'+str(i), linestyle='dashed')
    plt.plot(xs, spinzs[i], c='C'+str(i), linestyle='dotted')

plt.xlabel("Time index")
plt.ylabel("Spin components [Hz]")
plt.legend()
plt.show()
