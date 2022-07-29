import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) == 2:
    file_name = sys.argv[-1]
else:
    file_name = "2-params"

f = open(f"{file_name}-resolved.dat", 'r')
xs = []
ys = []
zs = []
for line in f.readlines():
    if line == "": continue
    x, y, z = line.split(" ")
    xs.append(float(x))
    ys.append(float(y))
    zs.append(float(z))
xs = np.array(xs)
ys = np.array(ys)
zs = np.array(zs)

plt.figure(figsize=(8, 8))
plt.plot(xs[zs>0], ys[zs>0])
plt.plot(xs[zs<=0], ys[zs<=0], linestyle='dashed')
plt.scatter(0, 0)
plt.scatter(xs[0], ys[0])
plt.axis('equal')
plt.xlim(-0.45e9, 0.45e9)
plt.ylim(-0.45e9, 0.45e9)
plt.show()
