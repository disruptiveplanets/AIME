import matplotlib.pyplot as plt
import numpy as np
import sys

if len(sys.argv) == 2:
    file_name = sys.argv[-1]
else:
    file_name = "2-params"

plt.figure(figsize=(8, 4))

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
time = np.arange(0, len(xs), 1) * 120/3600

plt.plot(time, xs, label=f'x')
plt.plot(time, ys, label=f'y')
plt.plot(time, zs, label=f'z')

plt.xlabel("Time (hr)")
plt.ylabel("Angular velocity (1/sec)")


plt.xlim(time[0], time[-1])
plt.legend()
plt.tight_layout()


plt.figure()
period = 2 * np.pi / np.sqrt(np.array(xs)**2 + np.array(ys)**2 + np.array(zs)**2) / 3600
plt.plot(time, period)
plt.xlabel("Time (hr)")
plt.ylabel("Period (hr)")
plt.tight_layout()
plt
plt.show()
