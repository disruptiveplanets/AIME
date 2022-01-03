import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(8, 4))

f = open("2-params-resolved-perfect.dat", 'r')
perfect_xs = []
perfect_ys = []
perfect_zs = []
for line in f.readlines():
    if line == "": continue
    x, y, z = line.split(" ")
    perfect_xs.append(float(x))
    perfect_ys.append(float(y))
    perfect_zs.append(float(z))

time = np.arange(0, len(perfect_xs), 1) * 120/3600

#plt.plot(time, perfect_xs, label=f'x perfect')
#plt.plot(time, perfect_ys, label=f'y perfect')
#plt.plot(time, perfect_zs, label=f'z perfect')

f = open("3-params-resolved.dat", 'r')
xs = []
ys = []
zs = []
for line in f.readlines():
    if line == "": continue
    x, y, z = line.split(" ")
    xs.append(float(x))
    ys.append(float(y))
    zs.append(float(z))

while len(xs) < len(perfect_xs):
    xs.append(xs[-1])
    ys.append(ys[-1])
    zs.append(zs[-1])
while len(xs) > len(perfect_xs):
    del xs[-1]
    del ys[-1]
    del zs[-1]

time = np.arange(0, len(zs), 1) * 120/3600

plt.plot(time, xs, label=f'x')
plt.plot(time, ys, label=f'y')
plt.plot(time, zs, label=f'z')

plt.xlim(time[0], time[-1])
plt.legend()
plt.tight_layout()

plt.figure()
plt.title("Residuals")
plt.plot(time, np.array(xs) - np.array(perfect_xs), label="x")
plt.plot(time, np.array(ys) - np.array(perfect_ys),  label="y")
plt.plot(time, np.array(zs) - np.array(perfect_zs), label="z")
plt.legend()
plt.tight_layout()

plt.show()
