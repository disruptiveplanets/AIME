import matplotlib.pyplot as plt
import numpy as np

plt.figure(figsize=(8, 4))

f = open("2-params-resolved.dat", 'r')
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

plt.plot(time, perfect_xs, label=f'x spin')
plt.plot(time, perfect_ys, label=f'y spin')
plt.plot(time, perfect_zs, label=f'z spin')
plt.title("Spin")
plt.legend()
plt.figure()

f = open("momentum-energy.dat", 'r')
momentum_xs = []
momentum_ys = []
momentum_zs = []
energies = []
for line in f.readlines():
    if line == "": continue
    if line[0] == "S": continue
    x, y, z, e = line.split(" ")[:-1]
    momentum_xs.append(float(x))
    momentum_ys.append(float(y))
    momentum_zs.append(float(z))
    energies.append(float(e))

time = np.arange(0, len(momentum_xs), 1) * 120/3600

plt.plot(time, momentum_xs, label='x')
plt.plot(time, momentum_ys, label='y')
plt.plot(time, momentum_zs, label='z')
plt.xlim(time[0], time[-1])
plt.legend()
plt.title("Momentum")
plt.tight_layout()
plt.savefig("momentum.pdf")

plt.figure()
plt.plot(time, energies, label='energy')
plt.xlim(time[0], time[-1])
plt.legend()
plt.tight_layout()
plt.savefig("energy.pdf")

plt.show()
