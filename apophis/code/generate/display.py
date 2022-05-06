import matplotlib.pyplot as plt
import numpy as np
import sys, os

plt.style.use("jcap")

if len(sys.argv) == 2:
    file_name = sys.argv[-1]
else:
    file_name = "2-params"

plt.figure(figsize=(6, 6))

if os.path.exists(f"{file_name}-resolved-perfect.dat"):
    f = open(f"{file_name}-resolved-perfect.dat", 'r')
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

    plt.plot(time, perfect_xs, label=f'x perfect')
    plt.plot(time, perfect_ys, label=f'y perfect')
    plt.plot(time, perfect_zs, label=f'z perfect')
    plt.scatter(time, perfect_xs, label=f'x perfect', s=2)
    plt.scatter(time, perfect_ys, label=f'y perfect', s=2)
    plt.scatter(time, perfect_zs, label=f'z perfect', s=2)
    plt.title(file_name)

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

if os.path.exists(f"{file_name}-resolved-perfect.dat"):
    while len(xs) < len(perfect_xs):
        xs.append(xs[-1])
        ys.append(ys[-1])
        zs.append(zs[-1])
    while len(xs) > len(perfect_xs):
        del xs[-1]
        del ys[-1]
        del zs[-1]

time = np.arange(-len(zs)/2, len(zs)/2, 1) * 120/3600

plt.plot(time, np.array(xs)*3600, label=f'$\omega_x$')
plt.plot(time, np.array(ys)*3600, label=f'$\omega_y$')
plt.plot(time, np.array(zs)*3600, label=f'$\omega_z$')

plt.xlabel("Time after perigee (hr)")
plt.ylabel("Angular velocity (rad/hr)")

plt.xlim(time[0], time[-1])
plt.legend()
plt.tight_layout()

if os.path.exists(f"{file_name}-resolved-perfect.dat"):
    plt.figure()
    plt.title(f"Residuals: {file_name}")
    plt.plot(time, np.array(xs) - np.array(perfect_xs), label="x")
    plt.plot(time, np.array(ys) - np.array(perfect_ys),  label="y")
    plt.plot(time, np.array(zs) - np.array(perfect_zs), label="z")
    plt.legend()
    plt.tight_layout()

plt.savefig("spin-data.png")
#plt.savefig("spin-data.pdf")

plt.show()
