import matplotlib.pyplot as plt
import numpy as np

def get_xyz(filename, length=None):
    f = open(filename, 'r')
    xs = []
    ys = []
    zs = []
    for line in f.readlines():
        if line == "": continue
        x, y, z = line.split(" ")
        xs.append(float(x))
        ys.append(float(y))
        zs.append(float(z))
    if length is not None:
        while len(xs) < length:
            xs.append(xs[-1])
            ys.append(ys[-1])
            zs.append(zs[-1])
        while len(xs) > length:
            del xs[-1]
            del ys[-1]
            del zs[-1]
    return np.array([xs, ys, zs])

perfect_resids2 = get_xyz("2-params-resolved-perfect.dat") - get_xyz("2-params-off2-resolved-perfect.dat")
perfect_resids3 = get_xyz("2-params-resolved-perfect.dat") - get_xyz("2-params-off3-resolved-perfect.dat")
resids2 = get_xyz("2-params-resolved.dat", length=len(perfect_resids2[0])) - get_xyz("2-params-off2-resolved.dat", length=len(perfect_resids2[0]))
resids3 = get_xyz("2-params-resolved.dat", length=len(perfect_resids3[0])) - get_xyz("2-params-off3-resolved.dat", length=len(perfect_resids3[0]))
time_resids = get_xyz("2-params-resolved-perfect.dat") - get_xyz("2-params-resolved.dat", length=len(perfect_resids3[0]))
delta2 = np.array([(perfect_resids2 - resids2)[i] / np.max(np.abs(perfect_resids2), axis=1)[i] for i in range(3)])
delta3 = np.array([(perfect_resids3 - resids3)[i] / np.max(np.abs(perfect_resids3), axis=1)[i] for i in range(3)])

time_perfect = np.arange(0, len(perfect_resids3[0]), 1) * 120/3600
time2 = np.arange(0, len(resids2[0]), 1) * 120/3600
time3 = np.arange(0, len(resids3[0]), 1) * 120/3600

plt.figure()
plt.title("Residuals")
plt.plot(time3, time_resids[0], label='x2')
plt.plot(time3, time_resids[1], label='y2')
plt.plot(time3, time_resids[2], label='z2')
plt.legend()

plt.figure()
plt.title("Residual difference")
plt.plot(time3, delta2[0], label='x2')
plt.plot(time3, delta2[1], label='y2')
plt.plot(time3, delta2[2], label='z2')
plt.plot(time3, delta3[0], label='x3')
plt.plot(time3, delta3[1], label='y3')
plt.plot(time3, delta3[2], label='z3')
ylim3 = np.max(np.abs(delta3[:,len(delta3[0]) // 2:])) * 1.05
ylim2 = np.max(np.abs(delta2[:,len(delta2[0]) // 2:])) * 1.05
ylim = max(ylim3, ylim2)
plt.ylim(-ylim, ylim)
plt.legend()
plt.show()