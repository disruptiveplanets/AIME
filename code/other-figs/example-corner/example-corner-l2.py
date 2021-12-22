import sys, corner, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

THETA_TRUE = [0.2, 0.07142857142857142, -0.14285714285714285]

plt.style.use("jcap")
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12

Ixx = 2/3.0 * THETA_TRUE[2] - 4 * THETA_TRUE[1] + 2/3.0
Iyy = 2/3.0 * THETA_TRUE[2] + 4 * THETA_TRUE[1] + 2/3.0
Izz = -4/3.0 * THETA_TRUE[2] + 2/3.0
asquared = 0.9 * (Iyy + Izz - Ixx)
bsquared = 0.9 * (Ixx + Izz - Iyy) / asquared
csquared = 0.9 * (Ixx + Iyy - Izz) / asquared
theta_true = (THETA_TRUE[0], np.sqrt(bsquared), np.sqrt(csquared))

samples = np.loadtxt("param-scan-1.14-0-samples.npy")

labels = ["$\phi_0$", "$b/a$", "$c/a$"]

new_samples = []
for phi, k22, k20 in samples.transpose():
    Ixx = 2/3.0 * k20 - 4 * k22 + 2/3.0
    Iyy = 2/3.0 * k20 + 4 * k22 + 2/3.0
    Izz = -4/3.0 * k20 + 2/3.0
    asquared = 0.9 * (Iyy + Izz - Ixx)
    bsquared = 0.9 * (Ixx + Izz - Iyy) / asquared
    csquared = 0.9 * (Ixx + Iyy - Izz) / asquared
    new_samples.append((phi, np.sqrt(bsquared), np.sqrt(csquared)))

fig = corner.corner(
    np.array(new_samples), labels=labels, truths=theta_true, truth_color=[0.3, 0.3, 1]
)

median = np.percentile(new_samples, 50, axis=0)
top = np.percentile(new_samples, 95, axis=0)
bottom = np.percentile(new_samples, 5, axis=0)

corner.overplot_lines(fig, median, color=[1, 0, 0, 1], linewidth=1)
corner.overplot_lines(fig, top, color=[1, 0, 0, 0.5], linestyle='dotted', linewidth=1)
corner.overplot_lines(fig, bottom, color=[1, 0, 0, 0.5], linestyle='dotted', linewidth=1)

fig.savefig("corner.png")
plt.show()
