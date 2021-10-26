import sys, corner, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

THETA_TRUE = [0.1, 0.0412789, -0.1, 0.024187, 0.0789512, -0.0521378, 0.058931, -0.013897, -0.087923149, 0.089423]

plt.style.use("jcap")
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12


samples = np.loadtxt("small-sigma-3-0-samples.dat")

labels = ["$\phi_0$", "$K_{22}$", "$K_{20}$", "$\\mathrm{Re}K_{33}$", "$\\mathrm{Im}K_{33}$",
    "$\\mathrm{Re}K_{32}$", "$\\mathrm{Im}K_{32}$", "$\\mathrm{Re}K_{31}$", "$\\mathrm{Im}K_{31}$",
    "$K_{30}$"]

fig = corner.corner(
    samples.transpose(), labels=labels, truths=THETA_TRUE, truth_color=[0.3, 0.3, 1]
)

median = np.percentile(samples, 50, axis=1)
top = np.percentile(samples, 95, axis=1)
bottom = np.percentile(samples, 5, axis=1)
print(median.shape)

corner.overplot_lines(fig, median, color=[1, 0, 0, 1], linewidth=1)
corner.overplot_lines(fig, top, color=[1, 0, 0, 0.5], linestyle='dotted', linewidth=1)
corner.overplot_lines(fig, bottom, color=[1, 0, 0, 0.5], linestyle='dotted', linewidth=1)

fig.savefig("corner-l3.png")
plt.show()
