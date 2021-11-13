# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
import os

plt.style.use("jcap")

names_center = []
names_edge = []
sigmas_center = []
sigmas_edge = []
file_names = os.listdir()
file_names.sort()
true_center = None
true_edge = None

param_names = ["$\alpha", "$K_{20}$", "$K_{22}$", "Re$K_{33}$", "Im$K_{33}$", "Re$K_{32}$", "Im$K_{32}$", "Re$K_{31}$", "Im$K_{31}$", "$K_{30}$"]

DIFF_THRESHOLD = [0.01, 0.01, 0.01]
SELECT_INDEX = 1

# Fill points
for bare in file_names:
    if not os.path.isdir(bare):
        continue
    f = open("{0}/{0}.txt".format(bare), 'r')
    max_j, max_l = f.readline().split(", ")
    max_j, max_l = (int(max_j), int(max_l))
    num_fits = [int(i) for i in f.readline().split(', ')]
    cadence = int(f.readline())
    impact_parameter = int(f.readline())
    radius = float(f.readline())
    speed = float(f.readline())
    spin = [float(x) for x in f.readline().split(',')]
    jlms = [float(x) for x in f.readline().split(',')]
    theta_true = [float(x) for x in f.readline().split(',')]
    theta_high = [float(x) for x in f.readline().split(',')]
    theta_low = [float(x) for x in f.readline().split(',')]
    sigma = float(f.readline())
    f.close()
    if abs(theta_true[1]) < 0.001:
        names_center.append(bare)
        sigmas_center.append(sigma)
        true_center = theta_true
    else:
        names_edge.append(bare)
        sigmas_edge.append(sigma)
        true_edge = theta_true
sigmas_center = np.array(sigmas_center)
sigmas_edge = np.array(sigmas_edge)

data_center = [None] * len(sigmas_center)
data_edge = [None] * len(sigmas_edge)
diff_center = [None] * len(sigmas_center)
diff_edge = [None] * len(sigmas_edge)
f = open("percentiles.dat", 'r')
for line in f.readlines():
    if line == "":
        continue
    line = line.split(": ")
    data = np.array([[float(f) for f in l.split(", ")] for l in line[1:]])
    name = line[0][:-14]
    if name in names_center:
        i = names_center.index(name)
        dist = np.abs(data[SELECT_INDEX][0] - true_center[SELECT_INDEX])
        if dist < DIFF_THRESHOLD[SELECT_INDEX]:
            if diff_center[i] is None or diff_center[i] > dist:
                data_center[i] = data
                diff_center[i] = dist
    elif name in names_edge:
        i = names_edge.index(name)
        dist = np.abs(data[SELECT_INDEX][0] - true_edge[SELECT_INDEX])
        if dist < DIFF_THRESHOLD[SELECT_INDEX]:
            if diff_edge[i] is None or diff_edge[i] > dist:
                data_edge[i] = data
                diff_edge[i] = dist
    else:
        raise Exception(f"Could not find name {name}")


# Plot pdfs
def plot_sigmas(names, sigmas, data, label_name, plot_index, color):
    very_high = [np.nan] * len(names)
    high = [np.nan] * len(names)
    med = [np.nan] * len(names)
    low = [np.nan] * len(names)
    very_low = [np.nan] * len(names)
    for i, name in enumerate(names):
        if data[i] is None: continue
        min_dist = None
        min_file = None
        very_high[i] = data[i][plot_index][1]
        high[i] = data[i][plot_index][2]
        med[i] = data[i][plot_index][3]
        low[i] = data[i][plot_index][4]
        very_low[i] = data[i][plot_index][5]

    very_high = np.array(very_high)
    high = np.array(high)
    med = np.array(med)
    low = np.array(low)
    very_low = np.array(very_low)

    plt.plot(sigmas, very_high-med, color=color, linewidth=1, linestyle='dashed')
    plt.plot(sigmas, high-med, color=color, linewidth=1)
    #plt.plot(sigmas, med, color=color)
    plt.plot(sigmas, low-med, color=color, linewidth=1)
    plt.plot(sigmas, very_low-med, color=color, linewidth=1, linestyle='dashed')
    plt.fill_between(sigmas, low-med, high-med, color=color, alpha=0.5, label=label_name)


if __name__ == "__main__":
    for i in range(3):#len(true_center)):
        plt.figure()
        plot_sigmas(names_center, sigmas_center, data_center, "center", i, "C0")
        plot_sigmas(names_edge, sigmas_edge, data_edge, "edge", i, "C1")
        plt.xlabel("$\sigma_\\theta$")
        plt.xscale('log')
        plt.xscale('log')
        plt.ylabel("$\sigma_{}$".format(i))
        plt.legend()
        plt.tight_layout()
        plt.savefig("sigmas-{}.png".format(i))
    plt.show()
