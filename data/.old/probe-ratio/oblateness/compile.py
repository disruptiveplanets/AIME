# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

plt.style.use("jcap")

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

oblateness_markers = {
    "Moon": 203e-6,
    "Earth": 1082e-6,
    #"Mars": 1960e-6,
    "Jupiter": 14696e-6,
    "Neptune": 3343e-6,}

percentiles = {}
name_index = {}
true_sigma = []
j2s = []

AXIS_SIZE = 12
LEGEND_SIZE = 12

N_DIM = None
N_PERCENTILES = None

def plot_best_fit(ax, xs, ys, scale):
    slope = np.cov(ys, xs)[0][1] / np.var(xs)
    yint = np.mean(ys) - slope * np.mean(xs)
    ax.plot(xs, (xs * slope + yint) * scale, color='k', linewidth=1, linestyle='dotted')
    return slope / yint

# Get percentiles
with open("percentiles.dat", 'r') as f:
    for line in f.readlines():
        if line == '': continue
        elements = line.split(':')
        name = elements[0]
        perc_array = []
        for percs in elements[1:]:
            perc_array.append([float(x) for x in percs.split(',')])
        perc_array = np.array(perc_array)
        N_DIM = perc_array.shape[0]
        N_PERCENTILES = perc_array.shape[1]
        percentiles[name] = perc_array

# Get true sigmas
for name in percentiles.keys():
    dir_name = name[:7]
    with open(f"{dir_name}/{dir_name}.txt", 'r') as f:
        max_j, max_l = f.readline().split(", ")
        max_j, max_l = (int(max_j), int(max_l))
        cadence = int(f.readline())
        perigee = int(f.readline())
        radius = float(f.readline())
        speed = float(f.readline())
        spin = [float(x) for x in f.readline().split(',')]
        jlms = [float(x) for x in f.readline().split(',')]
        theta_true = [float(x) for x in f.readline().split(',')]
        theta_high = [float(x) for x in f.readline().split(',')]
        theta_low = [float(x) for x in f.readline().split(',')]
        sigma = [float(d) for d in f.readline().split(',')]
    name_index[name] = int(dir_name[-2:])
    true_sigma = sigma[0]
    j2s.append(jlms[-1])

oblatenesses = -np.array(j2s) / 2

fig, axs = plt.subplots(figsize=(5, 7), ncols=1, nrows=3, sharex=True)
axs = axs.reshape(-1)

for i in range(N_DIM):
    param_data = np.zeros(len(oblatenesses) * N_PERCENTILES).reshape(N_PERCENTILES, len(oblatenesses))
    for f in percentiles.keys():
        param_data[:,name_index[f]] = percentiles[f][i]
    scale = 10**5 if i < 3 else 1
    axs[i].plot(oblatenesses, (param_data[1]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1)
    axs[i].plot(oblatenesses, (param_data[-1]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1)
    axs[i].fill_between(oblatenesses, (param_data[1]-param_data[0]) / true_sigma * scale, 
        (param_data[-1]-param_data[0]) / true_sigma * scale,  color=f"C{i}", alpha=0.3)

    axs[i].plot(oblatenesses, (param_data[2]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1)
    axs[i].plot(oblatenesses, (param_data[-2]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1)
    axs[i].fill_between(oblatenesses, (param_data[2]-param_data[0]) / true_sigma * scale,
        (param_data[-2]-param_data[0]) / true_sigma * scale, color=f"C{i}", alpha=0.3)

    axs[i].plot(oblatenesses, (param_data[3]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1, linestyle='dashed')

    axs[i].set_xscale('log')
    axs[i].set_ylabel(f"$\sigma({param_names[i]}) / \sigma_\\theta\\ (\\times 10^{{-5}})$", size=AXIS_SIZE)

    s1n = plot_best_fit(axs[i], oblatenesses, (param_data[-1]-param_data[0]) / true_sigma, scale)
    s1 = plot_best_fit(axs[i], oblatenesses, (param_data[1]-param_data[0]) / true_sigma, scale)
    s2n = plot_best_fit(axs[i], oblatenesses, (param_data[-2]-param_data[0]) / true_sigma, scale)
    s2 = plot_best_fit(axs[i], oblatenesses, (param_data[2]-param_data[0]) / true_sigma, scale)
    #sm = plot_best_fit(axs[i], oblatenesses, (param_data[3]-param_data[0]) / true_sigma, scale)
    
    print((s1n + s1 + s2n + s2) / 4)

    for ob in oblateness_markers.values():
        axs[i].axvline(x=ob, color='k', linewidth=1)

    if i == 0:
        for name, ob in oblateness_markers.items():
            axs[i].text(x=ob, y=17, verticalalignment='bottom', horizontalalignment='center', s=name, fontsize=12)

    if i == 2:
        axs[i].set_xlabel(f"$\epsilon$")

custom_lines = [Line2D([0], [0], color='k', lw=4, alpha=0.3),
                Line2D([0], [0], color='k', lw=4, alpha=0.6),
                Line2D([0], [0], color='k', lw=1, linestyle='dashed')]
fig.legend(custom_lines, ['95\%', '68\%', '50\%'], ncol=3, framealpha=1, loc="lower center", prop={'size': LEGEND_SIZE})
fig.tight_layout()
fig.subplots_adjust(bottom=0.15)
plt.savefig("oblateness.pdf")
plt.savefig("oblateness.png")
plt.show()