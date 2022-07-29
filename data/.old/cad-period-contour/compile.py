# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.colors import Colormap

plt.style.use("jcap")

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

percentiles = {}
name_index = {}

AXIS_SIZE = 12
LEGEND_SIZE = 12

N_DIM = None
N_PERCENTILES = None

SCALE_Y = 1.1
cadences = []
periods = []

def get_indices_from_name(name):
    return int(name[7:9]), int(name[10:12])

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

for name in percentiles.keys():
    dir_name = name[:12]
    with open(f"{dir_name}/{dir_name}.txt", 'r') as f:
        max_j, max_l = f.readline().split(", ")
        max_j, max_l = (int(max_j), int(max_l))
        cadence = float(f.readline())
        perigee = int(f.readline())
        radius = float(f.readline())
        speed = float(f.readline())
        spin = np.array([float(x) for x in f.readline().split(',')])
        jlms = [float(x) for x in f.readline().split(',')]
        theta_true = [float(x) for x in f.readline().split(',')]
        theta_high = [float(x) for x in f.readline().split(',')]
        theta_low = [float(x) for x in f.readline().split(',')]
        sigma = [float(d) for d in f.readline().split(',')]
    cadence_index, period_index = get_indices_from_name(name)
    if period_index == 0:
        cadences.append(cadence / 60)
    if cadence_index == 0:
        periods.append(2 * np.pi / np.sqrt(np.sum(spin * spin)) / 3600)

cadences = np.sort(cadences)
periods = np.sort(periods)


#fig, axs = plt.subplots(figsize=(6, 19), ncols=1, nrows=10,  sharex=True)
#axs = axs.reshape(-1)
fig = plt.figure(figsize=(6.6, 19))

for plot_index in range(N_DIM):
    offset = 1 if plot_index > 2 else 0
    ax = plt.subplot2grid((31, 1), (plot_index * 3 + offset, 0), rowspan=3)

    param_data = np.zeros((len(cadences), len(periods)))
    for f in percentiles.keys():
        cadence_index, period_index = get_indices_from_name(f)
        sigma = abs(percentiles[f][plot_index][2] - percentiles[f][plot_index][-2]) / 2
        param_data[cadence_index, period_index] = sigma

    scale = 1e5 if plot_index < 3 else 1e2

    #p = ax.pcolormesh(cadences, periods, param_data.transpose() * scale, vmin=0, cmap="Oranges_r")
    levels = np.linspace(0, np.percentile(param_data * scale, 90), 12)
    p = ax.contourf(cadences, periods, param_data.transpose() * scale, cmap="Oranges_r", levels=levels, extend='max')

    cbar = fig.colorbar(p, ax=ax)
    if plot_index < 3:
        cbar.set_label(f"${param_names[plot_index]}\ (\\times 10^{{5}})$")
    else:
        cbar.set_label(f"${param_names[plot_index]}\ (\\times 10^{{2}})$")

    ax.set_ylabel(f"$P_\omega$ (hr)")

    #axs[plot_index].set_xscale('log')
    #axs[plot_index].set_yscale('log')

    if plot_index == 9:
        ax.set_xlabel(f"$\Delta t$ (min)")
    else:
        ax.set_xticks([])
        
plt.savefig("cad-period-contour.pdf", bbox_inches="tight")
plt.savefig("cad-period-contour.png", bbox_inches="tight")

plt.show()
