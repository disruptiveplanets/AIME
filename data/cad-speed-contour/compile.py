# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl

plt.style.use("jcap")
mpl.rcParams["font.size"] = 14

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

percentiles = {}
name_index = {}

AXIS_SIZE = 12
LEGEND_SIZE = 12

N_DIM = None
N_PERCENTILES = None

SCALE_Y = 1.1
cadences = []
time_ratios = []

def get_indices_from_name(name):
    return int(name[7:9]), int(name[10:12])

with open("../../code/thresholds/cad-speed-contour.npy", 'rb') as f:
   uncs = np.load(f)[:,:,2] # Use median
print(uncs)

# Adjust uncs so that everything after the first large unc is large
def adjust(a):
    res = []
    for line in a:
        wheres = np.where(line > 1)[0]
        if len(wheres) == 0 or wheres[0] == len(line) - 1:
            res.append(line)
        else:
            new_line = np.copy(line)
            new_line[wheres[0]+1:] = 2
            res.append(new_line)
    return np.array(res)
#uncs = adjust(uncs)


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
        time_ratio = float(f.readline())
    cadence_index, period_index = get_indices_from_name(name)
    if period_index == 1:
        cadences.append(cadence / 60)
    if cadence_index == 1:
        time_ratios.append(time_ratio)

cadences = np.sort(cadences)
time_ratios = np.sort(time_ratios)


#fig, axs = plt.subplots(figsize=(6, 19), ncols=1, nrows=10,  sharex=True)
#axs = axs.reshape(-1)
fig = plt.figure(figsize=(6.6, 19))

for plot_index in range(N_DIM):
    offset = 1 if plot_index > 2 else 0
    ax = plt.subplot2grid((31, 1), (plot_index * 3 + offset, 0), rowspan=3)

    param_data = np.zeros((len(cadences), len(time_ratios)))
    for f in percentiles.keys():
        cadence_index, period_index = get_indices_from_name(f)
        sigma = abs(percentiles[f][plot_index][2] - percentiles[f][plot_index][-2]) / 2
        param_data[cadence_index, period_index] = sigma

    scale = 1e5 if plot_index < 3 else 1e2

    #p = ax.pcolormesh(cadences, time_ratios, param_data.transpose() * scale, vmin=0, cmap="Oranges_r")
    levels = np.linspace(0, np.percentile(param_data * scale, 90), 12)
    p = ax.contourf(cadences, time_ratios, param_data.transpose() * scale, cmap="PuBu_r", levels=levels, extend='max')
    ax.contour(cadences, time_ratios, uncs, colors='r', levels=[1e-4, 1e-3], linestyles=["dashed", "solid"])

    cbar = fig.colorbar(p, ax=ax)
    if plot_index < 3:
        cbar.set_label(f"$\sigma({param_names[plot_index]})\ (\\times 10^{{5}})$")
    else:
        cbar.set_label(f"$\sigma({param_names[plot_index]})\ (\\times 10^{{2}})$")

    ax.set_ylabel(f"$t_\\mathrm{{spin}}/t_\\mathrm{{orbit}}$")


    #axs[plot_index].set_xscale('log')
    #axs[plot_index].set_yscale('log')

    ax.axhline(y=1, color='k', linestyle='dashed', linewidth=1)

    if plot_index == 9:
        ax.set_xlabel(f"$\Delta t$ (min)")
    else:
        ax.set_xticks([])

    if plot_index != 0 and plot_index != 3:
        ax.set_yticklabels([0.5, 1.0, 1.5])
        ax.set_yticks([0.5, 1.0, 1.5])
        
plt.savefig("cad-speed-contour.pdf", bbox_inches="tight")
plt.savefig("cad-speed-contour.png", bbox_inches="tight")

# plt.figure()
# plt.pcolormesh(cadences, periods, uncs)

plt.show()
