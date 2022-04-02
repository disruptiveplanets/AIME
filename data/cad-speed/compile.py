# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

plt.style.use("jcap")

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

sigma_theta = None

AXIS_SIZE = 12
LEGEND_SIZE = 12

N_DIM = 10
N_PERCENTILES = None

SCALE_Y = 1.1
SPEEDS = [0.5, 9, 2]

cadences = []

# Get percentiles
def get_percentiles(include):
    global N_DIM, N_PERCENTILES
    percentiles = {}
    with open("percentiles.dat", 'r') as f:
        for line in f.readlines():
            if line == '': continue
            elements = line.split(':')
            name = elements[0]
            if f"cad-{include}-" not in name:
                continue
            perc_array = []
            for percs in elements[1:]:
                perc_array.append([float(x) for x in percs.split(',')])
            perc_array = np.array(perc_array)
            N_DIM = perc_array.shape[0]
            N_PERCENTILES = perc_array.shape[1]
            percentiles[name] = perc_array
    return percentiles


def select_styles(speed):
    if speed == 9:
        return 'solid'
    if speed == 2:
        return 'dashed'
    return 'dotted'


def show_plot(percentiles, axs, speed):
    global cadences
    # Get true sigmas
    new_cadences = []
    name_index = {}

    index = 0
    for name in percentiles.keys():
        name_index[name] = index
        index += 1

        if "-2-" in name or '-9-' in name:
            dir_name = name[:8]
        else:
            dir_name = name[:10]

        if '-9-' in name:
            continue # There is no file; don't bother to open it


        with open(f"{dir_name}/{dir_name}.txt", 'r') as f:
            max_j, max_l = f.readline().split(", ")
            max_j, max_l = (int(max_j), int(max_l))
            cadence = int(float(f.readline()))
            perigee = int(f.readline())
            radius = float(f.readline())
            v_excess = float(f.readline())
            spin = [float(x) for x in f.readline().split(',')]
            jlms = [float(x) for x in f.readline().split(',')]
            theta_true = [float(x) for x in f.readline().split(',')]
            theta_high = [float(x) for x in f.readline().split(',')]
            theta_low = [float(x) for x in f.readline().split(',')]
            sigma = [float(d) for d in f.readline().split(',')]
            sigma_theta = sigma[0]
            new_cadences.append(cadence / 60)

    if len(new_cadences) > 0:
        cadences = np.array(new_cadences)
    else:
        # Do not replace the cadence array with an empty one, in the case of scan-cadence.
        pass

    assert(len(cadences) > 0)


    for plot_index in range(N_DIM):
        param_data = np.zeros(len(cadences) * N_PERCENTILES).reshape(N_PERCENTILES, len(cadences))
        for f in percentiles.keys():
            param_data[:,name_index[f]] = percentiles[f][plot_index]
        scale = 1/np.max(np.abs(param_data[1]-param_data[0]))

        linestyle = select_styles(speed)
        color = f"C{plot_index}" if speed == 9 else 'k'

        axs[plot_index].plot(cadences, (param_data[1]-param_data[0]) * scale, color=color, linestyle=linestyle, linewidth=2)
        axs[plot_index].plot(cadences, (param_data[-1]-param_data[0]) * scale, color=color, linestyle=linestyle, linewidth=2)

        axs[plot_index].set_ylim(-SCALE_Y, SCALE_Y)

        axs[plot_index].set_ylabel(f"$\sigma({param_names[plot_index]})$", size=AXIS_SIZE)

        if plot_index == 9:
            axs[plot_index].set_xlabel(f"$\Delta t$ (min)")
        else:
            axs[plot_index].set_xticks([])

    custom_lines = []
    labels = []
    for speed in SPEEDS:
        linestyle = select_styles(speed)
        custom_lines.append(Line2D([0], [0], color='k', lw=1, linestyle=linestyle))
        labels.append(f"$t_\\mathrm{{spin}}/t_\\mathrm{{orbit}}={speed}$ hr")


    fig.legend(custom_lines, labels, ncol=3, loc='upper center', prop={'size': LEGEND_SIZE}, bbox_to_anchor=(0.5,0.91))



fig = plt.figure(figsize=(6.6, 19))
axes = []
for plot_index in range(N_DIM):
    offset = 1 if plot_index > 2 else 0
    axes.append(plt.subplot2grid((31, 1), (plot_index * 3 + offset, 0), rowspan=3))

for p in SPEEDS:
    show_plot(get_percentiles(p), axes, p)

plt.savefig(f"cad-speed.pdf", bbox_inches="tight")
plt.savefig(f"cad-speed.png", bbox_inches="tight")

plt.show()