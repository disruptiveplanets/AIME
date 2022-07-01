# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import matplotlib as mpl

AXIS_SIZE = 12
plt.style.use("jcap")
mpl.rcParams['ytick.labelsize'] = AXIS_SIZE

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

axis_names = {
    "scan-perigee": "$r_p$ (Earth radii)",
    "probe-s-theta": "$\sigma_\\theta$",
    "probe-s-rho": "$\sigma_P / P_\omega$",
    "scan-cadence": "$\Delta t$ (min)",
    "scan-period": "$P_\omega$ (hr)",
    "scan-am": "$a_\\mathcal{A}$ (m)",
    "scan-vex": "$v_\infty$ (km s$^{-1}$)",
    "observation-gap": "$T_\mathrm{gap}$ (hr)",
}
colors = [
    '#f00070', '#db8000', '#00b3cf',
    '#01821d', '#00c22a',
    '#ad8002', '#ffbf0d',
    '#870073', '#e600b3', 
    '#000ac9'
]
trues = {
    "scan-perigee": 5,
    "probe-s-theta": 1e-2,
    "probe-s-rho": 1e-7,
    "scan-cadence": 2,
    "scan-period": 9,
    "scan-am": 1000,
    "scan-vex": 6,
    "observation-gap": 0,
}
logs = {
    "scan-perigee": (False, False),
    "probe-s-theta": (True, True),
    "probe-s-rho": (True, True),
    "scan-cadence": (False, False),
    "scan-period": (False, False),
    "scan-am": (False, True),
    "scan-vex": (False, False),
    "observation-gap": (False, False),
}

specialized_y_labels = {
    "scan-perigee": (None, None, None, None, None, None, None, None, None, None),
    "probe-s-theta": ([1e2, 1e0, 1e-2], [1e1, 1e-1, 1e-3], [1e2, 1e0, 1e-2], [1e0, 1e-2, 1e-4], [1e-1, 1e-3, 1e-5], [1e0, 1e-2, 1e-4], [1e0, 1e-2, 1e-4], [1e1, 1e-1, 1e-3], [1e0, 1e-2, 1e-4], [1e1, 1e-1, 1e-3]),
    "probe-s-rho": (None, None, None, None, None, None, None, None, None, None),
    "scan-cadence": ([100, 0, -100], [10, 0, -10], None, None, None, None, None, None, None, None),
    "scan-period": (None, None, None, None, None, None, None, None, None, None),
    "scan-am": ([10, 6], [1, 0.6], [2, 1], None, None, None, None, None, None, None),
    "scan-vex": (None, None, None, None, None, None, None, None, None, None),
    "observation-gap": (None, None, None, None, None, None, None, None, None, None),
}

thresholds = {
    "scan-perigee":	(6.041275818201059, 18.20618262898073),
    "probe-s-theta": (0.032726900377081515, 0.6113839751941476),
    "scan-cadence":	(15.424913260494899, None),
    "scan-period": (7.591137362332508, None),
    "scan-am": (14.814169963587675, None),
    "scan-vex": (None, None),
    "observation-gap": (None, None),
    "probe-s-rho": (3.6172957141874776e-07, 8.228524489683316e-06),
}


scales = (1e7, 1e2)
LEGEND_SIZE = 12

N_DIM = 10

SCALE_Y = 1.1

FIG_SCALE = 2.2

def show_figs(plot_name, plot_name_index, num_columns):
    percentiles = {}
    name_index = {}
    xs = []

    # Get percentiles
    with open(f"../{plot_name}/percentiles.dat", 'r') as f:
        for line in f.readlines():
            if line == '': continue
            elements = line.split(':')
            name = elements[0]
            perc_array = []
            for percs in elements[1:]:
                perc_array.append([float(x) for x in percs.split(',')])
            perc_array = np.array(perc_array)
            assert(N_DIM == perc_array.shape[0])
            num_percentiles = perc_array.shape[1]
            percentiles[name] = perc_array

    # Get true sigmas
    index = 0
    for name in percentiles.keys():
        dir_name = name[:-14]
        with open(f"../{plot_name}/{dir_name}/{dir_name}.txt", 'r') as f:
            max_j, max_l = f.readline().split(", ")
            max_j, max_l = (int(max_j), int(max_l))
            cadence = float(f.readline())
            perigee = float(f.readline())
            radius = float(f.readline())
            speed = float(f.readline())
            spin = [float(x) for x in f.readline().split(',')]
            jlms = [float(x) for x in f.readline().split(',')]
            theta_true = [float(x) for x in f.readline().split(',')]
            theta_high = [float(x) for x in f.readline().split(',')]
            theta_low = [float(x) for x in f.readline().split(',')]
            sigma = [float(d) for d in f.readline().split(',')]
            try:
                extra_line = float(f.readline())
            except:
                extra_line = None

        name_index[name] = index
        index += 1
        true_sigma = sigma[0]
        if plot_name == "scan-perigee":
            xs.append(perigee)
        elif plot_name == "probe-s-theta":
            xs.append(sigma[0])
        elif plot_name == "probe-s-rho":
            xs.append(sigma[0] * sigma[1])
        elif plot_name == "scan-cadence":
            xs.append(cadence / 60)
        elif plot_name == "scan-am":
            xs.append(radius)
        elif plot_name == "scan-period":
            xs.append(2 * np.pi / np.sqrt(np.dot(spin, spin)) / 3600)
        elif plot_name == "scan-vex":
            xs.append(speed / 1000)
        elif plot_name == "observation-gap":
            xs.append(float(extra_line))
        else:
            raise Exception(f"Plot name {plot_name} not defined")

    xs = np.array(xs)
    with open(f"{plot_name}-x.npy", 'wb') as f:
        np.save(f, xs)

    for plot_index in range(N_DIM):
        offset = 1 if plot_index > 2 else 0
        ax = plt.subplot2grid((31, num_columns), (plot_index * 3 + offset, plot_name_index), rowspan=3)
        param_data = np.zeros(len(xs) * num_percentiles).reshape(num_percentiles, len(xs))
        for f in percentiles.keys():
            param_data[:,name_index[f]] = percentiles[f][plot_index]
        scale = scales[0] if plot_index < 3 else scales[1]


        ax.plot(xs, (param_data[1]-param_data[0]) * scale, color=colors[plot_index], linewidth=1)
        if not logs[plot_name][1]:
            ax.plot(xs, (param_data[-1]-param_data[0]) * scale, color=colors[plot_index], linewidth=1)
        ax.fill_between(xs, (param_data[1]-param_data[0]) * scale, 
            (param_data[-1]-param_data[0]) * scale,  color=colors[plot_index], alpha=0.3)

        ax.plot(xs, (param_data[2]-param_data[0]) * scale, color=colors[plot_index], linewidth=1)
        if not logs[plot_name][1]:
            ax.plot(xs, (param_data[-2]-param_data[0]) * scale, color=colors[plot_index], linewidth=1)
        ax.fill_between(xs, (param_data[2]-param_data[0]) * scale,
            (param_data[-2]-param_data[0]) * scale, color=colors[plot_index], alpha=0.3)

        if not logs[plot_name][1]:
            ax.plot(xs, (param_data[3]-param_data[0]) * scale, color=colors[plot_index], linewidth=1, linestyle='dashed')

        for thresh_index in range(len(thresholds[plot_name])):
            thresh = thresholds[plot_name][thresh_index]
            if thresh is None:
                continue
            if thresh_index == 0:
                ax.axvline(x=thresh, color='r', linewidth=1, linestyle='dashed')
            if thresh_index == 1:
                ax.axvline(x=thresh, color='r', linewidth=1, linestyle='solid')


        if plot_name_index == 0:
            if plot_index < 3:
                ax.set_ylabel(f"$\sigma({param_names[plot_index]}) (\\times 10^{{-7}})$", size=AXIS_SIZE)
            else:
                ax.set_ylabel(f"$\sigma({param_names[plot_index]}) (\\times 10^{{-2}})$", size=AXIS_SIZE)

        if logs[plot_name][0]:
            ax.set_xscale('log')
        if logs[plot_name][1]:
            ax.set_yscale('log')
        else:
            y_min_norm = np.min((param_data[-1]-param_data[0]) * scale)
            y_max_norm = np.max((param_data[1]-param_data[0]) * scale)
            ax.set_ylim(y_min_norm * SCALE_Y, y_max_norm * SCALE_Y)

        ax.set_xlim(np.min(xs), np.max(xs))
        ax.axvline(x=trues[plot_name], color='k', linewidth=1, linestyle='dotted')

        if specialized_y_labels[plot_name][plot_index] is not None:
            ax.set_yticks(ticks=specialized_y_labels[plot_name][plot_index])
            ax.set_yticks(ticks=[], minor=True)

        if plot_index == 9:
            ax.set_xlabel(axis_names[plot_name])
        else:
            ax.set_xticks([])
            

fig = plt.figure(figsize=((11-1.2)*FIG_SCALE, (8.5-1) * FIG_SCALE))

figs_to_show = ["scan-perigee", "scan-vex", "scan-am", "scan-period"]

for i, name in enumerate(figs_to_show):
    show_figs(name, i, len(figs_to_show))

custom_lines = [Line2D([0], [0], color='k', lw=4, alpha=0.3),
                Line2D([0], [0], color='k', lw=4, alpha=0.6),
                Line2D([0], [0], color='k', lw=1, linestyle='dashed')]
fig.legend(custom_lines, ['95\%', '68\%', '50\%'], ncol=3, loc='upper center', prop={'size': LEGEND_SIZE}, bbox_to_anchor=(0.5,0.91))

plt.savefig("all1.pdf", bbox_inches="tight")
plt.savefig("all1.png", bbox_inches="tight")
plt.show()




fig = plt.figure(figsize=((11-1.2)*FIG_SCALE, (8.5-1) * FIG_SCALE))

figs_to_show = ["probe-s-rho", "probe-s-theta", "scan-cadence", "observation-gap"]

for i, name in enumerate(figs_to_show):
    show_figs(name, i, len(figs_to_show))

custom_lines = [Line2D([0], [0], color='k', lw=4, alpha=0.3),
                Line2D([0], [0], color='k', lw=4, alpha=0.6),
                Line2D([0], [0], color='k', lw=1, linestyle='dashed')]
fig.legend(custom_lines, ['95\%', '68\%', '50\%'], ncol=3, loc='upper center', prop={'size': LEGEND_SIZE}, bbox_to_anchor=(0.5,0.91))

plt.savefig("all2.pdf", bbox_inches="tight")
plt.savefig("all2.png", bbox_inches="tight")
plt.show()
