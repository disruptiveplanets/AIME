# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import lmfit

plt.style.use("jcap")

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

percentiles = {}
name_index = {}
true_sigma = None
perigees = []

AXIS_SIZE = 12
LEGEND_SIZE = 12

N_DIM = None
N_PERCENTILES = None

SCALE_Y = 1.1

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
index = 0
for name in percentiles.keys():
    dir_name = name[:7]
    with open(f"{dir_name}/{dir_name}.txt", 'r') as f:
        max_j, max_l = f.readline().split(", ")
        max_j, max_l = (int(max_j), int(max_l))
        cadence = int(f.readline())
        perigee = float(f.readline())
        radius = float(f.readline())
        speed = float(f.readline())
        spin = [float(x) for x in f.readline().split(',')]
        jlms = [float(x) for x in f.readline().split(',')]
        theta_true = [float(x) for x in f.readline().split(',')]
        theta_high = [float(x) for x in f.readline().split(',')]
        theta_low = [float(x) for x in f.readline().split(',')]
        sigma = [float(d) for d in f.readline().split(',')]
    name_index[name] = index
    index += 1
    true_sigma = sigma[0]
    perigees.append(perigee)

perigees = np.array(perigees)

fig, axs = plt.subplots(figsize=(14, 8), ncols=3, nrows=4, sharex=True)
axs = axs.reshape(-1)
i = 0

for plot_index in range(N_DIM+1):
    if plot_index == 9:
        continue
    param_data = np.zeros(len(perigees) * N_PERCENTILES).reshape(N_PERCENTILES, len(perigees))
    for f in percentiles.keys():
        param_data[:,name_index[f]] = percentiles[f][i]
    scale = 1e5 if i < 3 else 1

    axs[plot_index].plot(perigees, (param_data[1]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1)
    axs[plot_index].plot(perigees, (param_data[-1]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1)
    axs[plot_index].fill_between(perigees, (param_data[1]-param_data[0]) / true_sigma * scale, 
        (param_data[-1]-param_data[0]) / true_sigma * scale,  color=f"C{i}", alpha=0.3)

    axs[plot_index].plot(perigees, (param_data[2]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1)
    axs[plot_index].plot(perigees, (param_data[-2]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1)
    axs[plot_index].fill_between(perigees, (param_data[2]-param_data[0]) / true_sigma * scale,
        (param_data[-2]-param_data[0]) / true_sigma * scale, color=f"C{i}", alpha=0.3)

    axs[plot_index].plot(perigees, (param_data[3]-param_data[0]) / true_sigma * scale, color=f"C{i}", linewidth=1, linestyle='dashed')

    def fit_func(x, power, slope):
        m = x[:len(x)//2]**power * slope
        return np.append(m, -m)
    model = lmfit.Model(fit_func)
    params = lmfit.Parameters()
    params.add('power', min=0.1, max=5, value=3)
    params.add('slope', min=0, max=100, value=1)
    data = np.append((param_data[1]-param_data[0]) / true_sigma * scale, (param_data[-1]-param_data[0]) / true_sigma * scale)
    result = model.fit(data, params, x=np.append(perigees, perigees), weights=np.ones(len(perigees)*2))
    power, slope = result.params["power"].value, result.params["slope"].value
    power_unc, slope_unc = result.params["power"].stderr, result.params["slope"].stderr
    print(param_names[i])
    print(power, "+/-", power_unc)
    print(slope, "+/-", slope_unc)
    print()

    if i != 9:
        axs[plot_index].plot(perigees, perigees**power * slope, color="k", linewidth=1, linestyle='dotted')
        axs[plot_index].plot(perigees, -perigees**power * slope, color="k", linewidth=1, linestyle='dotted')

    y_min_norm = np.min((param_data[-1]-param_data[0]) / true_sigma * scale)
    y_max_norm = np.max((param_data[1]-param_data[0]) / true_sigma * scale)
    axs[plot_index].set_ylim(y_min_norm * SCALE_Y, y_max_norm * SCALE_Y)

    if i < 3:
        axs[plot_index].set_ylabel(f"$\sigma({param_names[i]}) / \sigma_\\theta (\\times 10^{{-5}})$", size=AXIS_SIZE)
    else:
        axs[plot_index].set_ylabel(f"$\sigma({param_names[i]}) / \sigma_\\theta$", size=AXIS_SIZE)

    #axs[i].set_xscale('log')
    #axs[i].set_yscale('log')

    if plot_index in [6, 8, 10]:
        axs[plot_index].set_xlabel(f"$r_p$ (Earth radii)")
        
    i += 1

axs[9].remove()
axs[11].remove()

custom_lines = [Line2D([0], [0], color='k', lw=4, alpha=0.3),
                Line2D([0], [0], color='k', lw=4, alpha=0.6),
                Line2D([0], [0], color='k', lw=1, linestyle='dashed')]
fig.legend(custom_lines, ['95\%', '68\%', '50\%'], ncol=3, loc='lower right', prop={'size': LEGEND_SIZE})
fig.tight_layout()

line = plt.Line2D([0,1],[0.77, 0.77], transform=fig.transFigure, color="black")
fig.add_artist(line)

plt.savefig("perigee.pdf")
plt.savefig("perigee.png")
plt.show()
