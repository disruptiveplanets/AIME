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

fig = plt.figure(figsize=(6.6, 19))

for plot_index in range(N_DIM):
    offset = 1 if plot_index > 2 else 0
    ax = plt.subplot2grid((31, 1), (plot_index * 3 + offset, 0), rowspan=3)
    param_data = np.zeros(len(perigees) * N_PERCENTILES).reshape(N_PERCENTILES, len(perigees))
    for f in percentiles.keys():
        param_data[:,name_index[f]] = percentiles[f][plot_index]
    if plot_index < 1:
        scale = 1e5
    elif plot_index < 3:
        scale = 1e3
    else:
        scale = 1

    ax.plot(perigees, (param_data[1]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1)
    ax.plot(perigees, (param_data[-1]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1)
    ax.fill_between(perigees, (param_data[1]-param_data[0]) * scale, 
        (param_data[-1]-param_data[0]) * scale,  color=f"C{plot_index}", alpha=0.3)

    ax.plot(perigees, (param_data[2]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1)
    ax.plot(perigees, (param_data[-2]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1)
    ax.fill_between(perigees, (param_data[2]-param_data[0]) * scale,
        (param_data[-2]-param_data[0]) * scale, color=f"C{plot_index}", alpha=0.3)

    ax.plot(perigees, (param_data[3]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1, linestyle='dashed')

    def fit_func(x, power, slope):
        m = x[:len(x)//2]**power * slope
        return np.append(m, -m)
    model = lmfit.Model(fit_func)
    params = lmfit.Parameters()
    params.add('power', min=1, max=8, value=3)
    params.add('slope', min=0, max=100, value=1)

    # Mask all the data that is after the first point with sigma2 > 0.7
    overages = np.where(param_data[1]-param_data[0] > 0.7)[0]
    if len(overages) == 0:
        data_mask = np.arange(len(param_data[1]))
    else:
        data_mask = np.arange(overages[0])

    data = np.append((param_data[1]-param_data[0])[data_mask] * scale, (param_data[-1]-param_data[0])[data_mask] * scale)
    result = model.fit(data, params, x=np.append(perigees[data_mask], perigees[data_mask]), weights=np.ones(len(perigees[data_mask])*2))
    power, slope = result.params["power"].value, result.params["slope"].value
    power_unc, slope_unc = result.params["power"].stderr, result.params["slope"].stderr
    print(f"${param_names[plot_index]}$ & ${power} \pm {power_unc}$ \\\\")
    #print(slope, "+/-", slope_unc)

    ax.plot(perigees, perigees**power * slope, color="k", linewidth=1, linestyle='dotted')
    ax.plot(perigees, -perigees**power * slope, color="k", linewidth=1, linestyle='dotted')

    y_min_norm = np.min((param_data[-1]-param_data[0]) * scale)
    y_max_norm = np.max((param_data[1]-param_data[0]) * scale)
    ax.set_ylim(y_min_norm * SCALE_Y, y_max_norm * SCALE_Y)

    thresh = perigees[(np.abs(param_data[2]-param_data[0]) > 0.01) | np.abs((param_data[-2]-param_data[0]) > 0.01)]
    if len(thresh) > 0:
        ax.axvline(x=thresh[0], color='r', linewidth=1)

    if plot_index < 1:
        ax.set_ylabel(f"$\sigma({param_names[plot_index]}) (\\times 10^{{-5}})$", size=AXIS_SIZE)
    elif plot_index < 3:
        ax.set_ylabel(f"$\sigma({param_names[plot_index]}) (\\times 10^{{-3}})$", size=AXIS_SIZE)
    else:
        ax.set_ylabel(f"$\sigma({param_names[plot_index]})$", size=AXIS_SIZE)

    #ax.set_xscale('log')
    #ax.set_yscale('log')
    ax.set_xlim(np.min(perigees), np.max(perigees))
    ax.axvline(x=5, color='k', linewidth=1, linestyle='dashed')

    if plot_index == 9:
        ax.set_xlabel(f"$r_p$ (Earth radii)")
    else:
        ax.set_xticks([])
        
custom_lines = [Line2D([0], [0], color='k', lw=4, alpha=0.3),
                Line2D([0], [0], color='k', lw=4, alpha=0.6),
                Line2D([0], [0], color='k', lw=1, linestyle='dashed')]
fig.legend(custom_lines, ['95\%', '68\%', '50\%'], ncol=3, loc='upper center', prop={'size': LEGEND_SIZE}, bbox_to_anchor=(0.5,0.91))

plt.savefig("perigee.pdf", bbox_inches="tight")
plt.savefig("perigee.png", bbox_inches="tight")
plt.show()