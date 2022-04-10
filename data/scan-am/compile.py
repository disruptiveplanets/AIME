# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

plt.style.use("jcap")

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

percentiles = {}
name_index = {}
true_sigma = None
am = []

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
    dir_name = name[:5]
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
    name_index[name] = index
    index += 1
    true_sigma = sigma[0]
    am.append(radius)

am = np.array(am)

fig = plt.figure(figsize=(6.6, 19))

for plot_index in range(N_DIM):
    offset = 1 if plot_index > 2 else 0
    ax = plt.subplot2grid((31, 1), (plot_index * 3 + offset, 0), rowspan=3)
    param_data = np.zeros(len(am) * N_PERCENTILES).reshape(N_PERCENTILES, len(am))
    for f in percentiles.keys():
        param_data[:,name_index[f]] = percentiles[f][plot_index]
    scale = 1e7 if plot_index < 3 else 1e2

    ax.plot(am, (param_data[1]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1)
    ax.plot(am, (param_data[-1]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1)
    ax.fill_between(am, (param_data[1]-param_data[0]) * scale, 
        (param_data[-1]-param_data[0]) * scale,  color=f"C{plot_index}", alpha=0.3)

    ax.plot(am, (param_data[2]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1)
    ax.plot(am, (param_data[-2]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1)
    ax.fill_between(am, (param_data[2]-param_data[0]) * scale,
        (param_data[-2]-param_data[0]) * scale, color=f"C{plot_index}", alpha=0.3)

    #ax.plot(am, (param_data[3]-param_data[0]) * scale, color=f"C{plot_index}", linewidth=1, linestyle='dashed')

    y_min_norm = np.min((param_data[-1]-param_data[0]) * scale)
    y_max_norm = np.max((param_data[1]-param_data[0]) * scale)
    #ax.set_ylim(y_min_norm * SCALE_Y, y_max_norm * SCALE_Y)

    thresh = am[(np.abs(param_data[2]-param_data[0]) < 0.01) & np.abs((param_data[-2]-param_data[0]) < 0.01)]
    if len(thresh) > 0:
        ax.axvline(x=thresh[0], color='r', linewidth=1)

    if plot_index < 3:
        ax.set_ylabel(f"$\sigma({param_names[plot_index]}) (\\times 10^{{-7}})$", size=AXIS_SIZE)
        ax.set_ylim(0.1, y_max_norm * 2)
    else:
        ax.set_ylabel(f"$\sigma({param_names[plot_index]}) (\\times 10^{{-2}})$", size=AXIS_SIZE)

    #ax.set_xscale('log')
    ax.set_yscale('log')

    def model_fn(x, p, a):
        y0 = np.interp(1000, x, (param_data[2]-param_data[0])*scale)
        return y0 * np.exp(a * (1000**p-x**p))

    import lmfit
    model = lmfit.Model(model_fn)
    params = lmfit.Parameters()
    params.add("p", min=0, max=1, value=0.25)
    params.add("a", min=0, value=1)
    res = model.fit(x=am, data=(param_data[2]-param_data[0])*scale, params=params, weights=np.ones_like(am))
    #print(res.fit_report())
    print(f"{res.params['p'].value} +/- {res.params['p'].stderr}\t{res.params['a'].value} +/- {res.params['a'].stderr}")
    #ax.plot(am, model_fn(am, res.params["p"].value, res.params["a"].value), color='k', linestyle='dotted', linewidth=1)

    ax.set_xlim(np.min(am), np.max(am))
    ax.axvline(x=1000, color='k', linewidth=1, linestyle='dashed')

    if plot_index == 9:
        ax.set_xlabel(f"$a_m$ (m)")
    else:
        ax.set_xticks([])

custom_lines = [Line2D([0], [0], color='k', lw=4, alpha=0.3),
                Line2D([0], [0], color='k', lw=4, alpha=0.6)]
fig.legend(custom_lines, ['95\%', '68\%', '50\%'], ncol=3, loc='upper center', prop={'size': LEGEND_SIZE}, bbox_to_anchor=(0.5,0.91))

plt.savefig("am.pdf", bbox_inches="tight")
plt.savefig("am.png", bbox_inches="tight")
plt.show()
