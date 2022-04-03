# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

plt.style.use("jcap")

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

percentiles = {}
name_index = {}
true_sigmas = []

AXIS_SIZE = 12
LEGEND_SIZE = 12

N_DIM = None
N_PERCENTILES = None


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
    dir_name = name[:10]
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
    true_sigmas.append(sigma[0])
    ratio = sigma[1]

true_sigmas = np.array(true_sigmas)

fig = plt.figure(figsize=(6.6, 19))

for plot_index in range(N_DIM):
    offset = 1 if plot_index > 2 else 0
    ax = plt.subplot2grid((31, 1), (plot_index * 3 + offset, 0), rowspan=3)

    product = true_sigmas**2 * ratio
    param_data = np.zeros(len(product) * N_PERCENTILES).reshape(N_PERCENTILES, len(product))
    for f in percentiles.keys():
        param_data[:,name_index[f]] = percentiles[f][plot_index]
    scale = 10**5 if plot_index < 3 else 1

    ax.plot(product, (param_data[1]-param_data[0]) / true_sigmas * scale, color=f"C{plot_index}", linewidth=1)
    ax.plot(product, (param_data[-1]-param_data[0]) / true_sigmas * scale, color=f"C{plot_index}", linewidth=1)
    ax.fill_between(product, (param_data[1]-param_data[0]) / true_sigmas * scale, 
        (param_data[-1]-param_data[0]) / true_sigmas * scale,  color=f"C{plot_index}", alpha=0.3)

    ax.plot(product, (param_data[2]-param_data[0]) / true_sigmas * scale, color=f"C{plot_index}", linewidth=1)
    ax.plot(product, (param_data[-2]-param_data[0]) / true_sigmas * scale, color=f"C{plot_index}", linewidth=1)
    ax.fill_between(product, (param_data[2]-param_data[0]) / true_sigmas * scale,
        (param_data[-2]-param_data[0]) / true_sigmas * scale, color=f"C{plot_index}", alpha=0.3)

    ax.plot(product, (param_data[3]-param_data[0]) / true_sigmas * scale, color=f"C{plot_index}", linewidth=1, linestyle='dashed')

    if plot_index < 3:
        ax.set_ylabel(f"$\sigma({param_names[plot_index]}) / \sigma_\\theta\\ (\\times 10^{{-5}})$", size=AXIS_SIZE)
    else:
        ax.set_ylabel(f"$\sigma({param_names[plot_index]}) / \sigma_\\theta$", size=AXIS_SIZE)

    ax.set_xscale('log')

    #constant_val = (param_data[1]-param_data[0])[0] * scale / true_sigmas**2
    #ax.plot(product, np.ones_like(param_data[1]) * constant_val * true_sigmas, color='k', linewidth=1, linestyle='dotted')

    ax.axvline(x=1e-9, color='k', linewidth=1, linestyle='dashed')
    ax.set_xlim(np.min(product), np.max(product))

    thresh = product[(np.abs(param_data[2]-param_data[0]) > 0.01) | np.abs((param_data[-2]-param_data[0]) > 0.01)]
    #print(product)
    if len(thresh) > 0 and thresh[0] == product[0]:
        thresh = product[(np.abs(param_data[2]-param_data[0]) < 0.01) & np.abs((param_data[-2]-param_data[0]) < 0.01)]
    if len(thresh) > 0:
        ax.axvline(x=thresh[0], color='r', linewidth=1)
    
    if plot_index == 9:
        ax.set_xlabel(f"$\sigma_\\theta\sigma_\\rho$")
    else:
        ax.set_xticks([])

    print(f"{param_names[plot_index]}:\t mean:{np.mean((param_data[3]-param_data[0]) / true_sigmas)}\t"+
        f"95\% high: {np.mean((param_data[1]-param_data[0]) / true_sigmas)}\t 95\% low: {np.mean((param_data[-1]-param_data[0]) / true_sigmas)}\t"+
        f"68\% high: {np.mean((param_data[2]-param_data[0]) / true_sigmas)}\t 68\% low: {np.mean((param_data[-2]-param_data[0]) / true_sigmas)}")

custom_lines = [Line2D([0], [0], color='k', lw=4, alpha=0.3),
                Line2D([0], [0], color='k', lw=4, alpha=0.6),
                Line2D([0], [0], color='k', lw=1, linestyle='dashed'),
                Line2D([0], [0], color='k', lw=1, linestyle='dotted')]
fig.legend(custom_lines, ['95\%', '68\%', '50\%'], ncol=3, loc='upper center', prop={'size': LEGEND_SIZE}, bbox_to_anchor=(0.5,0.91))

plt.savefig("sigmas.pdf", bbox_inches="tight")
plt.savefig("sigmas.png", bbox_inches="tight")
plt.show()

print()
for i in range(N_DIM):
    param_data = np.zeros(len(true_sigmas) * N_PERCENTILES).reshape(N_PERCENTILES, len(true_sigmas))
    for f in percentiles.keys():
        param_data[:,name_index[f]] = percentiles[f][i]

    mean_95 = (np.mean((param_data[1]-param_data[0]) / true_sigmas) - np.mean((param_data[-1]-param_data[0]) / true_sigmas)) / 2
    mean_68 = (np.mean((param_data[2]-param_data[0]) / true_sigmas) - np.mean((param_data[-2]-param_data[0]) / true_sigmas)) / 2

    print(f"{param_names[i]}:\t mean:{np.mean((param_data[3]-param_data[0]) / true_sigmas)}\t95\%: {mean_95}\t 68\%: {mean_68}")
