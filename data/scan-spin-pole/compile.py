# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.interpolate import LinearNDInterpolator
from scipy.linalg import norm

plt.style.use("jcap")

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

percentiles = {}
name_index = {}
true_sigma = None
xyzs = []
theta_phis = []

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

index = 0
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
    name_index[name] = index
    sigma_theta = sigma[0]
    index += 1
    theta, phi = np.pi / 2 - np.arccos(spin[2] / norm(spin)), np.arctan2(spin[1], spin[0])
    theta_phis.append((theta, phi))
    xyzs.append(np.array(spin) / norm(spin))

theta_phis = np.array(theta_phis)
lons = np.linspace(-np.pi, np.pi, 180)
lats = np.linspace(-np.pi/2, np.pi/2, 90)
Lon, Lat = np.meshgrid(lons, lats)

fig, axs = plt.subplots(figsize=(14, 8), ncols=3, nrows=4, subplot_kw=dict(projection="mollweide"))
axs = axs.reshape(-1)

data = []
net_data = np.zeros(len(percentiles))
for i in range(N_DIM):
    data_line = []
    for f in percentiles.keys():
        scale = 1 if i >= 3 else 1e5;
        data_line.append(np.abs(percentiles[f][i][2] - percentiles[f][i][0]) * scale / sigma_theta)
    data.append(data_line)
    net_data += np.array(data_line) / np.mean(data_line) / N_DIM

net_data *= np.mean(data)

i = 0
for plot_index in range(N_DIM+1):
    if plot_index == 9:
        continue
    param_data = np.zeros(len(theta_phis) * N_PERCENTILES).reshape(N_PERCENTILES, len(theta_phis))
    for f in percentiles.keys():
        param_data[:,name_index[f]] = percentiles[f][i]

    interp = LinearNDInterpolator(xyzs, data[i])

    cart_array = []
    for lat in lats:
        cart_line = []
        for lon in lons:
            shrink = 0.8
            x, y, z = np.cos(lat) * np.cos(lon) * shrink, np.cos(lat) * np.sin(lon) * shrink, np.sin(lat) * shrink
            cart_line.append(interp(x, y, z))
        cart_array.append(cart_line)

    im = axs[plot_index].pcolormesh(Lon, Lat, cart_array)#, vmax=np.max(data), vmin=np.min(data))
    cbar = fig.colorbar(im, ax=axs[plot_index])
    if i < 3:
        cbar.set_label(f"$\sigma({param_names[i]}) / \sigma_\\theta$ ($\\times 10^{{-5}})$")
    else:
        cbar.set_label(f"$\sigma({param_names[i]}) / \sigma_\\theta$")

    if plot_index != 10:
        axs[plot_index].set_xticks([])
        axs[plot_index].set_yticks([])
    else:
        axs[plot_index].scatter(theta_phis[:,1], theta_phis[:,0], s=1, c='k')
        axs[plot_index].set_xlabel("$\\theta$")
        axs[plot_index].set_ylabel("$\\phi$")

    i += 1

axs[9].remove()
axs[11].remove()

custom_lines = [Line2D([0], [0], color='k', lw=4, alpha=0.3),
                Line2D([0], [0], color='k', lw=4, alpha=0.6),
                Line2D([0], [0], color='k', lw=1, linestyle='dashed')]
#fig.colorbar(im, orientation='horizontal')
fig.tight_layout()

line = plt.Line2D([0,1],[0.75, 0.75], transform=fig.transFigure, color="black")
fig.add_artist(line)

plt.savefig("pole.pdf")
plt.savefig("pole.png")



# Plot average
interp = LinearNDInterpolator(xyzs, net_data)
fig, ax = plt.subplots(ncols=1, nrows=1, subplot_kw=dict(projection="mollweide"))
cart_array = []
for lat in lats:
    cart_line = []
    for lon in lons:
        shrink = 0.8
        x, y, z = np.cos(lat) * np.cos(lon) * shrink, np.cos(lat) * np.sin(lon) * shrink, np.sin(lat) * shrink
        cart_line.append(interp(x, y, z))
    cart_array.append(np.array(cart_line))
print(cart_array)
im = ax.pcolormesh(Lon, Lat, cart_array)#, vmax=np.max(data), vmin=np.min(data))
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(f"$\overline{{\sigma}} / \sigma_\\theta$")
fig.tight_layout()

plt.savefig("avg-pole.pdf")
plt.savefig("avg-pole.png")


plt.show()
