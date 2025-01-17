# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
import matplotlib as mpl
from scipy.interpolate import LinearNDInterpolator
from scipy.linalg import norm
import os

plt.style.use("jcap")
mpl.rcParams["font.size"] = 14

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

percentiles = {}
name_index = {}
xyzs = []
theta_phis = []

N_DIM = None
N_PERCENTILES = None
PULL = False
DRAW_AVG_PLOTS = False

THRESHOLD = 0.001

if PULL:
    os.system("scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/thresholds/lumpy/scan-spin-pole.npy .")

with open("scan-spin-pole.npy", 'rb') as f:
    uncs = np.load(f)[:,2]

def show_true_point(ax, flip=False):
    spin = [0.00006464182, 0.00012928364, -0.00012928364]
    theta, phi = np.pi / 2 - np.arccos(spin[2] / norm(spin)), np.arctan2(spin[1], spin[0])
    if flip:
        theta *= -1
    ax.scatter(phi, theta, color='tab:orange', marker='*', s=32)

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
lons = np.linspace(-np.pi, np.pi, 60)
lats = np.linspace(-np.pi/2, np.pi/2, 30)
Lon, Lat = np.meshgrid(lons, lats)

fig, axs = plt.subplots(figsize=(14, 8), ncols=3, nrows=4, subplot_kw=dict(projection="mollweide"))
axs = axs.reshape(-1)

data = []
net_data = np.zeros(len(percentiles))
for i in range(N_DIM):
    data_line = []
    for f in percentiles.keys():
        scale = 1e3 if i >= 3 else 1e7;
        data_line.append(np.abs(percentiles[f][i][2] - percentiles[f][i][0]) * scale)
    data.append(data_line)
    net_data += np.array(data_line) / np.mean(data_line) / N_DIM

with open("projecteds.dat", 'rb') as f:
    projected_vecs = np.load(f)

interp_unc = LinearNDInterpolator(xyzs, uncs)

cart_array_unc = []
for pi, lat in enumerate(lats):
    cart_line_unc = []
    for pj, lon in enumerate(lons):
        p = projected_vecs[pi][pj]
        cart_line_unc.append(interp_unc(p[0], p[1], p[2]))
    cart_array_unc.append(np.array(cart_line_unc))
cart_array_unc = np.array(cart_array_unc)

i = 0
for plot_index in range(N_DIM+1):
    if plot_index == 9:
        continue
    param_data = np.zeros(len(theta_phis) * N_PERCENTILES).reshape(N_PERCENTILES, len(theta_phis))
    for f in percentiles.keys():
        param_data[:,name_index[f]] = percentiles[f][i]

    interp = LinearNDInterpolator(xyzs, data[i])
    interp_unc = LinearNDInterpolator(xyzs, uncs)

    cart_array = []
    for pi, lat in enumerate(lats):
        cart_line = []
        for pj, lon in enumerate(lons):
            p = projected_vecs[pi][pj]
            cart_line.append(interp(p[0], p[1], p[2]))
        cart_array.append(cart_line)

    im = axs[plot_index].pcolormesh(Lon, Lat, cart_array, cmap='Blues_r', shading="gouraud")

    axs[plot_index].contour(Lon, Lat, cart_array_unc, levels=[THRESHOLD], colors=['r'], linewidths=[1])
    # axs[plot_index].contour(Lon, Lat, cart_array_unc, levels=[1e-4], colors=['r'], linewidths=[1], linestyles=['dashed'])
    
    cbar = fig.colorbar(im, ax=axs[plot_index], extend='max')
    if i < 3:
        cbar.set_label(f"$\sigma({param_names[i]})$ ($\\times 10^{{-7}})$")
    else:
        cbar.set_label(f"$\sigma({param_names[i]})$ ($\\times 10^{{-3}})$")

    if plot_index == 10:
        axs[plot_index].scatter(theta_phis[:,1], theta_phis[:,0], s=1, c='limegreen')

    show_true_point(axs[plot_index])
    axs[plot_index].set_xticks([-3 * np.pi / 4, -np.pi / 2, -np.pi / 4, 0, np.pi / 4, np.pi / 2, 3 * np.pi / 4])
    axs[plot_index].set_yticks([np.pi/3, np.pi/6, 0, -np.pi/6, -np.pi/3])
    axs[plot_index].set_xticklabels(['']*7)
    axs[plot_index].set_yticklabels(['']*5)
    axs[plot_index].grid(True)
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

if DRAW_AVG_PLOTS:

    # Plot average    
    interp = LinearNDInterpolator(xyzs, net_data)
    interp_unc = LinearNDInterpolator(xyzs, uncs)
    cart_array = []
    for i, lat in enumerate(lats):
        cart_line = []
        for j, lon in enumerate(lons):
            p = projected_vecs[i][j]
            cart_line.append(interp(p[0], p[1], p[2]))
        cart_array.append(np.array(cart_line))

    fig, ax = plt.subplots(figsize=(6,5), ncols=1, nrows=1, subplot_kw=dict(projection="mollweide"))
    im = ax.pcolormesh(Lon, Lat, cart_array, cmap='Blues_r', shading="gouraud")
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal')
    cbar.set_label(f"$\overline{{\sigma}}$")
    ax.scatter(theta_phis[:,1], theta_phis[:,0], s=1, c='k')
    show_true_point(ax)
    ax.set_xlabel("$\\theta$")
    ax.set_ylabel("$\\phi$")
    ax.set_xticks([-3 * np.pi / 4, -np.pi / 2, -np.pi / 4, 0, np.pi / 4, np.pi / 2, 3 * np.pi / 4])
    ax.set_yticks([np.pi/3, np.pi/6, 0, -np.pi/6, -np.pi/3])
    ax.contour(Lon, Lat, np.array(cart_array_unc), levels=[1.0], colors=['r'], linewidths=[1])
    ax.contour(Lon, Lat, np.array(cart_array_unc), levels=[0.2], colors=['r'], linewidths=[1], linestyles=['dotted'])
    ax.grid(True)
    fig.tight_layout()
    plt.savefig("avg-pole-mollweide.pdf")
    plt.savefig("avg-pole-mollweide.png")

    fig, ax = plt.subplots(figsize=(6,5), ncols=1, nrows=1, subplot_kw=dict(projection="polar"))
    im = ax.pcolormesh(Lon, -Lat, cart_array, cmap='Blues_r', shading="gouraud")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(f"$\overline{{\sigma}}$")
    fig.tight_layout()
    ax.scatter(theta_phis[:,1], -theta_phis[:,0], s=1, c='k')
    ax.contour(Lon, -Lat, cart_array_unc, levels=[1.0], colors=['r'], linewidths=[1])
    ax.contour(Lon, -Lat, cart_array_unc, levels=[0.2], colors=['r'], linewidths=[1], linestyles=['dotted'])
    show_true_point(ax, flip=True)
    ax.set_yticklabels(['']*5)
    ax.grid(True)
    plt.savefig("avg-pole-polar.pdf")
    plt.savefig("avg-pole-polar.png")


plt.show()