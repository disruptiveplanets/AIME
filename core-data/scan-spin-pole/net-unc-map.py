from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.linalg import norm
import os

plt.style.use("jcap")

xyzs = []
theta_phis = []

PULL = False

if PULL:
    os.system("scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/thresholds/lumpy/scan-spin-pole.npy .")

with open("scan-spin-pole.npy", 'rb') as f:
    uncs = np.load(f)[:,2] # Use medians

def show_true_point(ax, flip=False):
    spin = [0.00006464182, 0.00012928364, -0.00012928364]
    theta, phi = np.pi / 2 - np.arccos(spin[2] / norm(spin)), np.arctan2(spin[1], spin[0])
    if flip:
        theta *= -1
    ax.scatter(phi, theta, color='tab:orange', marker='*', s=32)

for i in range(len(uncs)):
    dir_name = "pole-{:02}".format(i)
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
    sigma_theta = sigma[0]
    theta, phi = np.pi / 2 - np.arccos(spin[2] / norm(spin)), np.arctan2(spin[1], spin[0])
    theta_phis.append((theta, phi))
    xyzs.append(np.array(spin) / norm(spin))

theta_phis = np.array(theta_phis)
lons = np.linspace(-np.pi, np.pi, 60)
lats = np.linspace(-np.pi/2, np.pi/2, 30)
Lon, Lat = np.meshgrid(lons, lats)

fig, ax = plt.subplots(subplot_kw=dict(projection="mollweide"))

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
cart_array_unc = np.array(cart_array_unc) * 10**4

levels = np.linspace(0, 1.2, 13)
im = ax.contourf(Lon, Lat, cart_array_unc, cmap='Blues_r', levels=levels)

ax.contour(Lon, Lat, cart_array_unc, levels=[1], colors=['r'], linewidths=[1], linestyles=['dashed'])
    
cbar = fig.colorbar(im)
cbar.set_label("$\sigma_\\rho / \\rho$ ($\\times 10^{-4}$)")

show_true_point(ax)
ax.set_xticks([-3 * np.pi / 4, -np.pi / 2, -np.pi / 4, 0, np.pi / 4, np.pi / 2, 3 * np.pi / 4])
ax.set_yticks([np.pi/3, np.pi/6, 0, -np.pi/6, -np.pi/3])
ax.set_xticklabels(['']*7)
ax.set_yticklabels(['']*5)
ax.grid(True)

fig.tight_layout()

plt.savefig("unc-map.pdf")
plt.savefig("unc-map.png")
plt.show()