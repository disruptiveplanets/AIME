# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay
from scipy.linalg import norm

plt.style.use("jcap")

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

name_index = {}
true_sigma = None
xyzs = []
theta_phis = []

AXIS_SIZE = 12
LEGEND_SIZE = 12

N_DIM = None
N_PERCENTILES = None

# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
import os

TRUE_THETA = np.array([0.39269908169, 0.05200629, -0.2021978, 0, 0, 0, 0, 0, 0, 0])
DIST_DIMENSION = 1

plt.style.use("jcap")

file_names = os.listdir()
file_names.sort()
index = 0
for dir_name in file_names:
    if not os.path.isdir(dir_name):
        continue
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
    name_index[dir_name+".txt"] = index
    sigma_theta = sigma[0]
    index += 1
    theta, phi = np.pi / 2 - np.arccos(spin[2] / norm(spin)), np.arctan2(spin[1], spin[0])
    theta_phis.append((theta, phi))
    xyzs.append(np.array(spin) / norm(spin))
xyzs = np.array(xyzs)
tris = Delaunay(xyzs)

data = []
for bare in file_names:
    if not os.path.isdir(bare):
        continue
    min_dist = None
    min_trace = None
    for file in os.listdir(bare):
        if not file[-11:] == "samples.npy": continue
        with open("{}/{}".format(bare, file), 'rb') as f:
            array = np.load(f)
        
        flat_samples = array.reshape(-1, array.shape[-1])

        # Get means
        means = np.mean(flat_samples, axis=0)
        dist = (means - TRUE_THETA)[DIST_DIMENSION]
        if min_dist is None or dist < min_dist:
            min_trace=np.trace(np.cov(flat_samples.transpose()))
            min_dist = dist
    data.append(min_trace)
data = np.array(data)

theta_phis = np.array(theta_phis)
lons = np.linspace(-np.pi, np.pi, 180)
lats = np.linspace(-np.pi/2, np.pi/2, 90)
Lon, Lat = np.meshgrid(lons, lats)
data = np.sqrt(np.array(data) / sigma_theta**2)

# Plot average
with open("projecteds.dat", 'rb') as f:
    projected_vecs = np.load(f)
interp = LinearNDInterpolator(tris, data)
cart_array = []
for i, lat in enumerate(lats):
    cart_line = []
    for j, lon in enumerate(lons):
        p = projected_vecs[i][j]
        cart_line.append(interp(p[0], p[1], p[2]))
    cart_array.append(np.array(cart_line))

fig, ax = plt.subplots(figsize=(6,5), ncols=1, nrows=1, subplot_kw=dict(projection="mollweide"))
im = ax.pcolormesh(Lon, Lat, cart_array, vmin=0, cmap='Reds_r')#, vmax=np.max(data), vmin=np.min(data))
cbar = fig.colorbar(im, ax=ax, orientation='horizontal')
cbar.set_label(f"$\sqrt{{\mathrm{{tr}}(\Sigma) / \sigma_\\theta^2}}$")
ax.scatter(theta_phis[:,1], theta_phis[:,0], s=1, c='k')
ax.set_xlabel("$\\theta$")
ax.set_ylabel("$\\phi$")
ax.set_xticks([-3 * np.pi / 4, -np.pi / 2, -np.pi / 4, 0, np.pi / 4, np.pi / 2, 3 * np.pi / 4])
ax.set_yticks([np.pi/3, np.pi/6, 0, -np.pi/6, -np.pi/3])
ax.grid(True)
fig.tight_layout()
plt.savefig("trace-pole-mollweidee.pdf")
plt.savefig("trace-pole-mollweide.png")

fig, ax = plt.subplots(figsize=(6,5), ncols=1, nrows=1, subplot_kw=dict(projection="polar"))
im = ax.pcolormesh(Lon, Lat, cart_array, vmin=0, cmap='Reds_r')#, vmax=np.max(data), vmin=np.min(data))
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(f"$\sqrt{{\mathrm{{tr}}(\Sigma) / \sigma_\\theta^2}}$")
fig.tight_layout()
ax.scatter(theta_phis[:,1], theta_phis[:,0], s=1, c='k')
ax.set_yticklabels(['']*5)
ax.grid(True)
plt.savefig("trace-pole-polar.pdf")
plt.savefig("trace-pole-polar.png")


plt.show()
