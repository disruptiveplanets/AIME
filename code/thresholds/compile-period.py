# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
import matplotlib as mpl

plt.style.use("jcap")
mpl.rcParams["font.size"] = 14

percentiles = {}
cadences = []
periods = []

def get_indices_from_name(name):
    return int(name[8:10]), int(name[11:13])

with open("cad-period-sync-contour.npy", 'rb') as f:
   uncs = np.load(f)[:,:,2] # Use median

# Get percentiles
with open("../../data/cad-period-sync-contour/percentiles.dat", 'r') as f:
    for line in f.readlines():
        if line == '': continue
        elements = line.split(':')
        name = elements[0]
        perc_array = []
        for percs in elements[1:]:
            perc_array.append([float(x) for x in percs.split(',')])
        perc_array = np.array(perc_array)
        percentiles[name] = perc_array

for name in percentiles.keys():
    dir_name = name[:13]
    with open(f"../../data/cad-period-sync-contour/{dir_name}/{dir_name}.txt", 'r') as f:
        max_j, max_l = f.readline().split(", ")
        max_j, max_l = (int(max_j), int(max_l))
        cadence = float(f.readline())
        perigee = int(f.readline())
        radius = float(f.readline())
        speed = float(f.readline())
        spin = np.array([float(x) for x in f.readline().split(',')])
        jlms = [float(x) for x in f.readline().split(',')]
        theta_true = [float(x) for x in f.readline().split(',')]
        theta_high = [float(x) for x in f.readline().split(',')]
        theta_low = [float(x) for x in f.readline().split(',')]
        sigma = [float(d) for d in f.readline().split(',')]
    cadence_index, period_index = get_indices_from_name(name)
    if period_index == 1:
        cadences.append(cadence / 60)
    if cadence_index == 1:
        period = 2 * np.pi / np.linalg.norm(spin) / 3600.0
        periods.append(period)

cadences = np.sort(cadences)
periods = np.sort(periods)
print(cadences),
print(periods)

fig = plt.figure(figsize=(6.6, 4))
p = plt.contourf(cadences, periods, uncs, levels=10)
cbar = plt.colorbar(p)
plt.ylabel(f"$P_\omega$ (hr)")
plt.xlabel(f"$\Delta t$ (min)")
plt.savefig("cad-period-contour-net.pdf", bbox_inches="tight")
plt.savefig("cad-period-contour-net.png", bbox_inches="tight")

plt.show()
