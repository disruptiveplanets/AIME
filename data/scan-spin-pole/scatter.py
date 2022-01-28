# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import LinearNDInterpolator

LOW_ORDER_INDEX = 1
BASE_NAME = f"pole"
NUM_DIVISIONS = 64 # Run with 8 cores per process.
LOW_ORDER = [(0.39269908169, 0, -0.09766608), (0.39269908169, 0.05200629, -0.2021978)][LOW_ORDER_INDEX]
HIGH_ORDER = [0, 0, 0, 0, 0, 0, 0]

POLE_NORM = 2 * np.pi / (9*3600)

GOLDEN_RATIO = (1 + 5**0.5) / 2

def get_angle(index):
    return np.arccos(1 - 2*(index+0.5)/NUM_DIVISIONS), (2 * np.pi * index / GOLDEN_RATIO) % (2 * np.pi) - np.pi

def get_text(spin_x, spin_y, spin_z):
    return """0, 3
120
5
1000
6000
{}, {}, {}
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
0.78539816339, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
-0.78539816339, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5""".format(spin_x, spin_y, spin_z, LOW_ORDER[0], LOW_ORDER[1], LOW_ORDER[2], HIGH_ORDER[0], HIGH_ORDER[1],
    HIGH_ORDER[2], HIGH_ORDER[3], HIGH_ORDER[4], HIGH_ORDER[5], HIGH_ORDER[6])

pts = []
data = []
theta_phis = []
for pole_index in range(NUM_DIVISIONS):
    #f = open("../../staged/{}-{:02}.txt".format(BASE_NAME, pole_index), 'w')
    theta, phi = get_angle(pole_index)
    spin_x = POLE_NORM * np.sin(theta) * np.cos(phi)
    spin_y = POLE_NORM * np.sin(theta) * np.sin(phi)
    spin_z = POLE_NORM * np.cos(theta)
    d = np.cos(2*theta) * np.sin(2*phi)
    #pts.append([theta, phi])
    pts.append([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    theta_phis.append((np.pi / 2 - theta, phi))
    data.append(d)
    #f.write(get_text(spin_x, spin_y, spin_z))
    #f.close()
theta_phis = np.array(theta_phis)

lons = np.linspace(-np.pi, np.pi, 100)
lats = np.linspace(-np.pi/2, np.pi/2, 50)
Lon, Lat = np.meshgrid(lons, lats)
interp = LinearNDInterpolator(pts, data)

cart_array = []
for lat in lats:
    cart_line = []
    for lon in lons:
        shrink = 0.8
        x, y, z = np.cos(lat) * np.cos(lon) * shrink, np.cos(lat) * np.sin(lon) * shrink, np.sin(lat) * shrink
        cart_line.append(interp(x, y, z))
    cart_array.append(cart_line)

fig = plt.figure()
ax = fig.add_subplot(111, projection="mollweide")
im = ax.pcolormesh(Lon, Lat, cart_array)
ax.scatter(theta_phis[:,1], theta_phis[:,0], s=1, c='k')
fig.colorbar(im, orientation='horizontal')
plt.show()
