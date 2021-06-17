import matplotlib.pyplot as plt
import numpy as np
import colorsys

FILENAME = "params.ast"
f = open(FILENAME, 'r')
points = []
densities = []
brightnesses = []

xlim = [0, 0]
ylim = [0, 0]

for line in f.readlines():
    if line == ' ':
        continue
    c1x, c1y, c2x, c2y, c3x, c3y, d, n = line.split(' ')
    points.append(np.asarray([
        [float(c1x), float(c1y)],
        [float(c2x), float(c2y)],
        [float(c3x), float(c3y)]]))

    xlim[0] = min(xlim[0], float(c1x), float(c2x), float(c3x))
    xlim[1] = max(xlim[1], float(c1x), float(c2x), float(c3x))
    ylim[0] = min(ylim[0], float(c1y), float(c2y), float(c3y))
    ylim[1] = max(ylim[1], float(c1y), float(c2y), float(c3y))
    brightnesses.append(float(n))
    densities.append(float(d))

colors = []
for i in range(len(densities)):
    colors.append(colorsys.hls_to_rgb(densities[i] / np.max(densities),
    brightnesses[i] * 0.9, 1))


for i, point in enumerate(points):
    tri = plt.Polygon(point, color=colors[i])
    plt.gca().add_patch(tri)

plt.xlim(xlim[0]-0.1, xlim[1]+0.1)
plt.ylim(ylim[0]-0.1, ylim[1]+0.1)
plt.gca().set_aspect('equal', adjustable='box')
plt.savefig(FILENAME[:-4] + ".pdf")
plt.show()
