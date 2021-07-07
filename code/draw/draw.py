import matplotlib.pyplot as plt
import numpy as np
import colorsys, sys
import matplotlib as mpl

FILENAME = sys.argv[1]
print("Drawing filename " + FILENAME)
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

def get_density_fraction(density):
    return (density - np.min(densities)) / (np.max(densities) - np.min(densities))
def get_color(f, b=0.5):
    return colorsys.hls_to_rgb((1-f)*0.7, b * 0.9, 1)

colors = []
for i in range(len(densities)):
    colors.append(get_color(get_density_fraction(densities[i]), brightnesses[i]))


for i, point in enumerate(points):
    tri = plt.Polygon(point, color=colors[i])
    plt.gca().add_patch(tri)

def get_continuous_cmap():
    cmap_colors=[]
    for i in range(256):
        cmap_colors.append(get_color(i/256))
    return mpl.colors.ListedColormap(cmap_colors, N=256)

plt.xlim(xlim[0]-0.1, xlim[1]+0.1)
plt.ylim(ylim[0]-0.1, ylim[1]+0.1)
plt.gca().set_aspect('equal', adjustable='box')
plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=min(densities), vmax=max(densities)),
    cmap=get_continuous_cmap()), label="Density")
plt.savefig(FILENAME[:-4] + ".pdf")
plt.show()
