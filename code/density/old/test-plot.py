import numpy as np
import matplotlib.pyplot as plt

VERY_SMALL = 0.01
SKIP = 7
EXPAND_X=1.15

fig = plt.figure(figsize=(7, 7))
ax = fig.gca(projection='3d')
#ax.set_axis_off()
ax.grid(False)

line = np.linspace(-1, 1, 50)
levels = np.linspace(0, 1, 40)

densities = np.full((len(line), len(line), len(line)), np.nan)
for nx, x in enumerate(line):
    for ny, y in enumerate(line):
        for nz, z in enumerate(line):
            if x*x*2+y*y+z*z<1:
                densities[nx,ny,nz] = x*x+y*y+z*z

densities = np.array(densities)

for i in range(0, len(line), SKIP):
    z = line[i]
    ax.contourf(line, line, z+densities[:,:,i]*VERY_SMALL, zdir='z', levels=z+VERY_SMALL*levels)

ax.view_init(elev=15., azim=45)

ax.set_xlim3d(-EXPAND_X, EXPAND_X)
ax.set_ylim3d(-EXPAND_X, EXPAND_X)
ax.set_zlim3d(-1, 1)

ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

ax.set_xlabel("$x$")
ax.set_ylabel("$y$")
ax.set_zlabel("$z$")

plt.show()