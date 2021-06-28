import numpy as np
import matplotlib.pyplot as plt
import os

MARKER = "."

fig, ax1 = plt.subplots(figsize=(12, 4))
ax2 = ax1.twinx()

i = 0
for fname in os.listdir():
    if fname[:9] == "spin-mags":
        f = open(fname, 'r')
        xs = []
        ys = []
        zs = []
        os = []
        for line in f.readlines():
            if line == '':
                continue
            if line[0] == 'S':# Last line
                continue
            try:
                x, y, z, o = line.split(' ')
                os.append(float(o))
            except:
                x, y, z = line.split(' ')
            xs.append(float(x))
            ys.append(float(y))
            zs.append(float(z))

        f.close()

        ax1.plot(xs, linestyle='solid', color = "C"+str(i), linewidth=1)
        ax1.plot(ys, linestyle='dashed',color = "C"+str(i), linewidth=1)
        #ax1.plot(zs, linestyle='dotted', color = "C"+str(i), linewidth=1)
        ax1.plot([], [], label=fname, color = "C"+str(i), linewidth=1)

        if len(os) > 5:
            ax2.plot(os, color="C"+str(i), linestyle='dashdot')
        i+= 1

fig.tight_layout()
fig.legend(loc="lower left")
plt.show()
plt.savefig("img-spin-mags.png")
