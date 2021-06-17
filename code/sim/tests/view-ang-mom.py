import matplotlib.pyplot as plt
import numpy as np

lx = []
ly = []
lz = []

start = False

f = open("ang-mom.txt", 'r')
for line in f.readlines():
    if line == '' or line[0] == "S":
        continue
    if line == '\n':
        start = True
        continue
    if start:
        x, y, z = line[1:-2].split(' ')
        lx.append(float(x))
        ly.append(float(y))
        lz.append(float(z))

plt.xlabel("Time step")
plt.ylabel("Angular moment components")
plot_x = np.arange(len(lx))
plt.plot(plot_x, lx, label="x")
plt.plot(plot_x, ly, label="y")
plt.plot(plot_x, lz, label="z")
plt.legend()

plt.savefig("img-ang-mom.png")
plt.show()
