import matplotlib.pyplot as plt
import numpy as np

dots = []

start = False

f = open("nutation.txt", 'r')
for line in f.readlines():
    if line == '' or line[0] == "S" or line[0] == 'M':
        continue
    if line == '\n':
        start = True
        continue
    if start:
        dots.append(float(line))

plt.xlabel("Time step")
plt.ylabel("$\\L \\cdot \\omega")
plot_x = np.arange(len(dots))
plt.plot(plot_x, dots, label="x")

plt.savefig("img-nutation.png")
plt.show()
