# Goal: check that the usual deviation of delta / sigma is equal to one. It is.
from matplotlib import pyplot as plt
import numpy as np

plt.style.use("jcap")

percentiles = {}

# Get percentiles
with open("percentiles.dat", 'r') as f:
    for line in f.readlines():
        if line == '': continue
        elements = line.split(':')
        name = elements[0]
        perc_array = []
        for percs in elements[1:]:
            perc_array.append([float(x) for x in percs.split(',')])
        perc_array = np.array(perc_array)
        N_DIM = perc_array.shape[0]
        N_PERCENTILES = perc_array.shape[1]
        percentiles[name] = perc_array

trues = np.array([0.39269908169, 0.05200629, -0.2021978, 0, 0, 0, 0, 0, 0, 0])

ratios = []
for name, p in percentiles.items():
    delta = p[:,0] - trues
    sigma = (np.abs(p[:,0] - p[:,2]) + np.abs(p[:,0] - p[:,-2])) / 2
    ratios.append(delta / sigma)

print(np.mean(ratios, axis=0))
print(np.std(ratios, axis=0))

plt.plot(ratios)
plt.figure()
for r in np.array(ratios).transpose():
    plt.hist(r, histtype='step')

plt.figure()
plt.hist(np.array(ratios).reshape(-1), histtype='step')

plt.show()