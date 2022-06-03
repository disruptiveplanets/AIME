import numpy as np
from matplotlib import pyplot as plt

NUM_DRAWS = 10
FILE_NAMES = "asym-ell-"

scatter_uncs = []
scatter_vols = []

for i in range(NUM_DRAWS):
    with open(FILE_NAMES + str(i) + "-fe.npy", 'rb') as f:
        long_samples = np.load(f)
        long_means = np.mean(long_samples, axis=0)
        high_unc = np.percentile(long_samples, (100 + 68.27) / 2, axis=0) - long_means
        low_unc = long_means - np.percentile(long_samples, (100 - 68.27) / 2, axis=0)
        unc_ratio = (high_unc + low_unc) / 2 / long_means
    print(long_means.shape)
    with open(FILE_NAMES + str(i) + "-grids.npy", 'rb') as f:
        grids = np.load(f)

    total_volume = 0
    volumes = []
    for u, g in zip(unc_ratio, grids):
        volume = np.sum(g)
        total_volume += volume
        scatter_uncs.append(u)
        volumes.append(volume)
    for v in volumes:
        scatter_vols.append(v / total_volume * len(volumes))

plt.figure()
plt.scatter(scatter_vols, scatter_uncs)
plt.xlabel("Volume [compared to 1/16 of total volume]")
plt.ylabel("Density uncertainty ratio")
plt.savefig("scatter.png")
plt.show()