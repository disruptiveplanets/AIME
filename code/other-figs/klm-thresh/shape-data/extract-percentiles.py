# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
import os

DIST_DIMENSION = 1

plt.style.use("jcap")

file_names = os.listdir()
file_names.sort()

pf = open("percentiles.dat", 'w')

for bare in file_names:
    if not os.path.isdir(bare):
        continue

    min_dist = None
    min_text = None
    for file in os.listdir(bare):
        if not file[-11:] == "samples.npy": continue
        with open("{}/{}".format(bare, file), 'rb') as f:
            array = np.load(f)
        
        flat_samples = array.reshape(-1, array.shape[-1])

        # Get means
        means = np.mean(flat_samples, axis=0)
        data = []
        for i in range(flat_samples.shape[1]):
            data.append(", ".join([
                str(np.mean(flat_samples[:,i])),
                str(np.percentile(flat_samples[:,i], 50 + 95.45 / 2)),
                str(np.percentile(flat_samples[:,i], 50 + 68.27 / 2)),
                str(np.percentile(flat_samples[:,i], 50)),
                str(np.percentile(flat_samples[:,i], 50 - 68.27 / 2)),
                str(np.percentile(flat_samples[:,i], 50 - 95.45 / 2)),
            ]))

        min_text = file + ": " +  ": ".join(data)+"\n"

        #if min_text is None:
        #    continue
        pf.write(min_text)
        print(bare)

pf.close()
