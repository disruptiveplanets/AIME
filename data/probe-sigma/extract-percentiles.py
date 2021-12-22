# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
import os

plt.style.use("jcap")

file_names = os.listdir()
file_names.sort()

pf = open("percentiles.dat", 'w')

for bare in file_names:
    if not os.path.isdir(bare):
        continue
    for file in os.listdir(bare):
        if not file[-11:] == "samples.npy": continue
        array = np.loadtxt("{}/{}".format(bare, file))
        if len(array.shape) != 2: continue
        print(array.shape)
        data = []
        for i in range(array.shape[0]):
            data.append( ", ".join([
                str(np.mean(array[i])),
                str(np.percentile(array[i], 95.45)),
                str(np.percentile(array[i], 68.27)),
                str(np.percentile(array[i], 50)),
                str(np.percentile(array[i], 31.73)),
                str(np.percentile(array[i], 04.55)),
            ]))

        text = file + ": " +  ": ".join(data)+"\n"
        #print(text)
        pf.write(text)

pf.close()
