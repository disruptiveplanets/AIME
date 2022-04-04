import numpy as np

params = "gamma K22 K20 ReK33 ImK33 ReK32 ImK32 ReK31 ImK31 K30".split()
percentiles = {}

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

def show_ratios(name):
    print(name)
    name_sigma = (percentiles[name][:, 2] - percentiles[name][:, -2]) / 2
    earth = (percentiles["earth-1-samples.npy"][:, 2] - percentiles["earth-1-samples.npy"][:, -2]) / 2
    ratio = name_sigma / earth
    for p, r in zip(params, ratio):
        print(f"{p}\t{r}")

show_ratios("jupiter-same-orbit-0-samples.npy")