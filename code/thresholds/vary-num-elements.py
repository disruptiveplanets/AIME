import os
import numpy as np
import matplotlib.pyplot as plt

PULL = True
DIRECTORY = "fe"

if PULL:
    for i in range(9):
        os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/thresholds/fe/probe-s-theta-{i}.npy {DIRECTORY}/")

plt.figure()
plt.title("Spin pole uncertainty")
plt.xlabel("sigma theta")
plt.ylabel("Uncertainty ratio")

for i in range(9):
    fname = f"probe-s-theta-{i}.npy"
    if fname.endswith(".npy"):
        with open(DIRECTORY + "/" + fname, 'rb') as f:
            uncs = np.load(f)
        if len(fname[:-4]) < 8:
            print(fname[:-4], end='\t\t')
        else:
            print(fname[:-4], end='\t')
        
        plt.plot(uncs, label=fname[:-4])
plt.savefig("vary-num.png")
plt.show()