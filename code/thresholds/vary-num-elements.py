import os
import numpy as np
import matplotlib.pyplot as plt

PULL = False
DIRECTORY = "fe"

if PULL:
    for i in range(9):
        os.system(f"scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/thresholds/fe/probe-s-rho-{i}.npy {DIRECTORY}/")

plt.figure()
plt.title("Spin pole uncertainty")
plt.xlabel("sigma theta")
plt.ylabel("Uncertainty ratio")

for i in range(9):
    fname = f"probe-s-rho-{i}.npy"
    if not os.path.exists(DIRECTORY + "/" + fname):
        continue
    if fname.endswith(".npy"):
        with open(DIRECTORY + "/" + fname, 'rb') as f:
            uncs = np.load(f)
        if len(fname[:-4]) < 8:
            print(fname[:-4], end='\t\t')
        else:
            print(fname[:-4], end='\t')
        plt.scatter(np.arange(len(uncs)), uncs, label=f"{9-i} DOF")
        plt.plot(uncs)
plt.legend()
plt.savefig("vary-num.png")
plt.show()