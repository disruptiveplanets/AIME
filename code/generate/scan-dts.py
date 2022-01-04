import os
import numpy as np
import matplotlib.pyplot as plt

PERFECT_POW = -2
CALC_DATA = True
NUM_DTS = 50

def load_file(length):
    f = open("2-params-resolved.dat", 'r')
    perfect_xs = []
    perfect_ys = []
    perfect_zs = []
    for line in f.readlines():
        if line == "": continue
        x, y, z = line.split(" ")
        perfect_xs.append(float(x))
        perfect_ys.append(float(y))
        perfect_zs.append(float(z))
    f.close()

    if length is not None:
        while len(perfect_xs) < length:
            perfect_xs.append(perfect_xs[-1])
            perfect_ys.append(perfect_ys[-1])
            perfect_zs.append(perfect_zs[-1])
        while len(perfect_xs) > length:
            del perfect_xs[-1]
            del perfect_ys[-1]
            del perfect_zs[-1]

    return np.array([perfect_xs, perfect_ys, perfect_zs])

def get_results(dt, length=None):
    os.system(f'g++ -DASTEROIDS_MAX_J=0 -DASTEROIDS_MAX_K=2 -DDT={dt} -Wall -std=c++17 ../sim/*.cpp main.cpp -o "bin/generate" -O3')
    os.system("./bin/generate")
    return load_file(length)

def get_residuals():
    dts = 10**np.linspace(PERFECT_POW, 3, NUM_DTS)
    true_results = get_results(dts[0])
    dts = dts[1:]
    residuals = []
    for dt in dts:
        print("DT", dt)
        residuals.append(true_results - get_results(dt, length=len(true_results[0])))
    return np.array(residuals)

if __name__ == "__main__":
    if CALC_DATA:
        all_residuals = get_residuals()
        np.save("residuals.npy", all_residuals)

    
    dts = 10**np.linspace(PERFECT_POW, 3, NUM_DTS)[1:]

    with open("residuals.npy", 'rb') as f:
        residuals = np.load(f)

    colors = ["C0", "C1", "C2"]
    print(residuals.shape)
    print(np.percentile(residuals[:,0,:], 100, axis=1))
    print(np.percentile(residuals[:,1,:], 100, axis=1))
    print(np.percentile(residuals[:,2,:], 100, axis=1))
    
    for i in range(3):
        plt.plot(dts, np.percentile(residuals[:,i,:], 100, axis=1), color=colors[i])
        #plt.plot(dts, np.percentile(residuals[:,i,:], 90, axis=1), color=colors[i])
        #plt.plot(dts, np.percentile(residuals[:,i,:], 80, axis=1), color=colors[i])

    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("Delta t")
    plt.ylabel("Residual size")
    plt.savefig("dt-scan.png")
    plt.show()