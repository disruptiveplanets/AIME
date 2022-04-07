# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
import os
plt.style.use('jcap')

names = []
points = []
file_names = os.listdir()
file_names.sort()

DIFF_EXPECTED_FACTOR = [1e-2, 1e-5, 1e-4]
END_LENGTH = 100
FIG_SIZE = (5.333, 4)
KEY_INDICES = [24, 94, 101]


def get_ratios(k22, k20):
    Ixx = 2/3.0 * k20 - 4 * k22 + 2/3.0
    Iyy = 2/3.0 * k20 + 4 * k22 + 2/3.0
    Izz = -4/3.0 * k20 + 2/3.0
    asquared = 0.9 * (Iyy + Izz - Ixx)
    bsquared = 0.9 * (Ixx + Izz - Iyy) / asquared
    csquared = 0.9 * (Ixx + Iyy - Izz) / asquared
    a_per_c = 1 / np.sqrt(csquared)
    b_per_c = a_per_c * np.sqrt(bsquared)
    return a_per_c, b_per_c

def convert(thetas):
    out = []
    for entry in thetas:
        ac, bc = get_ratios(entry[1], entry[2])
        out.append([entry[0], ac, bc])
    return np.array(out)

def plot_pt():
    sym = get_ratios(0, -0.09766608)
    asym1 = get_ratios(0.05200629, -0.2021978)
    asym2 = get_ratios(-0.05200629, -0.2021978)
    plt.scatter([sym[0]], [sym[1]], color='k', marker='o')
    plt.scatter([asym1[0]], [asym1[1]], color='k', marker='s')
    plt.scatter([asym2[0]], [asym2[1]], color='k', marker='s')
# Fill points
for bare in file_names:
    if not os.path.isdir(bare):
        continue
    if 'compile' in bare:
        continue
    f = open("{0}/{0}.txt".format(bare), 'r')
    max_j, max_l = f.readline().split(", ")
    max_j, max_l = (int(max_j), int(max_l))
    cadence = int(f.readline())
    impact_parameter = int(f.readline())
    radius = float(f.readline())
    speed = float(f.readline())
    spin = [float(x) for x in f.readline().split(',')]
    jlms = [float(x) for x in f.readline().split(',')]
    theta_true = convert([[float(x) for x in f.readline().split(',')]])
    theta_high = [float(x) for x in f.readline().split(',')]
    theta_low = [float(x) for x in f.readline().split(',')]
    SIGMA = [float(d) for d in f.readline().split(',')]
    f.close()
    names.append(bare)
    points.append(theta_true)
points = np.array(points)

def covariance():
    side_length = int(-0.5 + np.sqrt(1 + 8 * len(names)) / 2)
    X = []
    Y = []
    corr12_data = []
    corr01_data = []
    corr02_data = []
    sigma0_data = []
    sigma1_data = []
    sigma2_data = []
    for i, k20 in enumerate(np.linspace(0, -0.25, side_length)):
        corr12_data_line = []
        corr01_data_line = []
        corr02_data_line = []
        sigma0_data_line = []
        sigma1_data_line = []
        sigma2_data_line = []
        X_line = []
        Y_line = []
        delta = 0.25 / side_length / 2 * (side_length - i - 1)
        for j, k22 in enumerate(np.linspace(-0.125+delta, 0.125+delta, side_length)):
            index = i * (i + 1) // 2 + j
            if i < j:
                cov = np.array([[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]])
            elif j == 0 or j == i or i == side_length - 1:
                cov = np.array([[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]])
            else:
                # Which data file is the given coords?
                data = [] # file, diffs, covs
                for file in os.listdir(names[index]):
                    if not file[-11:] == "samples.npy": continue
                    with open(f"{names[index]}/{file}", 'rb') as f:
                        #print(f"{names[index]}/{file}")
                        arrays = np.load(f).transpose()#[:,:,-END_LENGTH:]
                        if arrays.shape[1] == 0:
                            continue
                        arrays = arrays.reshape(arrays.shape[0], arrays.shape[1] * arrays.shape[2])
                        
                        arrays = convert(arrays.transpose()).transpose()

                        diffs = np.mean(arrays, axis=1) - points[index]
                        cov = np.cov(arrays)
                    data.append((file, diffs, cov))

                weights = [np.sum((data_row[1] / DIFF_EXPECTED_FACTOR) ** 2 ) for data_row in data]
                if len(weights) == 0:
                    cov = np.array([[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]])
                else:
                    min_index = np.argmin(weights)
                    _, _, cov = data[min_index]
            
            x, y = get_ratios(k22, k20)
            if np.isnan(x) or not np.isfinite(x):
                x = 0
            if np.isnan(y) or not np.isfinite(y):
                y = 0
                
            corr12 = cov[1][2] / np.sqrt(cov[1][1] * cov[2][2])
            corr01 = cov[0][1] / np.sqrt(cov[0][0] * cov[1][1])
            corr02 = cov[0][2] / np.sqrt(cov[0][0] * cov[2][2])
            corr12_data_line.append(corr12)
            corr01_data_line.append(corr01)
            corr02_data_line.append(corr02)
            sigma0_data_line.append(np.sqrt(cov[0][0]))
            sigma1_data_line.append(np.sqrt(cov[1][1]))
            sigma2_data_line.append(np.sqrt(cov[2][2]))

            if index in KEY_INDICES:
                print(f"{names[index]}: sigmas: {np.sqrt(cov[0][0])}, {np.sqrt(cov[1][1])}, {np.sqrt(cov[2][2])}. Corrs: {corr01}, {corr02}, {corr12}")

            X_line.append(x)
            Y_line.append(y)

        corr12_data.append(corr12_data_line)
        corr01_data.append(corr01_data_line)
        corr02_data.append(corr02_data_line)
        sigma0_data.append(sigma0_data_line)
        sigma1_data.append(sigma1_data_line)
        sigma2_data.append(sigma2_data_line)

        X_line = np.array(X_line)
        Y_line = np.array(Y_line)

        X.append(X_line)
        Y.append(Y_line)

    NUM_LEVELS = 21
    corr12_max = max(np.nanmax(corr12_data), -np.nanmin(corr12_data))
    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, corr12_data, levels=np.linspace(-corr12_max, corr12_max, NUM_LEVELS), cmap="BrBG_r")
    axc = plt.colorbar(c)
    axc.set_label("$\\textrm{Corr}(a/c, b/c)$")
    plot_pt()
    plt.xlabel("$a/c$")
    plt.ylabel("$b/c$")
    plt.xlim(1, 6.7)
    plt.ylim(1, 6.7)
    plt.tight_layout()
    plt.savefig("compile-figs/corrab.pdf")

    corr01_max = max(np.nanmax(corr01_data), -np.nanmin(corr01_data))
    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, corr01_data, levels=np.linspace(-corr01_max, corr01_max, NUM_LEVELS), cmap="BrBG_r")
    axc = plt.colorbar(c)
    axc.set_label("$\\textrm{Corr}(\\gamma_0 ,a/c)$")
    plot_pt()
    plt.xlabel("$a/c$")
    plt.ylabel("$b/c$")
    plt.xlim(1, 6.7)
    plt.ylim(1, 6.7)
    plt.tight_layout()
    plt.savefig("compile-figs/corr1a.pdf")

    corr02_max = max(np.nanmax(corr02_data), -np.nanmin(corr02_data))
    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, corr02_data, levels=np.linspace(-corr02_max, corr02_max, NUM_LEVELS), cmap="BrBG_r")
    axc = plt.colorbar(c)
    axc.set_label("$\\textrm{Corr}(\\gamma_0, b/c)$")
    plot_pt()
    plt.xlabel("$a/c$")
    plt.ylabel("$b/c$")
    plt.xlim(1, 6.7)
    plt.ylim(1, 6.7)
    plt.tight_layout()
    plt.savefig("compile-figs/corr1b.pdf")

    plt.figure(figsize=FIG_SIZE)
    flat_sig_0 = np.array(sigma0_data).reshape(-1) * 1e6
    c = plt.contourf(X, Y, np.array(sigma0_data) * 1e6, levels=np.linspace(0, np.max(np.sort(flat_sig_0[np.isfinite(flat_sig_0)])[:-7]), NUM_LEVELS), cmap='Purples_r')
    axc = plt.colorbar(c)
    axc.set_label("$\\sigma(\\gamma_0)$ ($\\times 10^{-6}$)")
    plot_pt()
    plt.xlabel("$a/c$")
    plt.ylabel("$b/c$")
    plt.xlim(1, 6.7)
    plt.ylim(1, 6.7)
    plt.tight_layout()
    plt.savefig("compile-figs/theta-1-ab-sigma.pdf")

    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, np.array(sigma1_data) * 10**6, levels=20, cmap="Purples_r")
    axc = plt.colorbar(c)
    axc.set_label("$\\sigma(a/c)$ ($\\times 10^{-6}$)")
    plot_pt()
    plt.xlabel("$a/c$")
    plt.ylabel("$b/c$")
    plt.xlim(1, 6.7)
    plt.ylim(1, 6.7)
    plt.tight_layout()
    plt.savefig("compile-figs/theta-a-sigma.pdf")

    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, np.array(sigma2_data) * 10**6, levels=20, cmap="Purples_r")
    axc = plt.colorbar(c)
    axc.set_label("$\\sigma(b/c)$ ($\\times 10^{-6}$)")
    plot_pt()
    plt.xlabel("$a/c$")
    plt.ylabel("$b/c$")
    plt.xlim(1, 6.7)
    plt.ylim(1, 6.7)
    plt.tight_layout()
    plt.savefig("compile-figs/theta-b-sigma.pdf")


if __name__ == "__main__":
    covariance()
    plt.show()
