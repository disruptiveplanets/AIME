# Goal: spawn a bunch of runs that cover the parameter space.
from re import S
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
KEY_INDICES = [24, 94, 101]


FIG_SIZE = (5.333, 4)

def plot_pt():
    plt.scatter([0], [-0.09766608], color='r', marker='o')
    plt.scatter([0.05200629], [-0.2021978], color='k', marker='^')
    plt.scatter([-0.05200629], [-0.2021978], color='k', marker='^')

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
    perigee = int(f.readline())
    radius = float(f.readline())
    speed = float(f.readline())
    spin = [float(x) for x in f.readline().split(',')]
    jlms = [float(x) for x in f.readline().split(',')]
    theta_true = [float(x) for x in f.readline().split(',')]
    theta_high = [float(x) for x in f.readline().split(',')]
    theta_low = [float(x) for x in f.readline().split(',')]
    SIGMA = [float(d) for d in f.readline().split(',')]
    f.close()
    names.append(bare)
    points.append(theta_true)
points = np.array(points)

'''plt.figure()
for i, point in enumerate(points):
    success_count = 0
    for file in os.listdir(names[i]):
        if file[-11:] == "samples.npy":
            success_count += 1
    plt.text(x=point[1], y=point[2], s=str(success_count))
plt.title("Distribution of parameter points")
plt.xlim(-0.15, 0.15)
plt.ylim(-0.28, 0.03)
plt.show()'''

# Plot pdfs
def plot_pdf(index):
    side_length = int(-0.5 + np.sqrt(1 + 8 * len(names)) / 2)
    fig = plt.figure(constrained_layout=True, figsize=(10, 7))
    gs = fig.add_gridspec(side_length, 2 * side_length)
    ax_min = 0
    ax_max = 0
    axs = []
    for i in range(side_length):
        for j in range(side_length):
            if i < j:
                continue
            ax = fig.add_subplot(gs[i, (side_length-i)+2*j-1:(side_length-i)+2*j+1])
            ax.set_axis_off()
            #ax.set_yticklabels([])
            axs.append(ax)

            name_index = j + i * (i+1) // 2

            data = []# file, means
            for file in os.listdir(names[name_index]):
                if not file[-11:] == "samples.npy": continue
                array = np.loadtxt("{}/{}".format(names[name_index], file))[index]
                means = np.mean(array)
                data.append((file, means))

            ### I took this apart
            sys.exit()

            # plot the minimizing array
            if min_dist is None or min_dist > DIFF_THRESHOLD[index]:
                continue
            array = np.loadtxt("{}/{}".format(names[name_index], min_file))[index] - points[name_index][index]
            _, bins, _ = ax.hist(array, bins=12, fill=False, density=True)
            ax_min = min(ax_min, np.min(bins))
            ax_max = max(ax_max, np.max(bins))
            ax.text(x=0.05, y=0.9, s="{}".format(name_index), transform=ax.transAxes)
            ax.axvline(x=0, color='b')

    for i, ax in enumerate(axs):
        ax.set_xlim(ax_min, ax_max)
        #if i < len(names) - side_length:
        #    ax.set_xticklabels([])

    fig.savefig("theta-{}-hists.pdf".format(index))
    fig.savefig("theta-{}-hists.png".format(index))

def covariance():
    side_length = int(-0.5 + np.sqrt(1 + 8 * len(names)) / 2)
    print(side_length)
    X = []
    Y = []
    cov12_data = []
    corr12_data = []
    cov01_data = []
    corr01_data = []
    cov02_data = []
    corr02_data = []
    sigma0_data = []
    sigma1_data = []
    sigma2_data = []
    for i, y in enumerate(np.linspace(0, -0.25, side_length)):
        cov12_data_line = []
        corr12_data_line = []
        cov01_data_line = []
        corr01_data_line = []
        cov02_data_line = []
        corr02_data_line = []
        sigma0_data_line = []
        sigma1_data_line = []
        sigma2_data_line = []
        X_line = []
        Y_line = []
        delta = 0.25 / side_length / 2 * (side_length - i - 1)
        for j, x in enumerate(np.linspace(-0.125+delta, 0.125+delta, side_length)):
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

                        diffs = np.mean(arrays, axis=1) - points[index]
                        cov = np.cov(arrays)
                    data.append((file, diffs, cov))

                weights = [np.sum((data_row[1] / DIFF_EXPECTED_FACTOR) ** 2 ) for data_row in data]
                if len(weights) == 0:
                    cov = np.array([[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]])
                else:
                    min_index = np.argmin(weights)
                    min_file, min_diffs, cov = data[min_index]
                    print(f"Min file: {min_file}; Min diffs: {min_diffs}")
                
            corr12 = cov[1][2] / np.sqrt(cov[1][1] * cov[2][2])
            corr01 = cov[0][1] / np.sqrt(cov[0][0] * cov[1][1])
            corr02 = cov[0][2] / np.sqrt(cov[0][0] * cov[2][2])
            cov12_data_line.append(cov[1][2])
            corr12_data_line.append(corr12)
            cov01_data_line.append(cov[0][1])
            corr01_data_line.append(corr01)
            cov02_data_line.append(cov[0][2])
            corr02_data_line.append(corr02)
            sigma0_data_line.append(np.sqrt(cov[0][0]))
            sigma1_data_line.append(np.sqrt(cov[1][1]))
            sigma2_data_line.append(np.sqrt(cov[2][2]))

            if index in KEY_INDICES:
                print(f"{names[index]}: sigmas: {np.sqrt(cov[0][0])}, {np.sqrt(cov[1][1])}, {np.sqrt(cov[2][2])}. Corrs: {corr01}, {corr02}, {corr12}")

            #if not np.isnan(corr):
            #    print(np.sqrt(cov[0][0]), names[index])

            X_line.append(x)
            Y_line.append(y)
        cov12_data.append(cov12_data_line)
        corr12_data.append(corr12_data_line)
        cov01_data.append(cov01_data_line)
        corr01_data.append(corr01_data_line)
        cov02_data.append(cov02_data_line)
        corr02_data.append(corr02_data_line)
        sigma0_data.append(sigma0_data_line)
        sigma1_data.append(sigma1_data_line)
        sigma2_data.append(sigma2_data_line)
        X.append(X_line)
        Y.append(Y_line)


    corr12_max = max(np.nanmax(corr12_data), -np.nanmin(corr12_data))
    NUM_LEVELS = 21
    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, corr12_data, levels=np.linspace(-corr12_max, corr12_max, NUM_LEVELS), cmap='BrBG_r')
    axc = plt.colorbar(c)
    axc.set_label("$\\textrm{Corr}(K_{22}, K_{20})$")
    plot_pt()
    plt.xlim(-0.11, 0.11)
    plt.ylim(-0.24, -0.03)
    plt.xlabel("$K_{22}$")
    plt.ylabel("$K_{20}$")
    plt.tight_layout()
    plt.savefig("compile-figs/corr23.pdf")
    plt.savefig("compile-figs/corr23.png")

    corr01_max = max(np.nanmax(corr01_data), -np.nanmin(corr01_data))
    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, np.array(corr01_data), levels=np.linspace(-corr01_max, corr01_max, NUM_LEVELS), cmap='BrBG_r')
    axc = plt.colorbar(c)
    axc.set_label("$\\textrm{Corr}(\\gamma_0, K_{22})$")
    plot_pt()
    plt.xlim(-0.11, 0.11)
    plt.ylim(-0.24, -0.03)
    plt.xlabel("$K_{22}$")
    plt.ylabel("$K_{20}$")
    plt.tight_layout()
    plt.savefig("compile-figs/corr12.pdf")
    plt.savefig("compile-figs/corr12.png")

    corr02_max = max(np.nanmax(corr02_data), -np.nanmin(corr02_data))
    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, np.array(corr02_data), levels=np.linspace(-corr02_max, corr02_max, NUM_LEVELS), cmap='BrBG_r')
    axc = plt.colorbar(c)
    axc.set_label("$\\textrm{Corr}(\\gamma_0, K_{20})$")
    plot_pt()
    plt.xlim(-0.11, 0.11)
    plt.ylim(-0.24, -0.03)
    plt.xlabel("$K_{22}$")
    plt.ylabel("$K_{20}$")
    plt.tight_layout()
    plt.savefig("compile-figs/corr13.pdf")
    plt.savefig("compile-figs/corr13.png")

    plt.figure(figsize=FIG_SIZE)
    flat_sig_0 = np.array(sigma0_data).reshape(-1) * 1e6
    c = plt.contourf(X, Y, np.array(sigma0_data) * 1e6, levels=np.linspace(0, np.max(np.sort(flat_sig_0[np.isfinite(flat_sig_0)])[:-7]), NUM_LEVELS), cmap='Purples_r')
    axc = plt.colorbar(c)
    plt.xlim(-0.11, 0.11)
    plot_pt()
    plt.ylim(-0.24, -0.03)
    axc.set_label("$\\sigma(\\gamma_0)$ ($\\times 10^{-6}$)")
    plt.xlabel("$K_{22}$")
    plt.ylabel("$K_{20}$")
    plt.tight_layout()
    plt.savefig("compile-figs/theta-1-sigma.pdf")
    plt.savefig("compile-figs/theta-1-sigma.png")

    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, np.array(sigma1_data) * 1e6, levels=np.linspace(0, np.nanmax(sigma1_data) * 1e6, NUM_LEVELS), cmap='Purples_r')
    axc = plt.colorbar(c)
    plt.xlim(-0.11, 0.11)
    plot_pt()
    plt.ylim(-0.24, -0.03)
    axc.set_label("$\\sigma(K_{22})$ ($\\times 10^{-6}$)")
    plt.xlabel("$K_{22}$")
    plt.ylabel("$K_{20}$")
    plt.tight_layout()
    plt.savefig("compile-figs/theta-2-sigma.pdf")
    plt.savefig("compile-figs/theta-2-sigma.png")

    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, np.array(sigma2_data) * 1e6, levels=np.linspace(0, np.nanmax(sigma2_data) * 1e6, NUM_LEVELS), cmap='Purples_r')
    axc = plt.colorbar(c)
    plt.xlim(-0.11, 0.11)
    plot_pt()
    plt.ylim(-0.24, -0.03)
    axc.set_label("$\\sigma(K_{20})$ ($\\times 10^{-6}$)")
    plt.xlabel("$K_{22}$")
    plt.ylabel("$K_{20}$")
    plt.tight_layout()
    plt.savefig("compile-figs/theta-3-sigma.pdf")
    plt.savefig("compile-figs/theta-3-sigma.png")


if __name__ == "__main__":
    covariance()
    plt.show()
    #plot_pdf(1)
    #plt.show()
    #plot_pdf(2)
    #plt.show()
