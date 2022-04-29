# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

plt.style.use('jcap')

SIDE_LENGTH = 17
FIG_SIZE = (5.333, 4)
NUM_LEVELS = 21

with open("sigmas.npy", 'rb') as f:
    sigmas = np.load(f)

# Plot pdfs
def show_sigmas(sigma_index):
    data = []
    X = []
    Y = []
    for i, y in enumerate(np.linspace(0, -0.25, SIDE_LENGTH)):
        data_line = []
        X_line = []
        Y_line = []
        delta = 0.25 / SIDE_LENGTH / 2 * (SIDE_LENGTH - i - 1)
        for j, x in enumerate(np.linspace(-0.125+delta, 0.125+delta, SIDE_LENGTH)):
            index = i * (i + 1) // 2 + j
            data_line.append(sigmas[index, sigma_index])
            X_line.append(x)
            Y_line.append(y)
        data.append(data_line)
        X.append(X_line)
        Y.append(Y_line)

    data_max = max(np.nanmax(data), -np.nanmin(data))
    levels = np.linspace(-data_max, data_max, NUM_LEVELS)
    print(levels)
    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, data, levels=levels, cmap='Oranges_r')
    axc = plt.colorbar(c)
    axc.set_label("$\\textrm{Corr}(K_{22}, K_{20})$")
    plt.xlim(-0.11, 0.11)
    plt.ylim(-0.24, -0.03)
    plt.xlabel("$K_{22}$")
    plt.ylabel("$K_{20}$")
    plt.tight_layout()
    plt.savefig("compile-figs/corr23.pdf")

if __name__ == "__main__":
    show_sigmas(3)
    plt.show()