# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

plt.style.use('jcap')

SIDE_LENGTH = 17
FIG_SIZE = (7, 4)
NUM_LEVELS = 21

param_names = ["\\gamma_0", "K_{22}", "K_{20}", "\Re K_{33}", "\Im K_{33}", "\Re K_{32}", "\Im K_{32}", "\Re K_{31}", "\Im K_{31}", "K_{30}"]

def plot_pt():
    plt.scatter([0], [-0.09766608], color='r', marker='o')
    plt.scatter([0.05200629], [-0.2021978], color='k', marker='^')
    plt.scatter([-0.05200629], [-0.2021978], color='k', marker='^')

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
            if j > i:
                data_line.append(np.nan)
            else:
                data_line.append(sigmas[index, sigma_index])
            X_line.append(x)
            Y_line.append(y)
        data.append(data_line)
        X.append(X_line)
        Y.append(Y_line)

    data_max = np.nanmax(data)
    levels = np.linspace(0, data_max, NUM_LEVELS)
    plt.figure(figsize=FIG_SIZE)
    c = plt.contourf(X, Y, data, levels=levels, cmap='Oranges_r')
    axc = plt.colorbar(c)
    axc.set_label(f"$\sigma({param_names[sigma_index]})$")
    plt.xlim(-0.125, 0.125)
    plt.ylim(-0.25, 0.0)
    plt.xlabel("$K_{22}$")
    plt.ylabel("$K_{20}$")
    plot_pt()
    plt.tight_layout()
    plt.savefig(f"sigma-{sigma_index}.png")

if __name__ == "__main__":
    show_sigmas(3)
    show_sigmas(4)
    show_sigmas(5)
    show_sigmas(6)
    show_sigmas(7)
    show_sigmas(8)
    show_sigmas(9)
    plt.show()