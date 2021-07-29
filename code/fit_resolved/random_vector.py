import numpy as np
import random

def vadd(a, b):
    return [a[i] + b[i] for i in range(len(a))]
def vmul(a, b):
    return [ai * b for ai in a]
def vcross(a, b):
    return [a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]]
def vdot(a, b):
    return np.sum([a[i] * b[i] for i in range(len(a))])
def vnorm(a):
    return np.sqrt(np.sum([ai * ai for ai in a]))


def randomize(y, sigma):
    newy = []
    yerr = []
    assert(len(y) % 3 == 0)
    for i in range(0, len(y), 3):
        norm_squared = y[i]**2 + y[i+1]**2 + y[i+2]**2
        errx = sigma * np.sqrt(1 - y[i]**2 / norm_squared)
        erry = sigma * np.sqrt(1 - y[i+1]**2 / norm_squared)
        errz = sigma * np.sqrt(1 - y[i+2]**2 / norm_squared)
        yerr.append(errx)
        yerr.append(erry)
        yerr.append(errz)
        newy.append(y[i] + errx * (1 - 2 * random.random()))
        newy.append(y[i+1] + erry * (1 - 2 * random.random()))
        newy.append(y[i+2] + errz * (1 - 2 * random.random()))

    return np.asarray(newy), np.asarray(yerr) / np.sqrt(3)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    spin = [1, -2, 3]
    sigma = 1
    xs = []
    ys = []
    zs = []
    err = None
    for i in range(100000):
        newspin, err = randomize(spin, sigma)
        xs.append(newspin[0])
        ys.append(newspin[1])
        zs.append(newspin[2])

    plt.hist(xs, label="x")
    plt.hist(ys, label="y")
    plt.hist(zs, label="z")
    plt.axvline(spin[0], color="k")
    plt.axvline(spin[1], color="k")
    plt.axvline(spin[2], color="k")
    plt.axvline(spin[0] - err[0], color="k")
    plt.axvline(spin[1] - err[1], color="k")
    plt.axvline(spin[2] - err[2], color="k")
    plt.axvline(spin[0] + err[0], color="k")
    plt.axvline(spin[1] + err[1], color="k")
    plt.axvline(spin[2] + err[2], color="k")
    plt.legend()
    plt.show()
