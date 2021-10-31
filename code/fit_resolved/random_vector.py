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

# Make an un_uniform randomizer
def randomize_rotate_uniform_err(theta, phi, sigma):
    return np.sqrt(np.array([
            1/2 * np.sinh(sigma**2) * (np.cos(theta)**2 * np.cos(phi)**2 + np.sin(phi)**2) + (np.cosh(sigma**2) - 1) * np.cos(phi)**2 * np.sin(theta)**2,
            1/2 * np.sinh(sigma**2) * (np.cos(theta)**2 * np.sin(phi)**2 + np.cos(phi)**2) + (np.cosh(sigma**2) - 1) * np.sin(phi)**2 * np.sin(theta)**2,
            1/2 * np.sinh(sigma**2) * np.sin(theta)**2 + (np.cosh(sigma**2) - 1) * np.cos(theta)**2]))

# Do the sophisticated, turning error
def randomize_rotate_uniform(y, sigma):
    newy = []
    yerr = []
    assert(len(y) % 3 == 0)
    for i in range(0, len(y), 3):
        norm = np.sqrt(y[i]**2 + y[i+1]**2 + y[i+2]**2)
        theta = np.arccos(y[i+2] / norm)
        phi = np.arctan2(y[i+1], y[i])
        rot_mat = np.matmul(
            np.array([[np.cos(-theta), 0, np.sin(-theta)], [0, 1, 0], [-np.sin(-theta), 0, np.cos(-theta)]]),
            np.array([[np.cos(-phi), -np.sin(-phi), 0], [np.sin(-phi), np.cos(-phi), 0], [0, 0, 1]])
        )
        tilt_phi = np.random.uniform() * 2 * np.pi
        tilt_theta = np.random.randn() * sigma
        untilt_vec = [np.sin(tilt_theta) * np.cos(tilt_phi), np.sin(tilt_theta) * np.sin(tilt_phi), np.cos(tilt_theta)]
        newvec = np.matmul(rot_mat.transpose(), untilt_vec) * norm
        errs = randomize_rotate_uniform_err(theta, phi, sigma)
        yerr.append(errs[0] * abs(y[i]))
        yerr.append(errs[1] * abs(y[i+1]))
        yerr.append(errs[2] * abs(y[i+2]))
        newy.append(newvec[0])
        newy.append(newvec[1])
        newy.append(newvec[2])

    return np.asarray(newy), np.asarray(yerr)

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


def randomize_flat(y, sigma):
    yerr = np.ones_like(y) * sigma
    return y + np.random.randn(len(y)) * sigma, yerr
