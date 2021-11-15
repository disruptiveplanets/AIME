import numpy as np
import random
import scipy.linalg

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
def randomize_rotate_uniform_err(spin, sigma):
    norm2 = spin[0]**2 + spin[1]**2 + spin[2]**2
    return 0.25 * (1 - np.exp(-sigma**2)) * (1 - 3 * np.exp(-sigma**2)) * np.array([
        [spin[0]**2, spin[0] * spin[1], spin[0] * spin[2]],
        [spin[0] * spin[1], spin[1]**2, spin[1] * spin[2]],
        [spin[0] * spin[2], spin[1] * spin[2], spin[2]**2]])\
    + norm2 / 4 * (1 - np.exp(-2 * sigma**2)) * np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])


# Do the sophisticated, turning error
def randomize_rotate_uniform(data, sigma):
    newy = []
    ycovs = []
    for x, y, z in data:
        norm = np.sqrt(x**2 + y**2 + z**2)
        theta = np.arccos(z / norm)
        phi = np.arctan2(y, x)
        rot_mat = np.matmul(
            np.array([[np.cos(-theta), 0, np.sin(-theta)], [0, 1, 0], [-np.sin(-theta), 0, np.cos(-theta)]]),
            np.array([[np.cos(-phi), -np.sin(-phi), 0], [np.sin(-phi), np.cos(-phi), 0], [0, 0, 1]])
        )
        tilt_phi = np.random.uniform() * np.pi
        tilt_theta = np.random.randn() * sigma
        untilt_vec = [np.sin(tilt_theta) * np.cos(tilt_phi), np.sin(tilt_theta) * np.sin(tilt_phi), np.cos(tilt_theta)]
        newvec = np.matmul(rot_mat.transpose(), untilt_vec) * norm
        covs = randomize_rotate_uniform_err(newvec, sigma)
        newy.append(newvec)
        #ycovs.append(covs)
        ycovs.append(scipy.linalg.pinvh(covs))
    return np.array(newy), np.array(ycovs)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    spin = [[1, -2, 3]]
    sigma = 0.1
    xs = []
    ys = []
    zs = []
    err = None
    for i in range(10000):
        newspin, err = randomize_rotate_uniform(spin, sigma)
        xs.append(newspin[0][0])
        ys.append(newspin[0][1])
        zs.append(newspin[0][2])

    err_vec = (np.sum(err[0]**2, axis=0))**(1/4)

    print(err)
    print(np.cov([xs, ys, zs]))
    plt.hist(xs, label="x")
    plt.hist(ys, label="y")
    plt.hist(zs, label="z")
    plt.axvline(spin[0], color="k")
    plt.axvline(spin[1], color="k")
    plt.axvline(spin[2], color="k")
    plt.axvline(spin[0] - err_vec[0], color="k")
    plt.axvline(spin[1] - err_vec[1], color="k")
    plt.axvline(spin[2] - err_vec[2], color="k")
    plt.axvline(spin[0] + err_vec[0], color="k")
    plt.axvline(spin[1] + err_vec[1], color="k")
    plt.axvline(spin[2] + err_vec[2], color="k")
    plt.legend()
    plt.show()
