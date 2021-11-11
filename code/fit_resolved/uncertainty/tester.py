import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.transform import Rotation as R

SIGMA = 0.1
SHAPE = 1000000

def randomize_rotate_uniform(vec):
    theta = np.arccos(vec[2] / np.linalg.norm(vec))
    phi = np.arctan2(vec[1], vec[0])
    rot_mat = np.matmul(
        np.array([[np.cos(-theta), 0, np.sin(-theta)], [0, 1, 0], [-np.sin(-theta), 0, np.cos(-theta)]]),
        np.array([[np.cos(-phi), -np.sin(-phi), 0], [np.sin(-phi), np.cos(-phi), 0], [0, 0, 1]])
    )
    vecs = []
    for i in range(SHAPE):
        theta = np.random.uniform() * 2 * np.pi
        phi = np.random.randn() * SIGMA
        displace = rot_mat.transpose()
        displace = np.matmul(displace, np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1]]))
        displace = np.matmul(displace, np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]]))
        displace = np.matmul(displace, rot_mat)
        vecs.append(np.matmul(displace, vec))
    return np.array(vecs)

def randomize_rotate_uniform_err(vec, sigma):
    #Old method, with no covariance
    '''return np.sqrt(np.array([
            1/2 * np.sinh(sigma**2) * (np.cos(theta)**2 * np.cos(phi)**2 + np.sin(phi)**2) + (np.cosh(sigma**2) - 1) * np.cos(phi)**2 * np.sin(theta)**2,
            1/2 * np.sinh(sigma**2) * (np.cos(theta)**2 * np.sin(phi)**2 + np.cos(phi)**2) + (np.cosh(sigma**2) - 1) * np.sin(phi)**2 * np.sin(theta)**2,
            1/2 * np.sinh(sigma**2) * np.sin(theta)**2 + (np.cosh(sigma**2) - 1) * np.cos(theta)**2]))'''
    # New method, with covariance.
    diag = np.linalg.norm(vec)**2/4 * (1 - np.exp(-2 * sigma**2))
    const = 1/4 * (1 - np.exp(-sigma**2)) * (1 - 3 * np.exp(-sigma**2))
    return np.array([
            [diag + const * vec[0]**2, const * vec[0] * vec[1], const * vec[0] * vec[2]],
            [const * vec[0] * vec[1], diag + const * vec[1]**2, const * vec[1] * vec[2]],
            [const * vec[0] * vec[2], const * vec[1] * vec[2], diag + const * vec[2]**2]])

def display_one(vec):
    vecs = randomize_rotate_uniform(vec).transpose()
    bins = np.linspace(-np.linalg.norm(vec), np.linalg.norm(vec), 100)
    plt.hist(vecs[0], label='x', bins=bins)
    plt.hist(vecs[1], label='y', bins=bins)
    plt.hist(vecs[2], label='z', bins=bins)
    errs = randomize_rotate_uniform_err(vec, SIGMA)
    print(np.cov(vecs))
    print(errs)
    plt.legend()
    plt.show()


def display_all():
    phis = np.linspace(0, 2 * np.pi, 20)
    thetas = np.linspace(0.0001, np.pi, 10)
    x_unc = []
    y_unc = []
    z_unc = []
    for t in thetas:
        x_line = []
        y_line = []
        z_line = []
        for p in phis:
            vec = np.array([np.sin(t) * np.cos(p), np.sin(t) * np.sin(p), np.cos(t)])
            vecs = randomize_rotate_uniform(vec).transpose()
            err = randomize_rotate_uniform_err(t, p, SIGMA)
            x_line.append(np.std(vecs[0]) - err[0])
            y_line.append(np.std(vecs[1]) - err[1])
            z_line.append(np.std(vecs[2]) - err[2])
        x_unc.append(x_line)
        y_unc.append(y_line)
        z_unc.append(z_line)
    plt.figure()
    c = plt.pcolor(phis, thetas, x_unc)
    plt.xlabel("phi")
    plt.ylabel("theta")
    plt.title("X")
    plt.colorbar(c)
    plt.figure()
    c = plt.pcolor(phis, thetas, y_unc)
    plt.xlabel("phi")
    plt.ylabel("theta")
    plt.title("Y")
    plt.colorbar(c)
    plt.figure()
    c = plt.pcolor(phis, thetas, z_unc)
    plt.xlabel("phi")
    plt.ylabel("theta")
    plt.title("Z")
    plt.colorbar(c)
    plt.show()


if __name__ == "__main__":
    display_one([0.1, 0.4, -0.532151])
    #display_all()
