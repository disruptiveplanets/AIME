from matplotlib import pyplot as plt
from scipy import linalg
import numpy as np

n = 3


sigmas = np.asarray([np.asarray([1.16596921e+08, 5.71074134e+10, 1.08431228e+10]),
       np.asarray([5.71078263e+10, 2.79705506e+13, 5.31083605e+12]),
       np.asarray([1.08432269e+10, 5.31084862e+12, 1.00838366e+12])])

def log_like(dxs, sigmas):
    return -np.sum([np.sum([dxs[i] * dxs[j] * sigmas[j][i] / 2 for i in range(n)]) for j in range(n)])

def populate(sigmas, count):
    evals, diagonalizer = linalg.eigh(sigmas)
    print(1 / evals)
    diagonal_points = 1/np.sqrt(np.abs(evals)) * (np.random.randn(count * n).reshape(count, n))
    global_points = np.asarray([np.matmul(diagonalizer, d) for d in diagonal_points])
    return global_points.transpose()

def populate_ball(sigmas, count):
    evals, diagonalizer = linalg.eigh(sigmas)
    weights = [np.sum([diagonalizer[i][j] / evals[j] for j in range(n)]) for i in range(n)]
    print(weights)
    global_points = np.sqrt(np.abs(weights)) * (np.random.randn(count * n).reshape(count, n))
    return global_points.transpose()

def make_x(xi, xj, ii, ij):
    out = [0]*n
    out[ii] = xi
    out[ij] = xj
    return out

side_x = np.linspace(-1, 1, 100)

for i in range(n):
    for j in range(i+1, n):
        plt.figure()
        data = [[log_like(make_x(xi, xj, i, j), sigmas) for xi in side_x] for xj in side_x]
        c = plt.pcolor(side_x, side_x, data)
        plt.xlabel(str(i))
        plt.ylabel(str(j))
        plt.colorbar(c)

        pts = populate(sigmas, 1000)
        plt.scatter(pts[i], pts[j], marker='.', alpha=0.5, s=10)
        pts = populate_ball(sigmas, 1000)
        plt.scatter(pts[i], pts[j], marker='.', alpha=0.5, s=10)
        plt.xlim(-1, 1)
        plt.ylim(-1, 1)

plt.show()
