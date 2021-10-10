from matplotlib import pyplot as plt
from scipy import linalg
from scipy import optimize
import numpy as np
from mpmath import mp

n = 2
quantize_num = 100
NUM_POINTS = 100000
LIMIT = 2


SIGMAS = [[10, 8], [8, 30]]

def log_like(dxs, sigmas):
    return -np.sum([np.sum([dxs[i] * dxs[j] * sigmas[j][i] / 2 for i in range(n)]) for j in range(n)])

def minus_log_like(dxs, sigmas):
    return -log_like(dxs, sigmas)

bfgs_min = optimize.minimize(minus_log_like, [0.1, 0.1], args=(SIGMAS), method='L-BFGS-B')
print(bfgs_min)
hess_inv = bfgs_min.hess_inv.todense()
print(linalg.inv(hess_inv))


def populate(sigmas, count):
    sigmas=np.asarray([np.asarray(s) for s in sigmas])
    evals, diagonalizer = mp.eigsy(mp.matrix(sigmas))
    diagonalizer = np.array(diagonalizer.tolist(),dtype=np.float32)
    evals = np.asarray([float(e) for e in evals])
    diagonal_points = 1/np.sqrt(np.abs(evals)) * (np.random.randn(count * n).reshape(count, n))
    global_points = np.asarray([np.matmul(diagonalizer, d) for d in diagonal_points])
    return global_points.transpose()

def populate_inv(hess_inv, count):
    evals, diagonalizer = mp.eigsy(mp.matrix(hess_inv))
    diagonalizer = np.array(diagonalizer.tolist(),dtype=np.float32)
    evals = np.asarray([float(e) for e in evals])
    if np.any(evals < 0):
        raise Exception("A hessian's evals were negative")
    diagonal_points = np.sqrt(evals) * (np.random.randn(count * n).reshape(count, n))
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

def quantize(points):
    grid = np.zeros((quantize_num, quantize_num))
    for p in points.transpose():
        i = (p + LIMIT) / (2 * LIMIT) * quantize_num
        i = np.minimum(quantize_num-1, np.maximum(0, i))
        grid[int(i[0]), int(i[1])] += 1
    return grid


side_x = np.linspace(-LIMIT, LIMIT, 100)

for i in range(n):
    for j in range(i+1, n):
        plt.figure()
        data = [[log_like(make_x(xi, xj, i, j), SIGMAS) for xi in side_x] for xj in side_x]
        c = plt.pcolor(side_x, side_x, data)
        plt.xlabel(str(i))
        plt.ylabel(str(j))
        plt.colorbar(c)

        pts = populate_inv(hess_inv, NUM_POINTS)
        plt.scatter(pts[i], pts[j], marker='.', alpha=0.5, s=10)
        #pts = populate_ball(sigmas, 1000)
        #plt.scatter(pts[i], pts[j], marker='.', alpha=0.5, s=10)
        plt.xlim(-LIMIT, LIMIT)
        plt.ylim(-LIMIT, LIMIT)

        plt.savefig("out.png")

        plt.figure()
        likes = []
        for x in np.linspace(LIMIT / quantize_num-LIMIT, LIMIT-LIMIT/ quantize_num, quantize_num):
            line = []
            for y in np.linspace(LIMIT / quantize_num-LIMIT, LIMIT-LIMIT/ quantize_num, quantize_num):
                line.append(np.exp(log_like([x, y], SIGMAS)))
            likes.append(np.array(line))
        likes = np.asarray(likes)
        grid = quantize(pts)

        print(np.sum(likes))
        rescale = NUM_POINTS / np.sum(likes)
        counts = np.rot90((grid - likes * rescale) / np.sqrt(grid), k=1)

        c = plt.imshow(counts)
        plt.colorbar(c)

        plt.figure()
        stripped_counts = counts.reshape(quantize_num * quantize_num)
        stripped_counts = stripped_counts[grid.reshape(quantize_num * quantize_num) > 15]
        stripped_counts = stripped_counts[np.isfinite(stripped_counts)]
        x = np.linspace(-5, 5, 100)
        bins = np.linspace(-5, 5, 50)
        plt.hist(stripped_counts, bins=bins, density=True)
        plt.plot(x, 1 / np.sqrt(2 * np.pi) * np.exp(-x**2/2))
        plt.xlim(-5, 5)

plt.show()
