import corner
import numpy as np
from scipy.special import lpmv, factorial
from scipy.linalg import norm
from multiprocessing import Pool
from matplotlib import pyplot as plt
from functools import reduce

A_M = 1000 / np.sqrt(5/3) # NOT IMPLEMENTED AS A CONSTRAINT
DIVISION = 50
MAX_RADIUS = 999
pos_array = np.arange(-MAX_RADIUS, MAX_RADIUS, DIVISION)

NUM_JOBS = 10000

def rlm(l, m, pos):
    r = np.sqrt(max(0, pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]))
    return lpmv(m, l, pos[2] / r) / factorial(l + m) * r**l * np.exp(1j * m * np.arctan2(pos[1], pos[0]))

def indicator(pos):
    return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] < A_M * A_M * 5/3

rlms = np.zeros((7, len(pos_array), len(pos_array), len(pos_array)), dtype=np.complex)
for m in range(-3, 4):
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if indicator([x, y, z]):
                    rlms[m+3,nx,ny,nz] = rlm(3, m, [x, y, z])
total_size = len(rlms[0].reshape(-1))

square_mask = np.zeros((len(pos_array), len(pos_array), len(pos_array)))
for nx, x in enumerate(pos_array):
    for ny, y in enumerate(pos_array):
        for nz, z in enumerate(pos_array):
            if indicator([x, y, z]):
                square_mask[nx,ny,nz] = norm([x, y, z])**2

def gen_density(seed):
    np.random.seed(seed)
    freqs = np.arange(1, 20)
    
    x = np.sum(np.random.random(len(freqs)) * np.cos(np.pi / MAX_RADIUS * np.outer(pos_array, freqs) + np.random.random(len(freqs)) * 2 * np.pi), axis=1)
    y = np.sum(np.random.random(len(freqs)) * np.cos(np.pi / MAX_RADIUS * np.outer(pos_array, freqs) + np.random.random(len(freqs)) * 2 * np.pi), axis=1)
    z = np.sum(np.random.random(len(freqs)) * np.cos(np.pi / MAX_RADIUS * np.outer(pos_array, freqs) + np.random.random(len(freqs)) * 2 * np.pi), axis=1)

    densities = np.abs(reduce(np.multiply.outer, [x, y, z]))
    
    #densities = np.abs(perlin3d.generate_perlin_noise_3d(rlms[0].shape, res=np.array(rlms[0].shape) // 8))
    #densities = np.random.random(size=total_size).reshape(rlms[0].shape)
    mass = np.sum(densities)
    radius = np.sqrt(np.sum(densities * square_mask) / mass)
    complex_klms = np.sum(densities * rlms, axis=(1, 2, 3)) / radius**3 / mass
    klms = [complex_klms[i//2].real if i % 2 == 0 else complex_klms[i//2].imag for i in range(7)]
    #if np.all(np.abs(klms) < 1):
    #    print(klms)
    return klms

seeds = [int(i) for i in np.random.random(NUM_JOBS) * 2**32]
with Pool() as pool:
    points = np.array(pool.map(gen_density, seeds))

with open("points.npy", 'wb') as f:
    np.save(f, points)

corner.corner(points, labels="Re(-3) Im(-3) Re(-2) Im(-2) Re(-1) Im(-1) 0".split(' '))
plt.savefig("corner.png")
plt.show()