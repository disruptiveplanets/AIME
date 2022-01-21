import sys
import numpy as np
from scipy.linalg import pinv, norm
from scipy.sparse.linalg import eigsh
from matplotlib import pyplot as plt
from setup import *
from multiprocessing import Pool


def compute_b(l, m, lp, mp):
    integral = 0
    for x in pos_array:
        for y in pos_array:
            for z in pos_array:
                if indicator([x, y, z]):
                    integral += rlm(l, m, [x, y, z]) * rlm(lp, mp, [x, y, z]).conj() 
    return integral * DIVISION**3

def wrapped_b(args):
    l, m = args
    line = []
    for lp in range(MAX_L+1):
        for mp in range(-lp, lp+1):
            line.append(compute_b(l, m, lp, mp))
    return line

def compute_v(l, m):
    integral = 0
    for x in pos_array:
        for y in pos_array:
            for z in pos_array:
                if indicator([x, y, z]):
                    integral += (x*x + y*y + z*z) * rlm(l, m, [x, y, z]).conj() 
    return integral * DIVISION**3

def get_clms():
    v = []
    
    args = []
    for l in range(MAX_L+1):
        for m in range(-l, l+1):
            v.append(compute_v(l, m))
            args.append((l, m))

    with Pool() as pool:
        B = pool.map(wrapped_b, args)
    B.append(v)
    A = np.array(B)

    data = np.append(complex_hlms, A_M**2)
    clms = np.matmul(pinv(A), data)
    return clms

def write_density(clms):
    densities = np.full((len(pos_array), len(pos_array), len(pos_array)), np.nan)
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if indicator([x, y, z]):
                    densities[nx,ny,nz] = 0
                    for l in range(MAX_L+1):
                        for m in range(-l, l+1):
                            densities[nx,ny,nz] += (rlm(l, m, [x, y, z]).conj() * clms[get_index(l, m)]).real
    return densities


if __name__ == "__main__":
    clms = get_clms()
    densities = write_density(clms)

    with open("data/harmonic.dat", 'wb') as f:
        np.save(f, densities)
    
    radius = get_radius(densities)
    print(f"Radius: {radius}\tTrue radius: {A_M}")
    for l in range(MAX_L+1):
        for m in range(-l, l+1):
            print(f"Klm: {get_hlms(l, m, densities, radius)}\tTrue Klm: {complex_hlms[get_index(l, m)]}")
    