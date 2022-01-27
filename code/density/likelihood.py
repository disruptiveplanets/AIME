import numpy as np
from scipy.linalg import pinv, norm
from setup import *

max_index = None
pos_cache = []

def get_number_elements():
    global pos_cache
    num = 0
    for x in pos_array:
        for y in pos_array:
            for z in pos_array:
                if indicator([x, y, z]):
                    pos_cache.append(np.array([x, y, z]))
                    num += 1
    return num

def get_a(row, col, n_elements):
    l, m = lm_from_index(row)
    if l == MAX_L + 1:
        return norm(pos_cache[col])**2
    return rlm(l, m, pos_cache[col])

def get_c(n_elements):
    a = []
    for row in range(len(complex_hlms) + 1):
        line = []
        for col in range(n_elements):
            line.append(get_a(row, col, n_elements))
        a.append(line)
    a = np.array(a)
    print(f"Generated a (shape {a.shape})")

    p = pinv(a)
    print(f"Generated pseudoinverse")

    k = np.append(np.array(complex_hlms), A_M**2) - np.matmul(a, np.ones(n_elements) / n_elements)
    c = np.matmul(p, k)

    return c + 1 / n_elements

def write_densities(c):
    densities = np.full((len(pos_array), len(pos_array), len(pos_array)), np.nan)
    num = 0
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if indicator([x, y, z]):
                    densities[nx, ny, nz] = c[num].real / DIVISION**3
                    num += 1
    return densities

if __name__ == "__main__":
    n_elements = get_number_elements()
    print("Generated cells")

    c = get_c(n_elements)
    print("Got densities")

    densities = write_densities(c)
    print(f"Wrote densities")

    with open("data/"+TAG+"-likelihood.dat", 'wb') as f:
        np.save(f, densities)
    
    radius = get_radius(densities)
    print(f"Radius: {radius}\tTrue radius: {A_M}")
    for l in range(MAX_L+1):
        for m in range(-l, l+1):
            print(f"Hlm: {get_hlms(l, m, densities)}\tTrue Hlm: {complex_hlms[get_index(l, m)]}")
    