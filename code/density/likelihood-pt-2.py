from re import S
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

def get_r_mat(n_elements):
    r_mat = np.zeros(((MAX_L + 1)**2 + 1, n_elements), dtype=np.complex)
    row = 0
    for l in range(MAX_L + 1):
        for m in range(-l, l+1):
            for col in range(n_elements):
                r_mat[row][col] = rlm(l, m, pos_cache[col]) / A_M**l
            row += 1
    for col in range(n_elements):
        r_mat[row][col] = np.sum(pos_cache[col]**2) / A_M**2
    return np.array(r_mat)

def get_a_mat(n_elements):
    a_mat = np.zeros(((MAX_L + 1)**2 + 1, n_elements), dtype=np.complex)
    row = 0
    for l in range(MAX_L + 1):
        for m in range(-l, l+1):
            for col in range(n_elements):
                a_mat[row][col] = complex_klms[row]
            row += 1
    for col in range(n_elements):
        a_mat[row][col] = 1
    return np.array(a_mat)

def get_element_densities(n_elements):
    key_mat = get_r_mat(n_elements) - get_a_mat(n_elements)
    ones = np.ones((n_elements, 1))
    densities = np.matmul(pinv(key_mat), np.matmul(key_mat, -ones)) + ones
    return densities

def write_densities(element_densities):
    densities = np.full((len(pos_array), len(pos_array), len(pos_array)), np.nan)
    num = 0
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if indicator([x, y, z]):
                    densities[nx, ny, nz] = element_densities[num].real
                    num += 1
    return densities

if __name__ == "__main__":
    n_elements = get_number_elements()
    print("Generated cells")


    element_densities = get_element_densities(n_elements)
    print("Got densities")

    densities = write_densities(element_densities)
    print(f"Wrote densities")

    with open("data/"+TAG+"-likelihood.dat", 'wb') as f:
        np.save(f, densities)
    
    mass = get_mass(densities)
    radius = get_radius(densities) / np.sqrt(mass)
    print(f"Radius: {radius / A_M}\tTrue radius: {1}")

    for l in range(MAX_L+1):
        print()
        for m in range(-l, l+1):
            print(f"Klm: {get_hlms(l, m, densities) / mass / A_M**l}\tTrue Klm: {complex_klms[get_index(l, m)]}")
    
    for percentile in range(0, 101, 10):
        print(f"Percentile {percentile}", np.nanpercentile(densities, percentile))
    print(np.nanmean(densities))
    