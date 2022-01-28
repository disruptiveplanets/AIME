from setup import *
from multiprocessing import Pool
import numpy as np
from scipy.linalg import inv, norm

NUM_DISTROS = 48

# Setup provides klm, radius, and indicator

def get_seed_points():
    seeds = []
    while len(seeds) < len(complex_hlms) + 1:
        new_seed = (np.random.random(3) - 0.5) * 2 * MAX_RADIUS
        if indicator(new_seed):
            seeds.append(new_seed)
    return np.array(seeds)

def get_seed_index(point, seeds):
    return np.argmin(norm(point - seeds, axis=1))

def get_kilms(l, m, seeds):
    # return an array with same length as seeds containing the klm contributions
    integrals = np.zeros(len(seeds), dtype=np.cdouble)
    for x in pos_array:
        for y in pos_array:
            for z in pos_array:
                if indicator([x, y, z]):
                    if l == MAX_L+1:
                        integrals[get_seed_index([x, y, z], seeds)] += x**2 + y**2 + z**2
                    else:
                        integrals[get_seed_index([x, y, z], seeds)] += rlm(l, m, [x, y, z])
    return integrals * DIVISION**3

def write_densities(rhos, seeds):
    densities = np.full((len(pos_array), len(pos_array), len(pos_array)), np.nan)
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if indicator([x, y, z]):
                    seed_index = get_seed_index([x, y, z], seeds)
                    densities[nx, ny, nz] = rhos[seed_index].real
    return densities

def get_voronoi_density(args):
    seeds, data = args
    kilms = []
    for l in range(MAX_L+1):
        for m in range(-l, l+1):
            kilms.append(get_kilms(l, m, seeds))
    kilms.append(get_kilms(MAX_L+1, 0, seeds))
    kilms = np.array(kilms)
    rhos = np.matmul(inv(kilms), data)
    return write_densities(rhos, seeds)
                
if __name__ == "__main__":
    data = [c for c in complex_hlms]
    data.append(A_M**2)
    data = np.array(data)
    args = []
    for i in range(NUM_DISTROS):
        args.append([get_seed_points(), data])

    with Pool() as pool:
        density_list = pool.map(get_voronoi_density, args)
    
    densities = np.mean(density_list, axis=0)

    with open("data/"+TAG+"-ensemble.dat", 'wb') as f:
        np.save(f, densities)
    
    radius = get_radius(densities)
    print(f"Radius: {radius}\tTrue radius: {A_M}")
    for l in range(MAX_L+1):
        for m in range(-l, l+1):
            print(f"Hlm: {get_hlms(l, m, densities)}\tTrue Hlm: {complex_hlms[get_index(l, m)]}")