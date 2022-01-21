import numpy as np
from scipy.special import lpmv, factorial
from scipy.optimize import root
from multiprocessing import Pool
from setup import *
#from functools import lru_cache

INTEGRAL_WIDTH = 200
NUM_THREADS = 24

realized_hlms = []
for l in range(MAX_L+1):
    for m in range(0, l+1):
        realized_hlms.append(complex_hlms[get_index(l, m)].real)
        if m != 0:
            realized_hlms.append(complex_hlms[get_index(l, m)].imag)
data = np.append(realized_hlms[1:], A_M**2)

def surface_integrate(theta_func, phi_func):
    integral_theta = 0
    integral_phi = 0
    for theta in np.linspace(0, np.pi, INTEGRAL_WIDTH):
        integral_theta += theta_func(theta) * np.sin(theta) * np.pi / INTEGRAL_WIDTH
    for phi in np.linspace(0, 2 * np.pi, 2 * INTEGRAL_WIDTH):
        integral_phi += phi_func(phi) * np.pi / INTEGRAL_WIDTH
    return integral_theta * integral_phi

def get_hlm_coeff(args):
    l, m, quantity_list = args
    assert(np.sum(quantity_list) == 3 + l)
    theta_func, phi_func = get_r_coeff(quantity_list)
    return 1 / (3 + l) * surface_integrate(
        lambda theta: theta_func(theta) * lpmv(m, l, np.cos(theta)) / factorial(l + m),
        lambda phi: phi_func(phi) * np.exp(1j * m * phi))

def get_radius_coeff(quantity_list):
    assert(np.sum(quantity_list) == 5)
    theta_func, phi_func = get_r_coeff(quantity_list)
    return 1 / 5 * surface_integrate(theta_func, phi_func )

def ylms_theta(index):
    l, m = lm_from_index(index)
    return lambda theta: lpmv(m, l, np.cos(theta))

def ylms_phi(index):
    l, m = lm_from_index(index)
    return lambda phi: np.exp(-1j * m * phi)

def get_pascal_coeff(quantity_list):
    return factorial(np.sum(quantity_list)) / factorial(np.sum(quantity_list) - np.sum(quantity_list!=0) + 1)

#@lru_cache(maxsize=None)
def make_quantity_lists(level):
    lists = [[a] for a in range(level+1)]
    for i in range(len(real_hlms)-1):
        new_lists = []
        for item in lists:
            thresh = np.sum(item)
            increment = 0
            while thresh + increment <= level:
                new_lists.append([l for l in item] + [increment])
                increment += 1
        lists = np.array(new_lists)
    return lists[np.sum(lists, axis=1)==level]

#@lru_cache(maxsize=None)
def get_r_coeff(quantity_list):
    prod_funcs_theta = []
    prod_funcs_phi = []
    for index, count in enumerate(quantity_list):
        for c in range(count):
            prod_funcs_theta.append(ylms_theta(index))
            prod_funcs_phi.append(ylms_phi(index))
    coeff = get_pascal_coeff(quantity_list)
    return (lambda theta: np.prod([f(theta) for f in prod_funcs_theta]) * coeff,
        lambda phi: np.prod([f(phi) for f in prod_funcs_phi]))

def generate_integrals():
    hlm_coeffs = []

    with Pool() as pool:
        for l in range(MAX_L+1):
            for m in range(-l, l+1):
                if m < 0:
                    hlm_coeffs.append([])
                    continue
                quantity_lists = make_quantity_lists(l+3)
                args = [(l, m, q) for q in quantity_lists]
                line = pool.map(get_hlm_coeff, args)
                hlm_coeffs.append(line)

        radius_coeffs = pool.map(get_radius_coeff,  make_quantity_lists(5))

    return hlm_coeffs, radius_coeffs


def solve_function(clms, hlm_coeffs, radius_coeffs):
    hlms = []
    for l in range(MAX_L+1):
        for m in range(0, l+1):
            hlm = 0
            for i, q in enumerate(make_quantity_lists(l+3)):
                this_term = hlm_coeffs[get_index(l, m)][i]
                for clm_index, count in enumerate(q):
                    this_term *= clms[clm_index]**count
                hlm += this_term
            hlms.append(hlm.real)
            if m != 0:
                hlms.append(hlm.imag)

    radius = 0
    for i, q in enumerate(make_quantity_lists(5)):
        this_term = radius_coeffs[i]
        for clm_index, count in enumerate(q):
            this_term *= clms[clm_index]**count
        radius += this_term
    
    return np.append(hlms[1:], radius.real) / hlms[0].real - data

def write_densities(clms):
    densities = np.full((len(pos_array), len(pos_array), len(pos_array)), np.nan)
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                r = np.sqrt(x*x + y*y + z*z)
                theta = np.arccos(z / r)
                phi = np.arctan2(y, x)
                limit_r = 0
                for l in range(MAX_L + 1):
                    for m in range(-l, l+1):
                        index = get_index(l, m)
                        limit_r += clms[index] * ylms_phi(index)(phi) * ylms_theta(index)(theta)
                if limit_r.real < 0:
                    return False, None
                if limit_r.real > r:
                    densities[nx,ny,nz] = 1
    return True, densities / (DIVISION**3 * np.nansum(densities))



if __name__ == "__main__":
    hlm_coeffs, radius_coeffs = generate_integrals()
    print("Calculated integrals")
    #print(solve_function([np.sqrt(5/3) * 1000, 0, 0, 0, 0, 0, 0, 0, 0], hlm_coeffs, radius_coeffs))

    seeds = []
    for i in range(NUM_THREADS):
        seeds.append(A_M * (np.random.random(len(real_hlms)) - 0.5))

    def multi_solve(seed):
        res = root(solve_function, seed, args=(hlm_coeffs, radius_coeffs))
        if res.success:
            success, densities = write_densities(res.x)
            if success:
                print(res.x)
                return densities
        return None

    #res = [multi_solve(s) for s in seeds]
    with Pool() as pool:
        res = pool.map(multi_solve, seeds)
    valids = []
    for d in res:
        if d is not None:
            valids.append(d)
    print(f"There were {len(valids)} valid densities with {NUM_THREADS} threads")

    densities = np.mean(valids, axis=0)
    with open("data/"+TAG+"-surface.dat", 'wb') as f:
        np.save(f, densities)
    
    radius = get_radius(densities)
    print(f"Radius: {radius}\tTrue radius: {A_M}")
    for l in range(MAX_L+1):
        for m in range(-l, l+1):
            print(f"Klm: {get_hlms(l, m, densities, radius)}\tTrue Klm: {complex_hlms[get_index(l, m)]}")
    