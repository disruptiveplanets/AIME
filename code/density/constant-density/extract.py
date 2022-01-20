import numpy as np
from scipy.special import lpmv, factorial
from scipy.optimize import root
from multiprocessing import Pool

INTEGRAL_WIDTH = 200

klms = [
    1.0,
    0, 0, 0,
    0.03, 0, 0, 0, -0.1,
]
radius = 1000

klms_and_radius_squared = np.array(klms)[1:] / klms[0]
klms_and_radius_squared = np.append(klms_and_radius_squared, radius**2)

MAX_L = int(np.sqrt(len(klms))) - 1


def surface_integrate(theta_func, phi_func):
    integral_theta = 0
    integral_phi = 0
    for theta in np.linspace(0, np.pi, INTEGRAL_WIDTH):
        integral_theta += theta_func(theta) * np.sin(theta) * np.pi / INTEGRAL_WIDTH
    for phi in np.linspace(0, 2 * np.pi, 2 * INTEGRAL_WIDTH):
        integral_phi += phi_func(phi) * np.pi / INTEGRAL_WIDTH
    return integral_theta * integral_phi

def rlm_precursor_theta(l, m, theta):
    return 1 / radius ** l / (3 + l) * lpmv(m, l, np.cos(theta)) / factorial(l + m)


def rlm_precursor_phi(l, m, phi):
    return np.exp(1j * m * phi)

def ylm_coeff_theta(lmps, theta):
    # Coefficients
    value = 1
    for lp, mp in lmps:
        value *= (-1)**mp * lpmv(mp, lp, np.cos(theta))
    return value

def ylm_coeff_phi(lmps, phi):
    # Coefficients
    value = 1
    for lp, mp in lmps:
        value *= np.exp(1j * mp * phi)
    return value

def get_clm_coeff(args):
    l, m, lmps = args
    return surface_integrate(lambda theta: rlm_precursor_theta(l, m, theta) * ylm_coeff_theta(lmps, theta),
                            lambda phi: rlm_precursor_phi(l, m, phi) * ylm_coeff_phi(lmps, phi))

def get_volume_coeff(lmps):
    return surface_integrate(lambda theta: 1/3 * ylm_coeff_theta(lmps, theta),
                            lambda phi: ylm_coeff_phi(lmps, phi))

def get_radius_coeff(lmps):
    return surface_integrate(lambda theta: 1/5 * ylm_coeff_theta(lmps, theta),
                            lambda phi: ylm_coeff_phi(lmps, phi))

def make_lmps(l):
    if l == -1:
        return [[]]
    these = []
    old_lmps = make_lmps(l - 1)
    for lmps in old_lmps:
        for lp in range(MAX_L+1):
            for mp in range(-lp, lp+1):
                these.append(lmps + [(lp, mp)])
    return these

def get_index(lp, mp):
    return lp * (lp + 1) // 2 + lp + mp

def generate_integrals():
    integrals = []#(lm, lpmp1, lpmp2, ...)
    num = 0
    for l in range(MAX_L+1):
        for m in range(-l, l+1):
            num += len(make_lmps(l + 2))
    print("Num integrals:", num)

    with Pool() as pool:
        for l in range(MAX_L+1):
            for m in range(-l, l+1):
                args = [(l, m, lmps) for lmps in make_lmps(l+2)]
                line = pool.map(get_clm_coeff, args)
                integrals.append(line)

        #volume_integrals = pool.map(get_volume_coeff, make_lmps(2))
        radius_integrals = pool.map(get_radius_coeff, make_lmps(4))
    return integrals, radius_integrals

def get_volume(clms, volume_integrals):
    term_sum = 0
    for lmps_index, lmps in enumerate(make_lmps(2)):
        prod = volume_integrals[lmps_index]
        for lp, mp in lmps:
            prod *= clms[get_index(lp, mp)]
        term_sum += prod
    return term_sum

def get_radius(clms, radius_integrals):
    term_sum = 0
    for lmps_index, lmps in enumerate(make_lmps(4)):
        prod = radius_integrals[lmps_index]
        for lp, mp in lmps:
            prod *= clms[get_index(lp, mp)]
        term_sum += prod
    return term_sum

def solve_function(clms, integrals, radius_integrals):
    res_klms = []
    for l in range(MAX_L+1):
        for m in range(-l, l+1):
            term_sum = 0
            for lmps_index, lmps in enumerate(make_lmps(l+2)):
                prod = integrals[get_index(l, m)][lmps_index]
                for lp, mp in lmps:
                    prod *= clms[get_index(lp, mp)]
                term_sum += prod
            res_klms.append(term_sum)

    real_klms = []
    for l in range(MAX_L+1):
        for m in range(0, l+1):
            real_klms.append(res_klms[get_index(l, m)].real)
            if m != 0:
                real_klms.append(res_klms[get_index(l, m)].imag)
    real_klms.append(get_radius(clms, radius_integrals).real)

    return np.array(real_klms)[1:] / real_klms[0] - np.array(klms_and_radius_squared)

if __name__ == "__main__":
    integrals, radius_integrals = generate_integrals()

    #print(solve_function([radius * np.sqrt(5/3), 0, 0, 0], integrals, radius_integrals))

    print("Calculated integrals")

    seeds = []
    for i in range(16):
        seeds.append(radius * (np.random.random(len(klms)) - 0.5))

    def multi_solve(seed):
        return root(solve_function, seed, args=(integrals, radius_integrals))

    #res = [multi_solve(s) for s in seeds]
    with Pool() as pool:
        res = pool.map(multi_solve, seeds)
    print("Performed solve")
    for r in res:
        print(r.x, r.fun)