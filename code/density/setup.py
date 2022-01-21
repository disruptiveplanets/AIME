import numpy as np
from scipy.special import lpmv, factorial

DIVISION = 100
MAX_RADIUS = 2001
pos_array = np.arange(-MAX_RADIUS, MAX_RADIUS, DIVISION)

KLMS = [
    1.0,
    0, 0, 0,
    0.03, 0, 0, 0, -0.1,
]

A_M = 1000

def indicator(pos):
    return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] < 1000**2

######################################################

def get_index(lp, mp):
    return lp **2 + mp + lp


def lm_from_index(index):
    l = int(np.sqrt(index))
    m = index - l**2 - l
    return (l, m)

MAX_L = int(np.sqrt(len(KLMS))) - 1

real_hlms = [k for k in KLMS]
complex_hlms = np.zeros_like(KLMS, dtype=np.cdouble)
for l in range(MAX_L+1):
    for m in range(-l, l+1):
        if m == 0:
            complex_hlms[get_index(l, m)] = KLMS[(l+1)**2 - 1]
        else:
            complex_hlms[get_index(l, m)] = KLMS[(l)**2 + 2 * l - 2 * abs(m)] + 1.0j * KLMS[(l)**2 + 2 * l - 2 * abs(m) + 1]
        if m < 0:
            complex_hlms[get_index(l, m)] = (-1)**m * complex_hlms[get_index(l, m)].conj()
        real_hlms[get_index(l, m)] *= A_M**l
        complex_hlms[get_index(l, m)] *= A_M**l

def rlm(l, m, pos):
    r = np.sqrt(max(0, pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]))
    return lpmv(m, l, pos[2] / r) / factorial(l + m) * r**l * np.exp(1j * m * np.arctan2(pos[1], pos[0]))

def get_hlms(l, m, densities, radius):
    hlm = 0
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if not np.isnan(densities[nx,ny,nz]):
                    hlm += densities[nx,ny,nz] * rlm(l, m, [x,y,z])
    return hlm * DIVISION**3

def get_radius(densities):
    rad = 0
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if not np.isnan(densities[nx,ny,nz]):
                    rad += densities[nx,ny,nz] * (x**2 + y**2 + z**2)
    return np.sqrt(rad * DIVISION**3)
