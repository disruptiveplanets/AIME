import numpy as np
from scipy.special import lpmv, factorial

DIVISION = 29
LENGTH = 700

RLM_EPSILON = 1e-20
NEGATIVE_VEC = None

MAX_RAD = 1500
GRID_LINE = np.arange(-MAX_RAD, MAX_RAD, DIVISION)

def get_indicator_map(am, k22, k20):
    a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22)
    b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22)
    c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20)
    if a > MAX_RAD or b > MAX_RAD or c > MAX_RAD:
        return None
    x,y,z = np.meshgrid(GRID_LINE, GRID_LINE, GRID_LINE)
    return np.array(x*x/(a*a) + y*y/(b*b) + z*z/(c*c) < 1, dtype=float)

def rlm(l,m,x,y,z):
    r = np.sqrt(np.maximum(RLM_EPSILON, x*x + y*y + z*z))
    return lpmv(m, l, z/r) / factorial(l + m) * r**l * np.exp(1j * m * np.arctan2(y, x))

def rlm_gen(l,m):
    return lambda x,y,z: rlm(l,m,x,y,z)

def map_np(fn):
    x,y,z = np.meshgrid(GRID_LINE, GRID_LINE, GRID_LINE)
    return fn(x,y,z)

def moment_field(max_l):
    n = (max_l+1)**2 + 1
    moments = np.zeros((n, len(GRID_LINE), len(GRID_LINE), len(GRID_LINE)), dtype=np.complex)
    i = 0
    for l in range(0, max_l+1):
        for m in range(-l, l+1):
            moments[i] = map_np(rlm_gen(l,m))
            i += 1
    moments[i] = map_np(lambda x,y,z: x*x + y*y + z*z)
    return moments


rlms = moment_field(3)

def get_max_k3m(klms):
    k22 = klms[0]
    k20 = klms[2]
    indicator_map = get_indicator_map(LENGTH, k22, k20)
    if indicator_map is None:
        return [np.nan] * 7
    real_constants = np.append([0, 0, 0, 0], klms)
    constants = []
    i = 0
    for l in range(0, 4):
        for m in range(-l, l+1):
            if m < 0:
                constants.append(real_constants[i] + 1j * real_constants[i+1])
                i += 2
            elif m == 0:
                constants.append(real_constants[i])
                i += 1
            else:
                other_index = (l)**2 + l - m
                constants.append(np.conj(constants[other_index])*(-1)**m)
    constants.append(0) # No radial component
    constants = np.array(constants)

    densities = np.zeros_like(indicator_map, dtype=np.complex)
    for i, c in enumerate(constants):
        densities += c * np.conj(rlms[i])
    densities *= indicator_map
    constants /= -np.min(densities)
    densities /= -np.min(densities)
    densities += np.ones_like(indicator_map)
    densities *= indicator_map
    densities /= np.sum(densities)

    klms = np.sum(rlms * densities, axis=(1,2,3)) * DIVISION**3
    klms /= klms[0]
    radius = np.sqrt(klms[-1])

    for l in range(0, 4):
        klms[l**2:(l+1)**2] /= radius**l
    
    # i = 0
    # for l in range(0, 4):
    #     for m in range(-l, l+1):
    #         print(f"({l}, {m})\t{klms[i]}")
    #         i += 1
    #     print()
    # print(f"a\t{radius}")

    return [abs(klms[9].real), abs(klms[9].imag), 
        abs(klms[10].real), abs(klms[10].imag), 
        abs(klms[11].real), abs(klms[11].imag), 
        abs(klms[12].real)]

def pick_k2m():
    k20 = 1
    k22 = 2
    while 2 * abs(k22) >= -k20:
        k20 = np.random.random() * -0.25
        k22 = (np.random.random()-0.5) * 0.25
    return k22, k20

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from multiprocessing import Pool

    plt.style.use('jcap')
    N_TRIALS = 500
    LOAD_DATA = True

    if not LOAD_DATA:
        seeds = []
        for _ in range(N_TRIALS):
            k22, k20 = pick_k2m()
            seeds.append(np.append([k22, 0, k20, 0, k22], np.random.randn(7)))

        with Pool() as pool:
           maxes = pool.map(get_max_k3m, seeds)
        with open('maxes.npy', 'wb') as f:
           np.save(f, maxes)
    else:
        with open('maxes.npy', 'rb') as f:
            maxes = np.load(f)

    def display(these_maxes, index):
        plt.figure()
        plt.hist(these_maxes, histtype='step', bins=np.linspace(0, np.nanmax(maxes), 20))
        plt.xlabel(f"Maximum $K_{{3{index}}}$")
        plt.ylabel("Count")
        plt.tight_layout()
        plt.savefig(f"klms-{index}.png")
    display(maxes[:,0:2], 3)
    display(maxes[:,2:4], 2)
    display(maxes[:,4:6], 1)
    display(maxes[:,6], 0)
    plt.show()