from stringprep import in_table_a1
import numpy as np
from scipy.special import ellipk
from scipy.optimize import root

DT = 0.1 # hours
MAX_THETA_SEED = 0.8

pPsiObs = 264.178
pPhiObs = 27.38547
i1 = 0.64
i2 = 0.97
i3 = 1

i_plus = 0.5 * (1/i1 + 1/i2)
i_minus = 0.5 * (1/i1 - 1/i2)
delta_phi = 2 * np.pi * pPsiObs / pPhiObs

def get_delta_phi(theta0, psi0, pPsi, l):
    integral = 0
    theta = theta0
    psi = psi0
    for t in np.arange(0, pPsi/2, DT):
        phidot = l * (i_plus - i_minus * np.cos(2 * psi))
        thetadot = l * i_minus * np.sin(theta) * np.sin(2 * psi)
        psidot = np.cos(theta) * (l / i3 - phidot)
        integral += phidot * DT
        theta += thetadot * DT
        psi += psidot * DT
    return 2 * integral

def get_energy(l, theta, psi):
    return 0.5 * l * l * (np.sin(theta)**2 * (np.sin(psi)**2 / i1 + np.cos(psi)**2 / i2) + np.cos(theta)**2 / i3)

def get_ksq(l, theta, psi):
    energy = get_energy(l, theta, psi)
    return (i2-i1) * (l * l / (2 * energy) - i3) / ((i1 - i3) * (i2 - l* l / (2 * energy)))

def get_pPsi(l, theta, psi):
    energy = get_energy(l, theta, psi)
    ksq = get_ksq(l, theta, psi)
    if ksq > 1:
        raise Exception("Incorrect axis choice")
    return 4 * np.sqrt(i1 * i2 * i3 / (2 * energy * (i1-i3) * (i2- l*l/(2 * energy)))) * ellipk(ksq)

def deviation(l, theta0, psi0):
    return (get_pPsi(l, theta0, psi0) - pPsiObs, get_delta_phi(theta0, psi0, pPsiObs, l) - delta_phi)

def get_with_seed(args):
    theta_guess, psi0 = args
    l_guess = 0.204#(i1 + i2 + i3) / 3 * 2 * np.pi / pPsiObs
    def obj_fn(x):
        try:
            return deviation(x[0], x[1], psi0)
        except Exception:
            return (np.inf, np.inf)
    result = root(obj_fn, (l_guess, theta_guess))
    if not result.success:
        return np.nan
    if np.any(result.fun > 1e-3):
        return np.nan
    l, theta0 = result.x
    theta0 = theta0 % np.pi
    if theta0 > np.pi / 2:
        theta0 = np.pi - theta0
    return theta0


def get_theta(psi0, pool):
    print(psi0)
    thetas = [np.nan]
    while np.isnan(np.nanmin(thetas)):
        seeds = np.random.random(8) * MAX_THETA_SEED
        thetas = pool.map(get_with_seed, zip(seeds, np.ones_like(seeds) * psi0))
    return np.nanmin(thetas)# Just find an entry without nans

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from multiprocessing import Pool
    psis = np.linspace(0, 2 * np.pi, 4000)
    with Pool() as pool:
        thetas = np.array([get_theta(psi0, pool) for psi0 in psis])
    with open("psis-thetas.txt", 'w') as f:
        np.savetxt(f, np.array([psis, thetas]))
    plt.plot(psis, thetas)
    plt.xlabel("Initial roll")
    plt.ylabel("Theta")
    plt.savefig("thetas.png")
    plt.show()