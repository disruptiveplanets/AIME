import sys
import matplotlib.pyplot as plt
from shutil import copyfile
import numpy as np
from numdifftools import Hessian
sys.path.append("../../fit_resolved")
import asteroids_0_2 as asteroids
import random_vector


plt.style.use("jcap")

SYMMETRIC = True
GAMMA_0 = 0.39269908169

if SYMMETRIC:
    TRUE_PARAMS = np.array([0, -0.09766608])
    SIGMA = 0.1

else:
    TRUE_PARAMS = np.array([0.05200629, -0.2021978])
    SIGMA = 0.1


NUM_TRIALS = 1000
SPIN_RESOLUTION = 20

SPIN = [0.00013712601, 0, -0.00013712601]
EARTH_RADIUS = 6_378_000
GM = 3.986004418e14
IMPACT_PARAMETER = 5 * EARTH_RADIUS
SPEED = 6000
RADIUS = 1000
JLMS = [1]
CADENCE = 120

data = np.array(asteroids.simulate(CADENCE, JLMS, TRUE_PARAMS, RADIUS, SPIN[0], SPIN[1], SPIN[2],
    GAMMA_0, IMPACT_PARAMETER, SPEED, GM, EARTH_RADIUS, -1, False)).reshape(-1, 3)
err_y, y_inv_covs = random_vector.randomize_rotate_uniform(data, SIGMA)

def hess_func(theta):
    try:
        model = np.array(asteroids.simulate(CADENCE, JLMS, theta, RADIUS, SPIN[0], SPIN[1], SPIN[2],
            GAMMA_0, IMPACT_PARAMETER, SPEED, GM, EARTH_RADIUS, -1, False)).reshape(-1, 3)
    except Exception:
        return np.inf
    chisq = 0
    for i in range(len(data)):
        chisq += np.matmul(data[i] - model[i], np.matmul(y_inv_covs[i], data[i] - model[i]))
    return -chisq / 2.0
    

def populate(evals, diagonalizer, count):
    spacing = 20 / np.sqrt(evals)
    diagonal_points = spacing * (np.random.randn(count * len(TRUE_PARAMS)).reshape(count, len(TRUE_PARAMS)))
    global_points = np.asarray([np.matmul(diagonalizer.transpose(), d) for d in diagonal_points])
    return global_points

def gen_data():
    HESSIAN = -Hessian(hess_func)(TRUE_PARAMS)
    print("Hessian", HESSIAN)

    evals, diagonalizer = np.linalg.eigh(HESSIAN)
    print(diagonalizer)
    print("Evals", evals)
    if np.any(evals < 0):
        raise Exception(f"Some eigenvalues were negative ({evals})")
    results = []
    trial = 0
    while trial < NUM_TRIALS:
        print(f"{trial}/{NUM_TRIALS}")
        if trial == 0:
            theta = TRUE_PARAMS
        else:
            theta = TRUE_PARAMS + populate(evals, diagonalizer, 1)[0]
        try:
            resolved_data = asteroids.simulate(CADENCE, JLMS, theta, RADIUS, SPIN[0], SPIN[1],
                SPIN[2], GAMMA_0, IMPACT_PARAMETER, SPEED, GM, EARTH_RADIUS, -1, False)
        except Exception:
            continue
        results.append(resolved_data) # x1, y1, z1, x2, y2, z2, ...
        trial += 1

    if SYMMETRIC:
        with open("data_runs-sym.npy", 'wb') as f:
            np.save(f, results)
    else:
        with open("data_runs-asym.npy", 'wb') as f:
            np.save(f, results)

def plot_data():
    with open("data_runs.npy", 'rb') as f:
        runs = np.load(f)
    time = np.arange(0, len(runs[0])//3) *  CADENCE / 3600.0
    time -= time[-1] / 2
    spins = np.linspace(np.min(runs), np.max(runs), SPIN_RESOLUTION + 1)
    spins_middle = np.array([(spins[i] + spins[i+1]) / 2 for i in range(SPIN_RESOLUTION)])
    x_hist = []
    y_hist = []
    z_hist = []
    for i in range(len(time)):
        counts_x, _ = np.histogram(runs[:, 3 * i], bins=spins)
        counts_y, _ = np.histogram(runs[:, 3 * i + 1], bins=spins)
        counts_z, _ = np.histogram(runs[:, 3 * i + 2], bins=spins)
        x_hist.append(counts_x)
        y_hist.append(counts_y)
        z_hist.append(counts_z)

    plt.figure()
    #plt.plot(time, runs[0, ::3], label="$\omega_x$", color="black")
    #plt.plot(time, runs[0, 1::3], label="$\omega_y$")
    #plt.plot(time, runs[0, 2::3], label="$\omega_z$")

    reds = np.array([[1, 0, 0, a] for a in np.linspace(0, 1, 4)])
    greens = np.array([[0, 1, 0, a] for a in np.linspace(0, 1, 4)])
    blues = np.array([[0, 0, 1, a] for a in np.linspace(0, 1, 4)])

    plt.contourf(time, spins_middle, np.array(x_hist).transpose(), colors=reds)
    plt.contourf(time, spins_middle, np.array(y_hist).transpose(), colors=greens)
    plt.contourf(time, spins_middle, np.array(z_hist).transpose(), colors=blues)

    plt.legend()
    plt.xlabel("Time to periapsis (h)")
    plt.ylabel("Angular velocity (rad / s)")
    plt.tight_layout()
    if SYMMETRIC:
        plt.savefig("spin-chaos-sym.pdf")
    else:
        plt.savefig("spin-chaos-asym.pdf")

def plot_sigma():
    if SYMMETRIC:
        with open("data_runs-sym.npy", 'rb') as f:
            runs = np.load(f)
    else:
        with open("data_runs-asym.npy", 'rb') as f:
            runs = np.load(f)
    time = np.arange(0, len(runs[0])//3) *  CADENCE / 3600.0
    time -= time[-1] / 2
    spins = np.linspace(np.min(runs), np.max(runs), SPIN_RESOLUTION + 1)
    spins_middle = np.array([(spins[i] + spins[i+1]) / 2 for i in range(SPIN_RESOLUTION)])
    lows_1 = []
    lows_2 = []
    lows_3 = []
    highs_1 = []
    highs_2 = []
    highs_3 = []

    for i in range(len(time)):
        low_x_1 = np.percentile(runs[:, 3 * i], 100-68.27)
        high_x_1 = np.percentile(runs[:, 3 * i], 68.27)
        low_y_1 = np.percentile(runs[:, 3 * i + 1], 100-68.27)
        high_y_1 = np.percentile(runs[:, 3 * i + 1], 68.27)
        low_z_1 = np.percentile(runs[:, 3 * i + 2], 100-68.27)
        high_z_1 = np.percentile(runs[:, 3 * i + 2], 68.27)

        low_x_2 = np.percentile(runs[:, 3 * i], 100-95.45)
        high_x_2 = np.percentile(runs[:, 3 * i], 95.45)
        low_y_2 = np.percentile(runs[:, 3 * i + 1], 100-95.45)
        high_y_2 = np.percentile(runs[:, 3 * i + 1], 95.45)
        low_z_2 = np.percentile(runs[:, 3 * i + 2], 100-95.45)
        high_z_2 = np.percentile(runs[:, 3 * i + 2], 95.45)

        low_x_3 = np.percentile(runs[:, 3 * i], 100-99.73)
        high_x_3 = np.percentile(runs[:, 3 * i], 99.73)
        low_y_3 = np.percentile(runs[:, 3 * i + 1], 100-99.73)
        high_y_3 = np.percentile(runs[:, 3 * i + 1], 99.73)
        low_z_3 = np.percentile(runs[:, 3 * i + 2], 100-99.73)
        high_z_3 = np.percentile(runs[:, 3 * i + 2], 99.73)

        lows_1.append((low_x_1, low_y_1, low_z_1))
        highs_1.append((high_x_1, high_y_1, high_z_1))
        lows_2.append((low_x_2, low_y_2, low_z_2))
        highs_2.append((high_x_2, high_y_2, high_z_2))
        lows_3.append((low_x_3, low_y_3, low_z_3))
        highs_3.append((high_x_3, high_y_3, high_z_3))

    lows_1 = np.transpose(lows_1) * 3600
    highs_1 = np.transpose(highs_1) * 3600
    lows_2 = np.transpose(lows_2) * 3600
    highs_2 = np.transpose(highs_2) * 3600
    lows_3 = np.transpose(lows_3) * 3600
    highs_3 = np.transpose(highs_3) * 3600

    plt.figure(figsize=(10, 4))
    plt.plot(time, runs[0, ::3] * 3600, color="black")
    plt.plot(time, runs[0, 1::3] * 3600, color="black")
    plt.plot(time, runs[0, 2::3] * 3600, color="black")
    #plt.scatter(time, err_y[:,0] * 3600, color="black")
    #plt.scatter(time, err_y[:,1] * 3600, color="black")
    #plt.scatter(time, err_y[:,2] * 3600, color="black")
    plt.plot([], [], label="$\omega_x$", color="C0")
    plt.plot([], [], label="$\omega_y$", color="C1")
    plt.plot([], [], label="$\omega_z$", color="C2")

    #reds = np.array([[1, 0, 0, a] for a in np.linspace(0, 1, 4)])
    #greens = np.array([[0, 1, 0, a] for a in np.linspace(0, 1, 4)])
    #blues = np.array([[0, 0, 1, a] for a in np.linspace(0, 1, 4)])
    plt.fill_between(time, lows_1[0], highs_1[0], color="C0", alpha=1)
    plt.fill_between(time, lows_2[0], highs_2[0], color="C0", alpha=0.5)
    plt.fill_between(time, lows_3[0], highs_3[0], color="C0", alpha=0.33)

    plt.fill_between(time, lows_1[1], highs_1[1], color="C1", alpha=1)
    plt.fill_between(time, lows_2[1], highs_2[1], color="C1", alpha=0.5)
    plt.fill_between(time, lows_3[1], highs_3[1], color="C1", alpha=0.33)

    plt.fill_between(time, lows_1[2], highs_1[2], color="C2", alpha=1)
    plt.fill_between(time, lows_2[2], highs_2[2], color="C2", alpha=0.5)
    plt.fill_between(time, lows_3[2], highs_3[2], color="C2", alpha=0.33)

    #plt.contourf(time, spins_middle, np.array(x_hist).transpose(), colors=reds)
    #plt.contourf(time, spins_middle, np.array(y_hist).transpose(), colors=greens)
    #plt.contourf(time, spins_middle, np.array(z_hist).transpose(), colors=blues)

    plt.legend(loc=(0.02, 0.6), fontsize=12)
    plt.xlabel("Time to periapsis (hr)")
    plt.ylabel("Angular velocity (rad / hr)")
    plt.xlim(-5, np.max(time))
    plt.tight_layout()
    if SYMMETRIC:
        plt.savefig("spin-chaos-sym.pdf")
        plt.savefig("spin-chaos-sym.png")
    else:
        plt.savefig("spin-chaos-asym.pdf")
        plt.savefig("spin-chaos-asym.png")

gen_data()
#plot_data()
plot_sigma()
plt.show()
