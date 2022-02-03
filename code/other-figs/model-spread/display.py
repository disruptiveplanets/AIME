import sys
import matplotlib.pyplot as plt
from shutil import copyfile
import numpy as np
sys.path.append("../../fit_resolved")
import asteroids_0_2
import random
from multiprocessing import Pool


plt.style.use("jcap")

SYMMETRIC = False
NUM_TRIALS = 1000
EARTH_RADIUS = 6_370_000
AMPLIFY = 700 * 3
FILE_NAME = "../../../data/probe-sigma/sigma-1-47/sigma-1-47"
GM = 3.986004418e14

with open(f"{FILE_NAME}.txt", 'r') as f:
    ASTEROIDS_MAX_J, ASTEROIDS_MAX_K = f.readline().split(', ')
    ASTEROIDS_MAX_J = int(ASTEROIDS_MAX_J)
    ASTEROIDS_MAX_K = int(ASTEROIDS_MAX_K)
    cadence = int(float(f.readline()))
    perigee = EARTH_RADIUS * float(f.readline())
    radius = float(f.readline())
    speed = float(f.readline())
    spin = [float(x) for x in f.readline().split(',')]
    jlms = [float(x) for x in f.readline().split(',')]
    theta_true = [float(x) for x in f.readline().split(',')][:3]
    theta_high = np.asarray([float(x) for x in f.readline().split(',')])[:3]
    theta_low = np.asarray([float(x) for x in f.readline().split(',')])[:3]
    sigma = [float(d) for d in f.readline().split(", ")]# theta, ratio

data = np.array(asteroids_0_2.simulate(cadence, jlms, theta_true, radius,
                spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS, 0, False))

with open(f"{FILE_NAME}-0-samples.npy", 'rb') as f:
    samples = np.load(f).reshape(-1, 10)

def gen_data():
    results = []
    seeds = [theta_true]
    for _ in range(NUM_TRIALS):
        seeds.append((samples[random.randint(0, len(samples))][:3] - theta_true) * AMPLIFY+ theta_true)
    for i, seed in enumerate(seeds):
        print(i)
        try:
            results.append(asteroids_0_2.simulate(cadence, jlms, seed, radius,
                spin[0], spin[1], spin[2], perigee, speed, GM, EARTH_RADIUS, 0, False))
        except Exception:
            continue
    if SYMMETRIC:
        with open("data_runs-sym.npy", 'wb') as f:
            np.save(f, results)
    else:
        with open("data_runs-asym.npy", 'wb') as f:
            np.save(f, results)

def plot_sigma():
    if SYMMETRIC:
        with open("data_runs-sym.npy", 'rb') as f:
            runs = np.load(f)
    else:
        with open("data_runs-asym.npy", 'rb') as f:
            runs = np.load(f)
    time = np.arange(0, len(runs[0])//3) *  cadence / 3600.0
    time -= time[-1] / 2
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
plot_sigma()
plt.show()
