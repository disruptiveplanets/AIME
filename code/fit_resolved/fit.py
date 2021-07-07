import numpy as np
import matplotlib.pyplot as plt
import emcee, asteroids, time
from multiprocessing import Pool

EARTH_RADIUS = 6370000
EARTH_MASS =5.972e24
CADENCE = 30 * 60
REGENERATE_DATA = True
N_WALKERS = 32
N_STEPS = 5000

MULTIPROCESSING = True

np.random.seed(123)
L = 1
n = 1
m = 1
n_clms = (L + 1)**2
n_densities = n * n * (m + 1) * (2 * m + 1) // m

spin = 0.000050189
impact_parameter = 50 * EARTH_RADIUS
speed = 4000


theta_true = (
    5.7893, # C0
    0.89312, -0.412, -0.908901, # C1
    0.89124, 1.01298, 1.29021, 0.978512, 0.8912, 1.1250, # densities
)
theta_start = (
    5.0,
    0.0, 0.0, 0.0,
    1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
)
theta_range = (
    [4.0, 8.0],
    [-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0],
    [0.5, 1.5], [0.5, 1.5], [0.5, 1.5], [0.5, 1.5], [0.5, 1.5], [0.5, 1.5]
)
theta_labels = []
for l in range(L+1):
    for m in range(-l, l+1):
        theta_labels.append("C" + str(l) + str(m))
for i in range(n_densities):
    theta_labels.append("d" + str(i))

def fit_function(theta):
    clms = theta[:n_clms]
    densities = theta[n_clms:]
    start = time.time()
    resolved_data = asteroids.simulate(CADENCE, L, n, m, clms, densities, spin,
        impact_parameter, speed, EARTH_MASS)
    #print(time.time() - start)
    return resolved_data

def get_list_from_file(filename):
    f = open(filename, 'r')
    l = [float(line) for line in f.readlines()]
    f.close()
    return l
def save_to_file(filename, l):
    f = open(filename, 'w')
    for entry in l:
        f.write(str(entry) + '\n')
    f.close()

# Generate some synthetic data from the model.
if REGENERATE_DATA:
    start = time.time()
    y = fit_function(theta_true)
    print("Took {} s".format(time.time() - start))
    save_to_file("simulated-data.dat", y)
else:
    y = get_list_from_file("simulated-data.dat")
x = np.arange(len(y))
yerr = np.abs(0.1 * np.random.rand(len(y)) * y)
y += yerr * np.random.randn(len(y))

plt.figure(figsize=(12, 4))
x_display = np.arange(len(y) / 3)
plt.errorbar(x_display, y[::3], yerr=yerr[::3], label = 'x', fmt='.')
plt.errorbar(x_display, y[1::3], yerr=yerr[1::3], label = 'y', fmt='.')
plt.errorbar(x_display, y[2::3], yerr=yerr[2::3], label = 'z', fmt='.')
plt.xlabel("Time (Cadences)")
plt.ylabel("Spin (rad/s)")
plt.legend()
plt.show()



def log_likelihood(theta, y, yerr):
    # Normal likelihood
    model = fit_function(theta)
    sigma2 = yerr ** 2
    return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(sigma2))

def log_prior(theta):
    for i, param in enumerate(theta):
        if param > max(theta_range[i]) or param < min(theta_range[i]):
            return -np.inf
    return 0.0

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, y, yerr)

pos = theta_start + 1e-4 * np.random.randn(N_WALKERS, len(theta_start))
nwalkers, ndim = pos.shape

save_filename = "asteroids.h5"
backend = emcee.backends.HDFBackend(save_filename)
backend.reset(nwalkers, ndim)

if MULTIPROCESSING:
    with Pool() as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,
            args=(x, y, yerr), backend=backend, pool=pool)
        sampler.run_mcmc(pos, N_STEPS, progress=True)
else:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,
        args=(x, y, yerr), backend=backend)
    sampler.run_mcmc(pos, N_STEPS, progress=True)

fig, axes = plt.subplots(len(theta_true), figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(labels[i])
    #ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");

plt.show()

tau = sampler.get_autocorr_time()
print(tau)

flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
print(flat_samples.shape)


import corner

fig = corner.corner(
    flat_samples, labels=labels, truths=theta_true
);
plt.show()


inds = np.random.randint(len(flat_samples), size=100)
for ind in inds:
    sample = flat_samples[ind]
    plt.plot(x0, fit_function(x0, sample), "C1", alpha=0.1)
plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
plt.plot(x0, fit_function(x0, theta_true), "k", label="truth")
plt.legend(fontsize=14)
plt.xlim(0, 10)
plt.xlabel("x")
plt.ylabel("y");
plt.show()

for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    print("{}: {} + {} - {} in SI units".format(labels[i], mcmc[1], q[0], q[1]))
