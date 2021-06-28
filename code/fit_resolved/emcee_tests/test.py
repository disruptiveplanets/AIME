import numpy as np
import matplotlib.pyplot as plt

np.random.seed(123)

# Choose the "true" parameters.
theta_true = (1.2412, 4.294, 0.534, 0.412)

theta_range = ((1.5, 0.5), (1.0, 5.0), (0.0, 1.0), (0.0, 1.0))
theta_start = theta_true
labels = ["a", "b", "c", "d"]

def fit_function(x, theta):
    a, b, c, d = theta
    return a * x * x + b * x + c + d

# Generate some synthetic data from the model.
N = 50
x = np.sort(10 * np.random.rand(N))
y = fit_function(x, theta_true)
yerr = np.abs(0.1 * np.random.rand(N) * y)
y += yerr * np.random.randn(N)

plt.errorbar(x, y, yerr=yerr, fmt=".k", capsize=0)
x0 = np.linspace(0, 10, 500)
plt.plot(x0, fit_function(x0, theta_true), "k")
plt.xlim(0, 10)
plt.xlabel("x")
plt.ylabel("y");

plt.show()

def log_likelihood(theta, x, y, yerr):
    # Normal likelihood
    model = fit_function(x, theta)
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
    return lp + log_likelihood(theta, x, y, yerr)

import emcee

pos = theta_start + 1e-4 * np.random.randn(32, len(theta_start))
nwalkers, ndim = pos.shape

sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability, args=(x, y, yerr))
sampler.run_mcmc(pos, 5000, progress=True);

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
