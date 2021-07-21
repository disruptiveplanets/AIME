import sys, corner, emcee, os
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) != 2:
    print("Please provide an h5 file to display")
    sys.exit()

filename = sys.argv[1]
bare_name = os.path.splitext(filename)[0]
reader = emcee.backends.HDFBackend(filename, read_only=True)

try:
    chain_length, nwalkers, ndim = reader.get_chain().shape
except AttributeError:
    print("Could not find the file {}.".format(filename))
    sys.exit()

#try:

tau = reader.get_autocorr_time()
burnin = int(2 * np.max(tau))
thin = int(0.5 * np.min(tau))
samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)
log_prior_samples = reader.get_blobs(discard=burnin, flat=True, thin=thin)

print("burn-in: {0}".format(burnin))
print("thin: {0}".format(thin))
print("flat chain shape: {0}".format(samples.shape))
print("flat log prob shape: {0}".format(log_prob_samples.shape))
print("flat log prior shape: {0}".format(log_prior_samples.shape))

'''except AutocorrError:
    print("Could not find autocorrelation time because the chain is too short.")
    samples = reader.get_chain(flat=True)
    log_prob_samples = reader.get_log_prob(flat=True)
    log_prior_samples = reader.get_blobs(flat=True)
    print(log_prob_samples)
    print(log_prior_samples)
    sys.exit()'''

all_samples = np.concatenate(
    (samples, log_prob_samples[:, None], log_prior_samples[:, None]), axis=1
)

labels = list(map(r"$\theta_{{{0}}}$".format, range(1, ndim + 1)))
labels += ["log prob", "log prior"]

corner.corner(all_samples, labels=labels);

plt.savefig(bare_name + ".png")
plt.show()













fig, axes = plt.subplots(len(theta_true), figsize=(10, 7), sharex=True)
samples = sampler.get_chain()
for i in range(ndim):
    ax = axes[i]
    ax.plot(samples[:, :, i], "k", alpha=0.3)
    ax.set_xlim(0, len(samples))
    ax.set_ylabel(theta_labels[i])
    #ax.yaxis.set_label_coords(-0.1, 0.5)

axes[-1].set_xlabel("step number");

plt.show()

tau = sampler.get_autocorr_time()
print(tau)

flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
print(flat_samples.shape)

fig = corner.corner(
    flat_samples, labels=theta_labels, truths=theta_true
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
    print("{}: {} + {} - {} in SI units".format(theta_labels[i], mcmc[1], q[0], q[1]))
