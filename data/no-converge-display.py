import sys, corner, emcee, os
import numpy as np
import matplotlib.pyplot as plt

theta_true = (
    0, 1.2e6, 1.1e5, -4.9e5,
)

dof = 1419

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

burnin = 0
thin = 1
theta_labels = list(map(r"$\theta_{{{0}}}$".format, range(1, ndim + 1)))
samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
log_prob_samples = reader.get_log_prob(discard=burnin, thin=thin)
log_prior_samples = reader.get_blobs(discard=burnin, flat=True, thin=thin)
print(log_prob_samples.shape)

print("burn-in: {0}".format(burnin))
print("thin: {0}".format(thin))
print("flat chain shape: {0}".format(samples.shape))
print("flat log prob shape: {0}".format(log_prob_samples.shape))
#print("flat log prior shape: {0}".format(log_prior_samples.shape))
fig = plt.figure(figsize=(10, 7),)
for i in range(log_prob_samples.shape[1]):
    plt.plot(-log_prob_samples[:,1] / dof, c='k', alpha=0.25)
plt.ylabel("Reduced chi squared")
plt.xlabel("Sample")
plt.ylim(0, 1e5)
plt.show()


burnin = -1
while burnin < 0:
    try:
        burnin = int(input("What is the burnin?"))
    except:
        pass

flat_samples = reader.get_chain(discard=burnin, thin=thin, flat=True)
print(flat_samples.shape)

fig = corner.corner(
    flat_samples, labels=theta_labels, truths=theta_true
);
plt.savefig(bare_name+".png")
plt.show()

for i in range(ndim):
    mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
    q = np.diff(mcmc)
    print("{}: {} + {} - {}. Compare {}".format(theta_labels[i], mcmc[1], q[0], q[1], theta_true[i]))
