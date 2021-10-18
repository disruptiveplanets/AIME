from matplotlib import pyplot as plt
import numpy as np
import emcee

reader = emcee.backends.HDFBackend("run-9.0-0.h5", read_only=True)

try:
    tau = reader.get_autocorr_time()
    burnin = int(2 * np.max(tau))
    thin = int(0.5 * np.min(tau))
except:
    print("Could not find autocorrelation time because the chain is too short.")
    thin = 1
    burnin = 1000

samples = reader.get_chain(discard=burnin, thin=thin)
log_prob_samples = reader.get_log_prob(discard=burnin, thin=thin)
log_prior_samples = reader.get_blobs(discard=burnin, thin=thin)

flat_samples = np.array([samples[:,:,i].flatten() for i in range(3)])
for i in range(3):
    plt.figure()
    plt.hist(flat_samples[i], bins=22, fill=False)
plt.show()
