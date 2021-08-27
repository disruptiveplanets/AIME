import sys, corner, emcee, os, asteroids
import numpy as np
import matplotlib.pyplot as plt

EARTH_RADIUS = 6370000
STANDARD_RESULTS_METHOD = False
REDCHI_THRESHOLD = 2
ASTEROIDS_MAX_K = 2

class Display:
    def __init__(self, bare_name):
        self.bare_name = bare_name
        self.reader = emcee.backends.HDFBackend(self.bare_name+".h5", read_only=True)
        self.samples = None
        self.theta_true = None
        self.true_results = None
        self.res = None
        try:
            self.chain_length, self.nwalkers, self.ndim = self.reader.get_chain().shape
        except AttributeError:
            raise Exception("Could not find the file {}.".format(self.bare_name+".h5"))
        self.theta_labels = list(map(r"$\theta_{{{0}}}$".format, range(1, self.ndim + 1)))

    def get_samples(self):
        if self.samples is not None:
            return
        try:
            tau = self.reader.get_autocorr_time()
            self.burnin = int(2 * np.max(tau))
            self.thin = int(0.5 * np.min(tau))
        except:
            print("Could not find autocorrelation time because the chain is too short.")
            self.thin = 1
            self.get_samples_burnin(0)
            self.burnin = 0
            self.show_params()
            plt.show()
            self.burnin = -1
            while self.burnin < 0:
                try:
                    self.burnin = int(input("What is the burnin? "))
                except:
                    pass

        self.get_samples_burnin(self.burnin)

        print("Burn-in: {}. Thin: {}".format(self.burnin, self.thin))

    def get_samples_burnin(self, burnin):
        self.samples = self.reader.get_chain(discard=burnin, thin=self.thin)
        self.log_prob_samples = self.reader.get_log_prob(discard=burnin, thin=self.thin)
        self.log_prior_samples = self.reader.get_blobs(discard=burnin, thin=self.thin)


    def show_params(self):
        self.get_samples()
        fig, axes = plt.subplots(self.ndim, figsize=(6.6, 4.6), sharex=True)
        samples = self.reader.get_chain()
        for i in range(self.ndim):
            ax = axes[i]
            ax.plot(self.samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(self.samples))
            ax.set_ylabel(self.theta_labels[i])
            #ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number");
        plt.savefig(self.bare_name+"-params.png")

    def show_redchi(self):
        self.get_params()
        self.get_samples()
        if self.true_results is None:
            self.true_results = self.run(self.theta_true)
        fig = plt.figure(figsize=(6.6, 4.6))
        redchiminmean = 0
        num_converged = 0
        for i in range(self.log_prob_samples.shape[1]):
            redchi = -self.log_prob_samples[:,i] / len(self.true_results)
            if redchi [-1] < REDCHI_THRESHOLD:
                num_converged += 1
            this_min = np.nanmin(redchi) / self.log_prob_samples.shape[1]
            redchiminmean += this_min if np.isfinite(this_min) else 0
            plt.plot(redchi, c='k', alpha=0.25)
        plt.ylabel("Reduced chi squared")
        plt.xlabel("Sample")
        plt.ylim(0, max(2, 2 * redchiminmean))
        plt.text(0.5, 0.5, "{} / {} walkers converged".format(num_converged, self.log_prob_samples.shape[1]),
        horizontalalignment='center', verticalalignment='center', transform = plt.gca().transAxes)
        plt.savefig(self.bare_name+"-redchi.png")

    def show_corner(self):
        self.get_samples()
        self.get_params()

        res = self.get_results()
        transpose_res = [np.asarray([l[i] for l in res]) for i in range(3)]

        flat_samples = self.samples.reshape(
            (self.samples.shape[0] * self.samples.shape[1], self.samples.shape[2]))
        fig = corner.corner(
            flat_samples, labels=self.theta_labels, truths=self.theta_true
        );

        corner.overplot_lines(fig, transpose_res[0], color='red')
        corner.overplot_lines(fig, transpose_res[0] + transpose_res[1], color='red', linestyle='dotted')
        corner.overplot_lines(fig, transpose_res[0] - transpose_res[2], color='red', linestyle='dotted')

        plt.savefig(self.bare_name+"-corner.png")

    def get_params(self):
        if self.theta_true is not None:
            return
        f = open(self.bare_name + ".dat", 'r')
        self.cadence = int(f.readline())
        self.impact_parameter = EARTH_RADIUS * int(f.readline())
        self.speed = float(f.readline())
        self.spin = [float(x) for x in f.readline().split(',')]
        self.jlms = [float(x) for x in f.readline().split(',')]
        self.theta_true = [float(x) for x in f.readline().split(',')]
        theta_start = [float(x) for x in f.readline().split(',')]
        theta_spread = [float(x) for x in f.readline().split(',')]
        theta_high = np.asarray([float(x) for x in f.readline().split(',')])
        theta_low = np.asarray([float(x) for x in f.readline().split(',')])
        sigma = float(f.readline()) * np.sqrt(self.spin[0]**2 + self.spin[1]**2 + self.spin[2]**2)
        f.close()

    def show_results(self):
        self.get_params()
        res = self.get_results()
        res_text = ""
        for i, (mean, plus, minus) in enumerate(res):
            res_text += "{}: {} +{}, -{}. True: {}\n".format(
                self.theta_labels[i], mean, plus, minus, self.theta_true[i])
        f = open(self.bare_name + '-results.txt', 'w')
        f.write(res_text)
        f.close()

    def get_results(self):
        if self.res is not None:
            return self.res
        self.get_samples()
        res = []

        if STANDARD_RESULTS_METHOD:
            mask = np.ones(self.log_prob_samples.shape[1], dtype=bool)
        else:
            # Only allow walkers which have converged. Definition of convergence: redchi
            self.get_params()
            if self.true_results is None:
                self.true_results = self.run(self.theta_true)

            redchi = -self.log_prob_samples[-1,:] / len(self.true_results)
            mask = redchi < REDCHI_THRESHOLD
            if np.sum(mask) == 0:
                mask = np.ones(self.log_prob_samples.shape[1], dtype=bool)

        print("Sampling results from {}/{} walkers".format(np.sum(mask), self.log_prob_samples.shape[1]))

        flat_samples = self.samples[:,mask,:].reshape(
            (self.samples.shape[0] * np.sum(mask), self.samples.shape[2]))

        for i in range(self.ndim):
            mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            res.append([mcmc[1], q[1], q[0]])# median, high bar, low bar

        self.res = res
        return res

    def run(self, theta):
        self.get_params()
        resolved_data = asteroids.simulate(self.cadence, self.jlms,
            theta[2:] if ASTEROIDS_MAX_K > 2 else theta[1:], # klm
            theta[2] if ASTEROIDS_MAX_K > 2 else 1, # radius
            self.spin[0], self.spin[1], self.spin[2], theta[0],
            self.impact_parameter, self.speed)
        return np.asarray(resolved_data)

    def show_compare(self):
        self.get_params()
        theta_results = [f[0] for f in self.get_results()]
        if self.true_results is None:
            self.true_results = self.run(self.theta_true)
        mean_res = self.run(theta_results)

        plt.figure(figsize=(12, 4))

        x_display = np.arange(len(self.true_results) / 3)
        plt.plot(x_display, self.true_results[::3], label = 'true x', alpha=0.5, color='C0')
        plt.plot(x_display, self.true_results[1::3], label = 'true y', alpha=0.5, color='C1')
        plt.plot(x_display, self.true_results[2::3], label = 'true z', alpha=0.5, color='C2')

        plt.plot(x_display, mean_res[::3], label = 'mean x', alpha=0.5, linestyle='dotted', color='C0')
        plt.plot(x_display, mean_res[1::3], label = 'mean y', alpha=0.5, linestyle='dotted', color='C1')
        plt.plot(x_display, mean_res[2::3], label = 'mean z', alpha=0.5, linestyle='dotted', color='C2')
        plt.xlabel("Time (Cadences)")
        plt.ylabel("Spin (rad/s)")
        plt.legend()
        plt.savefig(self.bare_name+"-compare.png")
