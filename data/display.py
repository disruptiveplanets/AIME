import sys, corner, emcee, os
import numpy as np
import matplotlib.pyplot as plt

EARTH_RADIUS = 6370000

class Display:
    def __init__(self, bare_name):
        self.bare_name = bare_name
        self.reader = emcee.backends.HDFBackend(self.bare_name+".h5", read_only=True)
        self.samples = None
        self.theta_true = None
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
            self.get_samples_burnin(0)
            self.burnin = 0
            self.thin = 1
            self.show_params()
            self.show_redchi()
            plt.show()
            self.burnin = -1
            while self.burnin < 0:
                try:
                    self.burnin = int(input("What is the burnin?"))
                except:
                    pass

        self.get_samples_burnin(self.burnin)

        print("burn-in: {0}".format(self.burnin))
        print("thin: {0}".format(self.thin))
        print("flat chain shape: {0}".format(self.samples.shape))
        print("flat log prob shape: {0}".format(self.log_prob_samples.shape))

    def get_samples_burnin(self, burnin):
        self.samples = self.reader.get_chain(discard=burnin, flat=True, thin=self.thin)
        self.log_prob_samples = self.reader.get_log_prob(discard=burnin, flat=True, thin=self.thin)
        self.log_prior_samples = self.reader.get_blobs(discard=burnin, flat=True, thin=self.thin)


    def show_params(self):
        self.get_samples()
        fig, axes = plt.subplots(self.ndim, figsize=(10, 7), sharex=True)
        samples = self.sampler.get_chain()
        for i in range(self.ndim):
            ax = axes[i]
            ax.plot(self.samples[:, :, i], "k", alpha=0.3)
            ax.set_xlim(0, len(self.samples))
            ax.set_ylabel(self.theta_labels[i])
            #ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number");
        plt.savefig(self.bare_name+"-params.png")

    def show_redchi(self):
        self.get_samples()
        fig = plt.figure(figsize=(10, 7),)
        for i in range(self.log_prob_samples.shape[1]):
            plt.plot(-self.log_prob_samples[:,i] / dof, c='k', alpha=0.25)
        plt.ylabel("Reduced chi squared")
        plt.xlabel("Sample")
        plt.ylim(0, 1e5)
        plt.savefig(self.bare_name+"-redchi.png")

    def show_corner(self):
        self.get_samples()
        self.get_params()
        fig = corner.corner(
            self.samples, labels=self.theta_labels, truths=self.theta_true
        );
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
        for i, (mean, plus, minus) in enumerate(res):
            print("{}: {} +{}, -{}. True: {}".format(
                self.theta_labels[i], mean, plus, minus, self.theta_true[i]))

    def get_results(self):
        self.get_samples()
        res = []
        for i in range(self.ndim):
            mcmc = np.percentile(self.samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            res.append([mcmc[1], q[0], q[1]])
        return res

    def run(self, theta):
        resolved_data = asteroids.simulate(cadence, jlms, theta[1:],
            spin[0], spin[1], spin[2], theta[0], impact_parameter, speed)
        return np.asarray(resolved_data)

    def show_compare(self):
        self.get_params()
        theta_results = [f[0] for f in self.get_results()]
        true_res = self.run(self.theta_true)
        mean_res = self.run(theta_results)

        plt.figure(figsize=(12, 4))

        x_display = np.arange(len(true_res) / 3)
        plt.plot(x_display, true_res[::3], label = 'true x', alpha=0.5, color='lightsalmon')
        plt.plot(x_display, true_res[1::3], label = 'true y', alpha=0.5, color='lightskyblue')
        plt.plot(x_display, true_res[2::3], label = 'true z', alpha=0.5, color='moccasin')

        plt.plot(x_display, mean_res[::3], label = 'mean x', alpha=0.5, linestyle='dotted', color='brown')
        plt.plot(x_display, mean_res[1::3], label = 'mean y', alpha=0.5, linestyle='dotted', color='steelblue')
        plt.plot(x_display, mean_res[2::3], label = 'mean z', alpha=0.5, linestyle='dotted', color='darkgoldenrod')
        plt.xlabel("Time (Cadences)")
        plt.ylabel("Spin (rad/s)")
        plt.legend()
        plt.savefig(self.bare_name+"-compare.png")
