TEST = False

SAVE_PDFS = False

import corner, emcee
if not TEST:
    import asteroids_0_3, asteroids_0_2, asteroids_2_3, asteroids_2_2, asteroids_3_3, asteroids_3_2
else:
    import test_loglike as asteroids_0_2
    import test_loglike as asteroids_0_3
    import test_loglike as asteroids_2_2
    import test_loglike as asteroids_2_3
    import test_loglike as asteroids_3_2
    import test_loglike as asteroids_3_3
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import pinvh

import matplotlib as mpl
mpl.rcParams["text.usetex"] = True
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["legend.framealpha"] = 0.5
mpl.rcParams["lines.linewidth"] = 2
mpl.rcParams["lines.markersize"] = 4
mpl.rcParams["font.size"] = 20
mpl.rcParams["legend.fontsize"] = 10


EARTH_RADIUS = 6370000
STANDARD_RESULTS_METHOD = False
REDCHI_THRESHOLD = 2
CENTRAL_RADIUS = EARTH_RADIUS
CENTRAL_MU = 5.972e24 * 6.674e-11

DATA_POINT_LIMIT = 642


AUTO_BURNIN = 100
DEFAULT_THIN = 10

class Display:
    def __init__(self, bare_name, h5_name):
        self.bare_name = bare_name
        self.h5_name = h5_name
        self.reader = emcee.backends.HDFBackend(self.h5_name+".h5", read_only=True)
        self.samples = None
        self.theta_true = None
        self.true_results = None
        self.res = None
        self.mask = None
        self.maxj = None
        self.maxk = None
        self.module = None
        try:
            self.chain_length, self.nwalkers, self.ndim = self.reader.get_chain().shape
        except AttributeError:
            raise Exception("Could not find the file {}.".format(self.h5_name+".h5"))
        #self.theta_labels = list(map(r"$\Delta\theta_{{{0}}}$".format, range(1, self.ndim + 1)))
        self.theta_labels = ["$\Delta\\gamma_0$", "$\Delta K_{22}$", "$\Delta K_{20}$"]
        if self.ndim >= 10:
            self.theta_labels += ["$\Re\Delta K_{33}$", "$\Im\Delta K_{33}$",
                                  "$\Re\Delta K_{32}$", "$\Im\Delta K_{32}$",
                                  "$\Re\Delta K_{31}$", "$\Im\Delta K_{31}$",
                                  "$\Delta K_{30}$"]

    def get_true_results(self):
        if self.true_results is not None:
            return
        self.true_results = np.load(f"{self.bare_name}-data.npy")
        self.true_uncs = np.load(f"{self.bare_name}-unc.npy")

    def get_samples(self):
        if self.samples is not None:
            return
        try:
            tau = self.reader.get_autocorr_time()
            self.burnin = int(2 * np.max(tau))
            self.thin = DEFAULT_THIN#int(0.5 * np.min(tau))
        except:
            print("Could not find autocorrelation time because the chain is too short.")
            self.burnin = AUTO_BURNIN
            self.thin = DEFAULT_THIN

        self.samples = self.reader.get_chain(discard=self.burnin, thin=self.thin)
        self.log_prob_samples = self.reader.get_log_prob(discard=self.burnin, thin=self.thin)
        self.log_prior_samples = self.reader.get_blobs(discard=self.burnin, thin=self.thin)

        print("Burn-in: {}. Thin: {}".format(self.burnin, self.thin))

    def show_params(self):
        self.get_samples()
        fig, axes = plt.subplots(self.ndim, figsize=(6.6, 6.6), sharex=True)
        for i in range(self.ndim):
            ax = axes[i]
            ax.plot(self.samples[:, :, i] - self.theta_true[i], "k", alpha=0.3)
            ax.set_xlim(0, len(self.samples))
            ax.set_ylabel(self.theta_labels[i])
            #ax.yaxis.set_label_coords(-0.1, 0.5)
        axes[-1].set_xlabel("step number");
        
        plt.savefig(self.h5_name+"-params.png")

    def show_redchi(self):
        self.get_params()
        self.get_samples()
        if self.true_results is None:
            self.get_true_results()
        fig = plt.figure(figsize=(6.6, 4.6))
        redchiminmean = 0
        num_converged = 0

        f = open(self.h5_name + '-redchis.txt', 'w')
        for i in range(self.log_prob_samples.shape[1]):
            redchi_list = -self.log_prob_samples[:,i] / len(self.true_results) / 3 * 2
            redchi = np.nanmin(redchi_list)
            f.write(str(redchi) + "\n")
            if redchi < REDCHI_THRESHOLD:
                num_converged += 1
            this_min = redchi / self.log_prob_samples.shape[1]
            redchiminmean += this_min if np.isfinite(this_min) else 0
            plt.plot(redchi_list, c='k', alpha=0.25)
        f.close()

        plt.ylabel("Reduced chi squared")
        plt.xlabel("Sample")
        plt.ylim(0, max(2, 2 * redchiminmean))
        plt.text(0.5, 0.5, "{} / {} walkers converged".format(num_converged, self.log_prob_samples.shape[1]),
        horizontalalignment='center', verticalalignment='center', transform = plt.gca().transAxes)
        plt.tight_layout()
        plt.savefig(self.h5_name+"-redchi.png")


    def get_mask(self):
        if self.mask is not None:
            return self.mask
        self.get_samples()

        if STANDARD_RESULTS_METHOD:
            mask =  np.ones(self.log_prob_samples.shape[1], dtype=bool)
        else:
            # Only allow walkers which have converged. Definition of convergence: redchi
            self.get_params()
            if self.true_results is None:
                self.get_true_results()

            redchi = -self.log_prob_samples[-1,:] / len(self.true_results) * 2
            mask = redchi < REDCHI_THRESHOLD
            if np.sum(mask) == 0:
                mask = np.ones(self.log_prob_samples.shape[1], dtype=bool)

        print("Sampling results from {}/{} walkers".format(np.sum(mask), self.log_prob_samples.shape[1]))
        self.mask = mask
        return mask

    def show_corner(self):
        self.get_samples()
        self.get_params()
        self.get_mask()

        flat_samples = self.samples[:,self.mask,:].reshape(
            (-1, self.samples.shape[2])) - np.array(self.theta_true)

        param_exps = np.ceil(-np.log10(np.maximum(np.abs(np.min(flat_samples, axis=0)), np.max(flat_samples, axis=0))))
        exp_labels = []
        for l, e in zip(self.theta_labels, param_exps):
            exp_labels.append(l + " ($\\times 10^{"+str(int(-e))+"}$)")

        corner.corner(
            flat_samples * 10**param_exps, labels=exp_labels, truths=np.zeros_like(self.theta_true)
        );

        #corner.overplot_lines(fig, (transpose_res[0]) * 10**param_exps, color='C1', linestyle='dashed')
        #corner.overplot_lines(fig, (transpose_res[0] + transpose_res[1]) * 10**param_exps, color='C1', linestyle='dotted')
        #corner.overplot_lines(fig, (transpose_res[0] - transpose_res[2]) * 10**param_exps, color='C1', linestyle='dotted')

        plt.savefig(self.h5_name+"-corner.png")
        if SAVE_PDFS:
            plt.savefig(self.h5_name+"-corner.pdf")
        return len(self.true_results)

    def get_params(self):
        if self.theta_true is not None:
            return
        f = open("../../staged/" + self.bare_name + ".txt", 'r')
        asteroids_max_j, asteroids_max_k = f.readline().split(', ')
        self.maxj = int(asteroids_max_j)
        self.maxk = int(asteroids_max_k)
        self.cadence = int(float(f.readline()))
        self.impact_parameter = EARTH_RADIUS * float(f.readline())
        self.radius = float(f.readline())
        self.speed = float(f.readline())
        self.spin = [float(x) for x in f.readline().split(',')]
        self.jlms = [float(x) for x in f.readline().split(',')]
        self.theta_true = [float(x) for x in f.readline().split(',')]
        theta_high = np.asarray([float(x) for x in f.readline().split(',')])
        theta_low = np.asarray([float(x) for x in f.readline().split(',')])
        self.sigma = [float(d) for d in f.readline().split(',')]
        last_line = f.readline()
        self.velocity_mul = 1 if last_line == '' else float(last_line)


        try:
            self.gap = float(f.readline())
        except Exception:
            self.gap = None

        f.close()

        if self.maxj == 0:
            if self.maxk == 2:
                self.module = asteroids_0_2
            elif self.maxk == 3:
                self.module = asteroids_0_3

        elif self.maxj == 2:
            if self.maxk == 2:
                self.module = asteroids_2_2
            elif self.maxk == 3:
                self.module = asteroids_2_3

        elif self.maxj == 3:
            if self.maxk == 2:
                self.module = asteroids_3_2
            elif self.maxk == 3:
                self.module = asteroids_3_3
        else:
            raise Exception("Could not find module for simulation.")

    def show_results(self):
        self.get_params()
        res = self.get_results()
        res_text = ""
        for i, (mean, plus, minus) in enumerate(res):
            res_text += "{}: {} +{}, -{}. True: {}\n".format(
                self.theta_labels[i], mean, plus, minus, self.theta_true[i])
        f = open(self.h5_name + '-results.txt', 'w')
        f.write(res_text)
        f.close()

    def get_results(self):
        if self.res is not None:
            return self.res
        self.get_samples()
        self.get_mask()
        res = []

        flat_samples = self.samples[:,self.mask,:].reshape(
            (self.samples.shape[0] * np.sum(self.mask), self.samples.shape[2]))

        for i in range(self.ndim):
            mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
            q = np.diff(mcmc)
            res.append([mcmc[1], q[1], q[0]])# median, high bar, low bar

        self.res = res
        return res

    def run(self, theta):
        self.get_params()
        try:
            resolved_data = self.module.simulate(self.cadence, self.jlms, theta,
                self.radius, self.spin[0], self.spin[1], self.spin[2],
                self.impact_parameter, self.speed, CENTRAL_MU, CENTRAL_RADIUS, 0, False, self.velocity_mul)
        except:
            print("Coordinates are invalid ({})".format(theta))
            return None
        vector_data = np.asarray(resolved_data).reshape(-1, 3)

        if self.gap is not None:
            # Time of data in hours
            times = np.arange(len(vector_data)) * self.cadence / 3600
            times -= times[-1] / 2

            # Generate mask
            uncut_mask = abs(times) > self.gap
            chop_index = 0
            while chop_index < len(uncut_mask)/2:
                if np.sum(uncut_mask[chop_index:-(1+chop_index)]) > DATA_POINT_LIMIT:
                    chop_index += 1
                else:
                    break
            uncut_mask[:chop_index] = False
            uncut_mask[-(1+chop_index):] = False
            
            return vector_data[uncut_mask]
        else:
            return vector_data

    def show_compare(self):
        self.get_params()
        theta_results = [f[0] for f in self.get_results()]
        if self.true_results is None:
            self.true_results = self.get_true_results()
        mean_res = list(self.run(theta_results))
        
        while len(mean_res) > len(self.true_results):
            del mean_res[-1]

        while len(mean_res) < len(self.true_results):
            mean_res.append(mean_res[-1])

        mean_res = np.array(mean_res)
        uncertainties = np.array([2 * np.sqrt(np.diagonal(pinvh(a))) for a in self.true_uncs])

        fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(9, 6), sharex=True)

        x_display = np.arange(len(self.true_results)) * self.cadence / 3600.0
        x_display -= np.median(x_display)
        ax1.scatter(x_display, self.true_results[:,0] * 3600, label = 'true x', alpha=0.5, color='C0', s=1)
        ax1.scatter(x_display, self.true_results[:,1] * 3600, label = 'true y', alpha=0.5, color='C1', s=1)
        ax1.scatter(x_display, self.true_results[:,2] * 3600, label = 'true z', alpha=0.5, color='C2', s=1)

        ax1.fill_between(x_display, self.true_results[:,0] * 3600 + uncertainties[:,0] * 3600,
                self.true_results[:,0] * 3600 - uncertainties[:,0] * 3600, color="C0", alpha=0.2)
        ax1.fill_between(x_display, self.true_results[:,1] * 3600 + uncertainties[:,1] * 3600,
                self.true_results[:,1] * 3600 - uncertainties[:,1], color="C1", alpha=0.2)
        ax1.fill_between(x_display, self.true_results[:,2] * 3600 + uncertainties[:,2] * 3600,
                self.true_results[:,2] * 3600 - uncertainties[:,2] * 3600, color="C2", alpha=0.2)

        if mean_res is not None:
            ax1.plot(x_display, mean_res[:,0] * 3600, label = 'mean x', color='C0')
            ax1.plot(x_display, mean_res[:,1] * 3600, label = 'mean y', color='C1')
            ax1.plot(x_display, mean_res[:,2] * 3600, label = 'mean z', color='C2')
        ax1.set_ylabel("Spin (rad/hr)")
        ax1.legend()

        if mean_res is not None:
            ax2.scatter(x_display, mean_res[:,0] * 3600 - self.true_results[:,0] * 3600, color='C0', s=1)
            ax2.scatter(x_display, mean_res[:,1] * 3600 - self.true_results[:,1] * 3600, color='C1', s=1)
            ax2.scatter(x_display, mean_res[:,2] * 3600 - self.true_results[:,2] * 3600, color='C2', s=1)
            ax2.fill_between(x_display, uncertainties[:,0] * 3600, -uncertainties[:,0] * 3600, color="C0", alpha=0.2)
            ax2.fill_between(x_display, uncertainties[:,1] * 3600, -uncertainties[:,1] * 3600, color="C1", alpha=0.2)
            ax2.fill_between(x_display, uncertainties[:,2] * 3600, -uncertainties[:,2] * 3600, color="C2", alpha=0.2)

        ax2.set_ylabel("Residuals (rad/hr)")
        ax2.set_xlabel("Time to perigee (hours)")
        ax2.set_xlim(np.min(x_display), np.max(x_display))
        plt.tight_layout()

        plt.savefig(self.h5_name+"-compare.png")
        if SAVE_PDFS:
            plt.savefig(self.h5_name+"-compare.pdf")


if __name__ == "__main__":
    d = Display("minimizer/run-3.0/run-3.0")
    d.thin = DEFAULT_THIN
    d.get_samples_burnin(4000)
    d.show_compare()
    plt.show()
