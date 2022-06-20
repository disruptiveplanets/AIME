from mcmc_core import MCMCMethod, MIN_DENSITY, MAX_DENSITY, N_FITTED_MOMENTS, N_CONSTRAINED
import numpy as np
import grids
from scipy.linalg import inv

VERY_LARGE_SLOPE = 1e30


class FiniteElement(MCMCMethod):
    # Parameters are theta_short

    def set_up(self, asteroid, generate):
        if generate:
            masks = grids.get_grids_centroid(self.n_all, asteroid.grid_line, asteroid.indicator_map, asteroid.indicator)[1]
            with open(asteroid.name + "-grids.npy", 'wb') as f:
                np.save(f, masks)
        else:
            with open(asteroid.name + "-grids.npy", 'rb') as f:
                masks = np.load(f)

        rlms = asteroid.moment_field(surface_am=asteroid.surface_am)

        rlm_mat_complex = np.einsum("iabc,jabc->ji", masks, rlms) * asteroid.division**3

        # Set order
        rlm_mat = np.zeros((N_FITTED_MOMENTS + N_CONSTRAINED, self.n_all))
        # Unconstrained
        rlm_mat[0, :] = rlm_mat_complex[8,:].real # K22
        rlm_mat[1, :] = rlm_mat_complex[6,:].real # K20
        rlm_mat[2, :] = rlm_mat_complex[15,:].real # R K33
        rlm_mat[3, :] = rlm_mat_complex[15,:].imag # I K33
        rlm_mat[4, :] = rlm_mat_complex[14,:].real # R K32
        rlm_mat[5, :] = rlm_mat_complex[14,:].imag # I K32
        rlm_mat[6, :] = rlm_mat_complex[13,:].real # R K31
        rlm_mat[7, :] = rlm_mat_complex[13,:].imag # I K31
        rlm_mat[8, :] = rlm_mat_complex[12,:].real # K30
        # Constrained
        rlm_mat[9, :] = rlm_mat_complex[ 3,:].real # R K11
        rlm_mat[10, :] = rlm_mat_complex[3,:].imag # I K11
        rlm_mat[11, :] = rlm_mat_complex[2,:].real # K10
        rlm_mat[12, :] = rlm_mat_complex[7,:].real # R K21
        rlm_mat[13, :] = rlm_mat_complex[7,:].imag # I K21
        rlm_mat[14, :] = rlm_mat_complex[8,:].imag # I K22
        rlm_mat[15, :] = rlm_mat_complex[0,:].real # K00
        radius_vec = rlm_mat_complex[-1, :].real # radius

        rlm_fixed_inv = inv(rlm_mat[np.arange(self.n_free, self.n_all),:][:,np.arange(self.n_free, self.n_all)])
        rlm_cross = rlm_mat[np.arange(self.n_free, self.n_all),:][:,np.arange(0, self.n_free)]
        rlm_prod = rlm_fixed_inv @ rlm_cross

        self.rlm_mat = rlm_mat
        self.radius_vec = radius_vec
        self.rlm_fixed_inv = rlm_fixed_inv
        self.rlm_prod = rlm_prod
        self.masks = masks

    def short_name(self):
        return "fe"

    def get_theta_long(self, theta_short):
        # Get the densities consistent with making the mass 1 and com 0 and rotation
        return np.append(theta_short, self.rlm_fixed_inv[:,-1] - self.rlm_prod @ theta_short)

    def pick_parameters(self, local_rng):
        return [(local_rng.random() * 4 + 0.5) * self.mean_density for _ in range(self.n_free)]

    def log_prior(self, theta_long):
        mask_too_small = theta_long < self.mean_density * MIN_DENSITY
        if np.any(mask_too_small):
            return VERY_LARGE_SLOPE * np.sum((theta_long / self.mean_density - MIN_DENSITY) * mask_too_small)
        mask_too_big = theta_long > self.mean_density * MAX_DENSITY
        if np.any(mask_too_big):
            return -VERY_LARGE_SLOPE * np.sum((np.max(theta_long) / self.mean_density - MAX_DENSITY) * mask_too_big)
        return 0.0

    def get_klms(self, theta_long):
        unscaled_klms = self.rlm_mat @ theta_long
        radius_sqr = self.radius_vec @ theta_long
        scaled_klms = np.array(unscaled_klms) / radius_sqr
        scaled_klms[-1] *= radius_sqr # Do not scale mass term
        # Mass is already normalized
        return scaled_klms

    def get_map(self, means, unc, asteroid):
        densities = np.einsum("ijkl,i->jkl", self.masks, means)
        unc_ratios = np.einsum("ijkl,i->jkl", self.masks, unc)
        densities = densities / np.nanmean(densities)

        densities[~asteroid.indicator_map] = np.nan
        unc_ratios[~asteroid.indicator_map] = np.nan

        return densities, unc_ratios