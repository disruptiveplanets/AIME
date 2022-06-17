from mcmc_core import MCMCMethod, MIN_DENSITY, MAX_DENSITY, N_FITTED_MOMENTS, N_CONSTRAINED
import numpy as np
from scipy.special import lpmv, factorial

VERY_LARGE_SLOPE = 1e30
LMS = [
    # (l, m, is_real)
    (2, 2, True),
    (2, 0, True),
    (3, 3, True),
    (3, 3, False),
    (3, 2, True),
    (3, 2, False),
    (3, 1, True),
    (3, 1, False),
    (3, 0, True)]
VOLUME_SCALE = 4 / 3 * np.pi * (5/3)**(3/2)
RLM_EPSILON = 1e-20

# Made for spherical lumps

# Theta long: [lump_a, lump_mass | lump_a, lump_mass, lump_pos | ..., shell_mass, lump_pos]

MAX_LOG_PRIOR_LUMP = 1

def rlm(l,m,x,y,z):
    r = np.sqrt(np.maximum(RLM_EPSILON, x*x + y*y + z*z))
    return lpmv(m, l, z/r) / factorial(l + m) * r**l * np.exp(1j * m * np.arctan2(y, x))

class FiniteElement(MCMCMethod):
    # Parameters are theta_short

    def set_up(self, asteroid, generate):
        self.N = 1
        moment_field = self.asteroid.moment_field()
        bulk_am_sqr = np.sum(moment_field[-1] * asteroid.indicator_map) * DIVISION**3
        complex_klms = np.einsum("iabc,abc->i", moment_field[:-1,:,:,:], asteroid.indicator_map)
        self.shell_com = [
            complex_klms[3].real,
            complex_klms[3].imag,
            complex_klms[2].real
        ] ## CONFIRM
        self.shell_volume = np.sum(asteroid.indicator_map)
        self.surface_am = asteroid.surface_am
        self.shell_free_klms = np.array([
            complex_klms[4].real,
            complex_klms[6].real,
            complex_klms[15].real,
            complex_klms[15].imag,
            complex_klms[14].real,
            complex_klms[14].imag,
            complex_klms[13].real,
            complex_klms[13].imag,
            complex_klms[12].real,
        ])
        self.asteroid = asteroid
        assert(self.N >= 1) # Must have at least one lump
        assert(self.N <= MAX_LOG_PRIOR_LUMP) # Log prior can only handle so many

    def short_name(self):
        return "lumpy"

    def get_theta_long(self, theta_short):
        mass_sum = theta_short[1]
        for i in range(self.N):
            mass_sum += theta_short[3 + 5 * i]
        
        pos_sum = self.shell_com * (1 - mass_sum)
        for i in range(self.N - 1):
            pos_sum += theta_short[(4 + i * 5):(7 + i * 5)] * theta_short[3 + i * 5]
        return np.array(theta_short, [self.shell_volume - mass_sum, -pos_sum / theta_short[3]])

    def pick_parameters(self, local_rng):
        params = np.array([np.random.random() * self.surface_am, np.random.random()])
        for i in range(self.N - 1):
            params = np.append(params, [
                np.random.random() * self.surface_am,
                np.random.random() * self.shell_volume,
                np.random.random() * self.surface_am,
                np.random.random() * self.surface_am,
                np.random.random() * self.surface_am])

    def log_prior(self, theta_long):
        # What are the pairs that are overlapping?
        shell_density = theta_long[-4] / self.shell_volume
        lump_densities = [theta_long[1] / theta_long[0]**3 / VOLUME_SCALE]
        for i in range(0, self.N - 1):
            lump_densities.append(theta_long[3 + i * 5] / theta_long[2 + i * 5]**3 / VOLUME_SCALE)
        intersections = []

        # No lumps
        intersections.append(shell_density)
        # Single lumps
        for density in lump_densities:
            intersections.append(density + shell_density)

        return np.any(np.array(intersections) < MIN_DENSITY) or np.any(np.array(intersections) > MAX_DENSITY)

    def get_klms(self, theta_long):
        denom = self.surface_am ** 2 * theta_long[-4] # Shell moi
        denom += theta_long[1] * (theta_long[0]**2 + theta_long[-1]**2 + theta_long[-2]**2 + theta_long[-3]**2)# Lump 1 moi
        for i in range(self.N - 1):
            denom += theta_long[5 * i + 3] * (theta_long[5 * i + 2]**2 + theta_long[5 * i + 4]**2 + theta_long[5 * i + 5]**2 + theta_long[5 * i + 6]**2)
        # Surface
        klms = theta_long[-4] * self.shell_free_klms
        # Lumps
        for i in range(self.N):
            lump_mass = theta_long[1] if i == 0 else theta_long[3 + 5 * i]
            lump_pos = theta_long[-3:] if i == 0 else theta_long[(4 + i * 5):(7 + i * 5)]
            for j, (l, m, is_real) in enumerate(LMS):
                rlm_val = rlm(l, m, lump_pos)
                klms[j] += lump_mass / self.surface_am**l * (rlm_val.real if is_real else rlm_val.imag)

        return klms / denom * self.surface_am**2

    def get_map(self, means, unc, asteroid):
        densities = np.ones_like(self.asteroid.indicator_map) * means[-4] / self.shell_volume
        X, Y, Z = np.meshgrid(self.asteroid.gridline, self.asteroid.gridline, self.asteroid.gridline)

        for i in range(self.N):
            lump_pos = means[-3:] if i == 0 else means[(4 + i * 5):(7 + i * 5)]
            lump_length = means[0] if i == 0 else means[2 + 5 * i]
            lump_mass = means[1] if i == 0 else means[3 + 5 * i]
            lump_density = lump_mass / lump_length / VOLUME_SCALE
            densities[(X - lump_pos[0]) ** 2 + (Y - lump_pos[1]) ** 2 + (Z - lump_pos[2]) ** 2 < (lump_length**2 * 5 / 3)] += lump_density

        densities[np.isnan(self.asteroid.indicator_map)] = np.nan
        return densities