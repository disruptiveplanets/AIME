import numpy as np
from scipy.special import lpmv, factorial
from scipy.linalg import norm
from display import make_gif, make_slices
import warnings, os

def rlm(l,m,x,y,z):
    r = np.sqrt(np.maximum(0, x*x + y*y + z*z))
    return lpmv(m, l, z / r) / factorial(l + m) * r**l * np.exp(1j * m * np.arctan2(y, x))

def rlm_gen(l,m):
    return lambda x,y,z: rlm(l,m,x,y,z)

class Method:
    def __init__(self, asteroid):
        self.asteroid = asteroid
        self.unc = None
        self.d = None
        self.densities = None
        self.density_uncs = None

    def get_a(self):
        raise NotImplementedError()
        
    def get_b(self, pos):
        raise NotImplementedError()

    def solve(self):
        if self.unc is None:
            a = self.get_a()
            self.d = a @ self.asteroid.data
            self.unc = a @ self.asteroid.sigma_data @ a.transpose()


    def get_density(self, x,y,z):
        out = (self.get_b(x,y,z) @ self.d).real
        if type(out) != np.float64:
            return out[0].real
        return out.real

    def get_unc(self, x,y,z):
        b = self.get_b(x,y,z)
        if len(x.shape) > 0:
            b = b.reshape(-1, b.shape[-1])
            return np.array([np.sqrt(e @ self.unc @ e)[0] for e in b]).reshape(x.shape).real
        else:
            out = np.sqrt(b @ self.unc @ b.transpose()).real
            if type(out) != np.float64:
                return out[0].real
            return out.real

    def map_density(self):
        if self.densities is None:
            self.densities = self.asteroid.map(self.get_density, dtype=float)
        return self.densities

    def map_unc(self):
        if self.density_uncs is None:
            self.density_uncs = self.asteroid.map(self.get_unc, dtype=float)
        return self.density_uncs

    def save_density(self, fname):
        with open(fname, 'wb') as f:
            np.save(f, self.map_density())

    def save_unc(self, fname):
        with open(fname, 'wb') as f:
            np.save(f, self.map_unc())

    def check(self):
        rlms_field = self.asteroid.moment_field()
        index = 0
        for l in range(0, self.asteroid.max_l + 1):
            for m in range(-l, l+1):
                rlms_field[index] /= self.asteroid.am ** l
                index += 1
        rlms_field[index] /= self.asteroid.am**2
        rlms = np.sum(rlms_field * self.map_density() * self.asteroid.division**3, axis=(1,2,3))
        for i, r in enumerate(rlms):
            print("Expected  {:.5f}\t Got {:.5f} \t Difference {:.2g}".format(self.asteroid.data[i], r, abs(self.asteroid.data[i]-r)))

    def display(self, fname, duration=2):
        asteroid_name = fname.split("/")[1]
        if not os.path.isdir(f"figs/{asteroid_name}"):
            os.mkdir(f"figs/{asteroid_name}")

        warnings.filterwarnings("ignore")

        display_densities = np.copy(self.densities)
        display_densities[~self.asteroid.indicator_map] = np.nan
        display_uncs = np.copy(self.density_uncs)
        display_uncs[~self.asteroid.indicator_map] = np.nan

        display_uncs /= np.abs(display_densities)
        display_densities /= np.nanmean(display_densities)

        make_slices(display_densities, self.asteroid.grid_line, "$\\rho$", 'plasma', f"{fname}-d")
        make_gif(display_densities, self.asteroid.grid_line, "$\\rho$", 'plasma', f"{fname}-d.gif", duration)
        make_slices(display_uncs, self.asteroid.grid_line, "$\\sigma_\\rho$", 'Greys_r', f"{fname}-u")
        make_gif(display_uncs, self.asteroid.grid_line, "$\\sigma_\\rho$", 'Greys_r', f"{fname}-u.gif", duration)

        warnings.filterwarnings("default")

    def reload(self, density_name, unc_name):
        with open(density_name, 'rb') as f:
            self.densities = np.load(f)
        with open(unc_name, 'rb') as f:
            self.density_uncs = np.load(f)


class Asteroid:
    def __init__(self, sample_file, am, division, max_radius, indicator):
        self.max_l = None
        self.am = am
        self.division = division
        self.grid_line = np.arange(-max_radius, max_radius, division)
        self.indicator = indicator
        self.indicator_map = self.get_indicator_map()
        if sample_file != "":
            self.data, self.sigma_data = self.load_samples(sample_file)
        self.moments = None

    def load_samples(self, fname):
        with open(fname, 'rb') as f:
            samples = np.load(f)
        samples = samples.reshape(-1, samples.shape[-1]).transpose()
        if samples.shape[0] not in [3, 10]:
            raise Exception("Only l=2 or l=3 are supported")
        if samples.shape[0] == 3:
            self.max_l = 2
        elif samples.shape[0] == 10:
            self.max_l = 3
        klms = np.zeros((samples.shape[0]+6, samples.shape[1]), dtype=complex)
        #             00
        #         1-1 10 11
        #     2-2 2-1 20 21 22
        # 3-3 3-2 3-1 30 31 32 33
        klms[0] = 1
        klms[1] = 0
        klms[2] = 0
        klms[3] = 0
        klms[4] = samples[1]
        klms[5] = 0
        klms[6] = samples[2]
        klms[7] = 0
        klms[8] = samples[1]
        ms = np.array([0, -1, 0, 1, -2, -1, 0, 1, 2])
        if self.max_l > 2:
            klms[9] = -samples[3] + 1j * samples[4]
            klms[10] = samples[5] - 1j * samples[6]
            klms[11] = -samples[7] + 1j * samples[8]
            klms[12] = samples[9]
            klms[13] = klms[11].conj()
            klms[14] = klms[10].conj()
            klms[15] = klms[9].conj()
            ms = np.append(ms, [-3, -2, -1, 0, 1, 2, 3])

        delta_gamma = samples[0] - np.mean(samples[0])
        hybrids = klms * np.exp(-1j * np.outer(ms, delta_gamma))
        klm_means = np.mean(hybrids, axis=1)
        klm_cov = np.cov(hybrids)
        added_one_cov = np.hstack((klm_cov, np.zeros((klm_cov.shape[0], 1))))
        added_cov = np.vstack((added_one_cov, np.zeros((1, klm_cov.shape[0]+1))))

        return np.append(klm_means, 1), added_cov


    def get_indicator_map(self):
        x,y,z = np.meshgrid(self.grid_line, self.grid_line, self.grid_line)
        return self.indicator(x,y,z)


    def map(self, fn, dtype=np.complex):
        result = np.zeros((len(self.grid_line), len(self.grid_line), len(self.grid_line)), dtype=dtype)
        for nx, x in enumerate(self.grid_line):
            for ny, y in enumerate(self.grid_line):
                for nz, z in enumerate(self.grid_line):
                    if self.indicator_map[nx,ny,nz]:
                        result[nx, ny, nz] = fn(x,y,z)
        return result

    def map_np(self, fn):
        x,y,z = np.meshgrid(self.grid_line, self.grid_line, self.grid_line)
        return fn(x,y,z) * self.indicator(x,y,z)


    def columnate(self, a):
        # Return a, but formatted as a column where indicator is 1
        return a.reshape(-1)[self.indicator_map.reshape(-1)]


    def moment_field(self, max_l=None):
        if max_l is None:
            max_l = self.max_l
        if self.moments is None:
            n = (max_l+1)**2 + 1
            self.moments = np.zeros((n, len(self.grid_line), len(self.grid_line), len(self.grid_line)), dtype=np.complex)
            i = 0
            for l in range(0, max_l+1):
                for m in range(-l, l+1):
                    self.moments[i] = self.map_np(rlm_gen(l,m))
                    i += 1
            self.moments[i] = self.map_np(lambda x,y,z: x*x + y*y + z*z)
        return self.moments


class Indicator:
    def sph(am):
        return lambda x,y,z: x*x + y*y + z*z < am*am*5/3

    def ell(am, k22, k20):
        b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22)
        a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22)
        c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20)
        return lambda x,y,z: x*x/(a*a) + y*y/(b*b) + z*z/(c*c) < 1

    def ell_x_shift(am, k22, k20, x_shift):
        b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22)
        a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22)
        c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20)
        return lambda x,y,z: (x - x_shift)**2/(a*a) + y*y/(b*b) + z*z/(c*c) < 1

    def tet(am, tet_shrink=1):
        tet_corners = 1.82688329031 * am * np.array([
            (1, 0, -1/np.sqrt(2)*tet_shrink),
            (-1, 0, -1/np.sqrt(2)*tet_shrink),
            (0, 1, 1/np.sqrt(2)*tet_shrink),
            (0, -1, 1/np.sqrt(2)*tet_shrink)])
        tet_norms = []
        for i in range(len(tet_corners)):
            v1 = tet_corners[i]
            v2 = tet_corners[(i+1)%4]
            v3 = tet_corners[(i+2)%4]
            normal = np.cross(v2-v1, v3-v1)
            tet_norms.append(normal / norm(normal))
        d = [np.dot(tet_norms[i], tet_corners[i]) for i in range(4)]
        return lambda x,y,z: np.all([np.sum(np.array([x, y, z]).transpose() * tet_norms[i], axis=-1) / d[i] < 1 for i in range(4)], axis=0)

    def dumbbell(am):
        db_rad = am / np.sqrt(19/20)
        return lambda x,y,z: np.minimum(np.sum(np.array([x - db_rad/2,y,z])**2, axis=0), np.sum(np.array([x + db_rad/2,y,z])**2, axis=0)) < db_rad*db_rad
        