from re import L
import numpy as np
from scipy.special import lpmv, factorial
from scipy.linalg import norm, eigh, pinv
from scipy.spatial.transform import Rotation
from display import make_gif, make_slices
import warnings, os

RLM_EPSILON = 1e-20
NEGATIVE_VEC = None


def rlm(l,m,x,y,z):
    r = np.sqrt(np.maximum(RLM_EPSILON, x*x + y*y + z*z))
    return lpmv(m, l, z/r) / factorial(l + m) * r**l * np.exp(1j * m * np.arctan2(y, x))

def rlm_gen(l,m):
    return lambda x,y,z: rlm(l,m,x,y,z)

class Method:
    def __init__(self, asteroid, finite_element, final_uncertainty=True):
        self.asteroid = asteroid
        self.unc = None
        self.d = None
        self.densities = None
        self.density_uncs = None
        self.finite_element = finite_element
        self.final_uncertainty = final_uncertainty
        self.klm_error = None

    def get_method_name(self):
        raise NotImplementedError

    def get_bare_name(self):
        return f"{self.asteroid.name}/{self.get_method_name()}"

    def get_a(self):
        raise NotImplementedError
        
    def get_b(self, pos):
        # Return number for finite_element, vector otherwise
        raise NotImplementedError

    def get_c(self):
        # Default value.
        return np.zeros_like(self.asteroid.data)

    def get_loc_density(self, k2ms, x, y, z):
        raise NotImplementedError

    def solve(self):
        if self.unc is None:
            a = self.get_a()
            self.d = a @ (self.asteroid.data - self.get_c())

            if self.asteroid.sigma_data is not None:
                if self.finite_element:
                    self.unc = np.einsum('ij,jk,ki->i', a, self.asteroid.sigma_data, a.transpose().conjugate())# Diagonal entries
                else:
                    self.unc = a @ self.asteroid.sigma_data @ a.transpose().conjugate()

    def get_density(self, x,y,z):
        if self.finite_element:
            out = self.d[self.get_b(x,y,z)]
        elif not self.final_uncertainty:
            out = self.get_loc_density(self.d, x, y, z)
        else:
            out = np.dot(self.get_b(x,y,z), self.d)

        if type(out) != np.float64 and type(out) != np.complex128:
            return out[0].real
        return out.real

    def get_unc(self, x,y,z):
        if self.final_uncertainty:
            b = self.get_b(x,y,z)
            if self.finite_element:
                out = np.sqrt(np.abs(self.unc[b].real)).real
                return out

            else:
                if len(x.shape) > 0:
                    b = b.reshape(-1, b.shape[-1])
                    return np.array([np.sqrt(e @ self.unc @ e)[0] for e in b]).reshape(x.shape).real
                else:
                    out = np.sqrt(np.abs((b.transpose() @ self.unc @ b.conjugate()).real))
                    if type(out) != np.float64:
                        return out[0].real
                    return out.real

    def map_density(self):
        if self.densities is None:
            self.densities = self.asteroid.map(self.get_density, dtype=float)
        return self.densities

    def map_unc(self):
        if self.density_uncs is None:
            self.density_uncs = self.asteroid.map(self.get_unc, dtype=float) / np.abs(self.map_density())
        return self.density_uncs

    def save_density(self):
        fname = f"data/{self.get_bare_name()}-d.npy"
        if not os.path.exists(f"data/{self.asteroid.name}"):
            print("Make")
            os.mkdir(f"data/{self.asteroid.name}")
        with open(fname, 'wb') as f:
            np.save(f, self.map_density())

    def save_unc(self):
        fname = f"data/{self.get_bare_name()}-u.npy"
        if not os.path.exists(f"data/{self.asteroid.name}"):
            os.mkdir(f"data/{self.asteroid.name}")
        with open(fname, 'wb') as f:
            np.save(f, self.map_unc())

    def check(self, display=True):
        rlms_field = self.asteroid.moment_field()
        calc_rlms = np.zeros((self.asteroid.max_l + 1)**2 + 1, dtype=np.complex)
        index = 0
        densities = self.map_density()
        for l in range(0, self.asteroid.max_l + 1):
            for m in range(-l, l+1):
                calc_rlms[index] = np.sum(rlms_field[index] * densities) / self.asteroid.am ** l * self.asteroid.division**3
                index += 1
        calc_rlms[index] = np.sum(rlms_field[index] * densities) / self.asteroid.am ** 2 * self.asteroid.division**3

        # Calculate chisqr
        uncertain_params = []
        for l in range(self.asteroid.max_l + 1):
            for m in range(-l, l+1):
                if l == 0 or l == 1:
                    uncertain_params.append(False)
                elif l == 2 and (abs(m) == 1 or m == -2):
                    uncertain_params.append(False)
                else:
                    uncertain_params.append(True)
        uncertain_params = np.append(uncertain_params, False)
        dof = (self.asteroid.max_l + 1)**2 - 7

        deltas = (calc_rlms - self.asteroid.data)[uncertain_params]
        nonzero_cov = self.asteroid.sigma_data[uncertain_params,:][:,uncertain_params]
        self.klm_error = abs(deltas.transpose().conj() @ pinv(nonzero_cov) @ deltas / dof)

        for i, r in enumerate(calc_rlms):
            l = int(np.sqrt(i))
            m = i - l**2 - l
            if display:
                print("({}, {})\tExpected  {:.5f}\t Got {:.5f} \t Difference {:.2g}".format(l, m, self.asteroid.data[i], r, abs(self.asteroid.data[i]-r)))
        if display:
            print("Reduced chi squared:", self.klm_error)

    def display(self, duration=5):
        fname = f"figs/{self.get_bare_name()}"
        asteroid_name = fname.split("/")[1]
        if not os.path.isdir(f"figs/{asteroid_name}"):
            os.mkdir(f"figs/{asteroid_name}")
        if self.klm_error is None:
            self.check(False)

        warnings.filterwarnings("ignore")

        display_densities = np.copy(self.map_density())
        display_densities[~self.asteroid.indicator_map] = np.nan
        if self.final_uncertainty:
            display_uncs = np.copy(self.map_unc())
            display_uncs[~self.asteroid.indicator_map] = np.nan

        display_densities /= np.nanmean(display_densities)

        true_densities = self.asteroid.get_true_densities()
        
        if true_densities is not None:
            true_densities[~self.asteroid.indicator_map] = np.nan
            true_densities /= np.nanmean(true_densities)
            if self.final_uncertainty:
                ratios = (true_densities - display_densities) / (display_densities * display_uncs)
                make_slices(ratios, self.asteroid.grid_line, "$\\Delta\\sigma$", 'coolwarm', f"{fname}-r", self.klm_error, percentile=95, balance=True)
                make_gif(ratios, self.asteroid.grid_line, "$\\Delta\\sigma$", 'coolwarm', f"{fname}-r.gif", duration=duration, percentile=95, balance=True)

                for p in 10 * np.arange(11):
                    print(f"Ratio percentile {p} = {np.nanpercentile(ratios, p)}")
            difference = (true_densities - display_densities) / true_densities

        print("Plotting density")
        make_slices(display_densities, self.asteroid.grid_line, "$\\rho$", 'plasma', f"{fname}-d", self.klm_error)
        make_gif(display_densities, self.asteroid.grid_line, "$\\rho$", 'plasma', f"{fname}-d.gif", duration)
        
        if self.final_uncertainty:
            print("Plotting uncertainty")
            make_slices(display_uncs, self.asteroid.grid_line, "$\\sigma_\\rho / \\rho$", 'Greys_r', f"{fname}-u", self.klm_error, 95)
            make_gif(display_uncs, self.asteroid.grid_line, "$\\sigma_\\rho / \\rho$", 'Greys_r', f"{fname}-u.gif", duration, 95)
        if true_densities is not None:
            print("Plotting differences")
            make_slices(difference, self.asteroid.grid_line, "$\\Delta\\rho$", 'PuOr_r', f"{fname}-s", self.klm_error, 95, balance=True)
            make_gif(difference, self.asteroid.grid_line, "$\\Delta\\rho$", 'PuOr_r', f"{fname}-s.gif", duration, 95, balance=True)

        
        warnings.filterwarnings("default")

    def reload(self, density_name, unc_name):
        with open(density_name, 'rb') as f:
            self.densities = np.load(f)
        with open(unc_name, 'rb') as f:
            self.density_uncs = np.load(f)


class Asteroid:
    def __init__(self, name, sample_file, am, division, max_radius, indicator, true_shape):
        self.name = name
        self.max_l = None
        self.am = am
        self.true_shape = true_shape
        self.true_densities = None
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
            klms[15] = samples[3] + 1j * samples[4]
            klms[14] = samples[5] + 1j * samples[6]
            klms[13] = samples[7] + 1j * samples[8]
            klms[12] = samples[9]
            klms[11] = -klms[13].conj()
            klms[10] = klms[14].conj()
            klms[9] = -klms[15].conj()
            ms = np.append(ms, [-3, -2, -1, 0, 1, 2, 3])

        delta_gamma = samples[0] - np.mean(samples[0])
        hybrids = klms * np.exp(-1j * np.outer(ms, delta_gamma))
        klm_means = np.mean(hybrids, axis=1)

        klm_cov = np.cov(hybrids)
        added_one_cov = np.hstack((klm_cov, np.zeros((klm_cov.shape[0], 1))))
        added_cov = np.vstack((added_one_cov, np.zeros((1, klm_cov.shape[0]+1))))

        return np.append(klm_means, 1), added_cov

    def get_true_densities(self):
        if self.true_densities is None:
            if self.true_shape is None:
                return None
            self.true_densities = self.map_np(self.true_shape)

        return self.true_densities


    def get_indicator_map(self):
        x,y,z = np.meshgrid(self.grid_line, self.grid_line, self.grid_line)
        return self.indicator(x,y,z)


    def map(self, fn, dtype=np.complex):
        x,y,z = np.meshgrid(self.grid_line, self.grid_line, self.grid_line)
        result = np.zeros_like(x, dtype=dtype)
        for ni in range(len(self.grid_line)):
            for nj in range(len(self.grid_line)):
                for nk in range(len(self.grid_line)):
                    if self.indicator_map[ni,nj,nk]:
                        result[ni, nj, nk] = fn(x[ni,nj,nk], y[ni,nj,nk], z[ni,nj,nk])
        return result

    def map_np(self, fn):
        x,y,z = np.meshgrid(self.grid_line, self.grid_line, self.grid_line)
        return fn(x,y,z) * self.indicator_map


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
        a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22)
        b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22)
        c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20)
        return lambda x,y,z: x*x/(a*a) + y*y/(b*b) + z*z/(c*c) < 1

    def ell_y_shift(am, k22, k20, y_shift):
        a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22)
        b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22)
        c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20)
        return lambda x,y,z: (x)**2/(a*a) + (y - y_shift)**2/(b*b) + z*z/(c*c) < 1

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
        

class TrueShape:
    def uniform():
        return lambda x, y, z: 1.0

    def in_(am):
        return lambda x, y, z: np.exp(-(x*x + y*y + z*z)/am**2)

    def out(am):
        return lambda x, y, z: np.exp((x*x + y*y + z*z)/am**2)

    def in_sph(am):
        r2 = 5/3 * am * am
        return lambda x, y, z: np.exp(-0.5 * (x*x + y*y + z*z)/r2)

    def out_sph(am):
        r2 = 5/3 * am * am
        return lambda x, y, z: np.exp(0.5 * (x*x + y*y + z*z)/r2)

    def blob(am, k22, k20):
        blob_displacement = 500
        blob_rad = 300
        density_factor = 5
        def lump(x, y, z):
            o = np.ones_like(x, dtype=float)
            o[(x)**2 + (y-blob_displacement)**2 + z**2 < blob_rad**2] += density_factor
            return o
        return lump

    def rot_blob(am, k22, k20):
        blob_displacement = 500
        blob_rad = 300
        density_factor = 5
        b = (4/3)**(1/3) * blob_rad # To keep volume of blob and rot blob the same
        a = b/2
        c = 3*b/2
        blob_rot_mat = Rotation.from_euler('zyz', (np.pi/4, np.pi/4, np.pi/4)).as_matrix()
        def lump(x, y, z):
            o = np.ones_like(x, dtype=float)
            diff = np.einsum('ij,jklm', blob_rot_mat, np.array([x,y-blob_displacement, z])) # Multiply all vectors by the rotation matrix
            o[diff[0]**2/a**2 + diff[1]**2/b**2 + diff[2]**2/c**2 < 1] += density_factor
            return o
        return lump