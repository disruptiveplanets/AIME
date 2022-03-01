import numpy as np
from core import Method, rlm
from scipy.linalg import pinv
from scipy.optimize import minimize, fsolve
from scipy.special import factorial
from scipy.spatial.transform import Rotation

class Lumpy(Method):
    def __init__(self, asteroid):
        print("Lumpy model")
        super().__init__(asteroid, False, final_uncertainty=False)
        self.N = 1
        self.rot_dim = 3
        if self.rot_dim == 1:
            self.lms = []
            self.lms_dof = 1
        elif self.rot_dim == 2:
            self.lms = []
            self.lms_dof = 4
            raise NotImplementedError
        elif self.rot_dim == 3:
            self.lms = [(2, -2), (2, -1), (2, 0), (2, 1), (2, 2)]
            self.lms_dof = 6
        else:
            raise Exception(f"For the lumpy model, d must be 1, 2, or 3. Not {self.rot_dim}")
        self.dof = (self.lms_dof + 1 + 3) * self.N - 3

        self.shell_lms = self.calc_shell()
        self.x0 = np.array([(self.shell_lms[1]-self.shell_lms[3]).real, (1j*(self.shell_lms[1]+self.shell_lms[3])).real, self.shell_lms[2].real])
        self.intersection = None
        self.xs_result, self.masses_result = self.get_positions()
        print("X result", self.xs_result)
        print("Mass result", self.masses_result)


    def calc_shell(self):
        rlms = self.asteroid.moment_field()
        klms = np.zeros((self.asteroid.max_l + 1)**2 + 1, dtype=np.complex)
        i = 0
        for l in range(self.asteroid.max_l + 1):
            for m in range(-l, l+1):
                klms[i] = np.sum(rlms[i]) / self.asteroid.am**l * self.asteroid.division**3
                i += 1
        klms[i] = np.sum(rlms[i]) / self.asteroid.am**2 * self.asteroid.division**3
        self.volumes = [klms[0]]
        return klms / klms[0] # Normalize to mass = 1 for now


    def get_positions(self):
        def minim(args):
            masses = np.append(args[-1], args[::4][:-1])
            poses = np.array([args[1::4], args[2::4], args[3::4]]).transpose().reshape(-1,3) * self.asteroid.am
            m0 = 1 - np.sum(masses)
            x1 = -(m0 * self.x0 + np.dot(masses[1:], poses)) / masses[0]
            if np.any(np.abs(x1) > 2):
                # Enforce boundary conditions on x1
                return np.inf
            masses = np.append(m0, masses)
            poses = np.vstack((self.x0, x1, poses))
            m = self.get_m(poses)
            c = self._get_c(poses, masses)
            np.set_printoptions(linewidth=150)
            vec = (m @ pinv(m) - np.diag(np.ones_like(self.asteroid.data, dtype=float))) @ (self.asteroid.data - c)
            return np.sum(vec * vec.conj()).real

        bounds = [(-1,1), (-2, 2), (-2, 2), (-2, 2)] * self.N
        start_x = [0.1, 0, 0, 0] * self.N
        result = minimize(minim, start_x[:-3], method="Nelder-Mead", bounds=bounds [:-3])
        print(result)

        print(minim([0.3]))

        masses = np.append(result.x[-1], result.x[::4][:-1])
        poses = np.array([result.x[1::4], result.x[2::4], result.x[3::4]]).transpose().reshape(-1,3) * self.asteroid.am
        m0 = 1 - np.sum(masses)
        x1 = -(m0 * self.x0 + np.dot(masses[1:], poses)) / masses[0]
        masses = np.append(m0, masses)
        poses = np.vstack((self.x0, x1, poses))
        return poses, masses


    def unrotate(self, k2ms):
        if self.rot_dim == 1:
            return np.diag(np.ones(3)), 0, 0
        def mod_wigner_d(m, mp, alpha, beta, gamma):
            smin = max(0, mp-m)
            smax = min(2 + mp, 2 - m)
            d = 0
            for s in range(smin, smax + 1):
                d += (-1)**(m-mp+s) * np.cos(beta/2) ** (4 + mp - m - 2*s) * np.sin(beta/2)**(m-mp+2*s) \
                    / (factorial(2+mp-s)*factorial(2-m-s)*factorial(m-mp+s)*factorial(s))
            return d * np.exp(-1j * m * alpha) * np.exp(-1j * mp * gamma) * factorial(2+mp)*factorial(2-mp)

        def object_fn(args):
            dmat = np.array([[mod_wigner_d(m, mp, args[0], args[1], args[2]).conj() for m in range(-2, 3)] for mp in range(-2, 3)]).transpose()
            rotated = dmat @ k2ms



######### THE ABOVE VERY WELL MAY BE BACKWARDS> MAYBE DMAT SHOULD NOT BE TRANSPOSED



            return rotated[1].real, rotated[1].imag, rotated[0].imag

        result = fsolve(object_fn, x0=(1, 1, 1))
        rot_mat = Rotation.from_euler('zyz', [result[2], result[1], result[0]]).as_matrix()


######### THE ABOVE COULD ALSO BE WRONG
        dmat = np.array([[mod_wigner_d(m, mp, result[0], result[1], result[2]).conj() for m in range(-2, 3)] for mp in range(-2, 3)]).transpose()
        rotated = dmat @ k2ms
######## REMEMBER TO CHANGE THIS IF YOU CHANGE THE ABOVE

        return rot_mat, rotated[0].real, rotated[2].real


    def get_intersections(self, k2ms):
        if self.intersection is None:
            print(self.d)

            self.intersection = []
            for i in range(1, self.N + 1):
                rot_mat, k22, k20 = self.unrotate(k2ms[:-1])
                a = np.sqrt(5/3 * k2ms[-1]) * np.sqrt(1 - 2 * k20 + 12 * k22)
                b = np.sqrt(5/3 * k2ms[-1]) * np.sqrt(1 - 2 * k20 - 12 * k22)
                c = np.sqrt(5/3 * k2ms[-1]) * np.sqrt(1 + 4 * k20)
                self.volumes.append(a * b * c * np.pi * 4 / 3)
                def check_intersection(x, y, z):
                    rot_vec = rot_mat @ np.array([x, y, z])
                    return rot_vec[0]**2 / (a*a) + rot_vec[1]**2 / (b*b) + rot_vec[2]**2 / (c*c) < 1
                self.intersection.append(check_intersection)
        return self.intersection
        

    def get_a(self):
        #print("M:", self.get_m(self.xs_result))
        return pinv(self.get_m(self.xs_result))


    def get_loc_density(self, k2ms, x, y, z):
        intersections = self.get_intersections(k2ms)
        density = 0
        for i in range(1, self.N + 1):
            if intersections[i-1](x, y, z):
                density += self.masses_result[i] / self.volumes[i]
        density += self.masses_result[0] / self.volumes[0]
        return density
        

    def get_b(self, x,y,z):
        raise NotImplementedError

    
    def get_c(self):
        print("C:", self._get_c(self.xs_result, self.masses_result))
        return self._get_c(self.xs_result, self.masses_result)


    def _get_c(self, xs, masses):
        c = np.append([1, 0, 0, 0], np.zeros((self.asteroid.max_l + 1)**2 - 3, dtype=np.complex))
        c_index = 4
        # Skip l = 0 and 1, which is boring
        for l in range(2, self.asteroid.max_l + 1):
            for m in range(-l, l + 1):
                # Shift shell
                c[c_index] += masses[0] * self.shell_lms[c_index] # Shells are already shifted & normalized

                # Shift masses
                for i in range(1, self.N + 1):
                    c[c_index] += rlm(l, m, xs[i,0], xs[i,1], xs[i,2]) * masses[i] * (-1)**(l) / self.asteroid.am**l
                c_index += 1

        # Shift shell
        c[c_index] += self.shell_lms[-1] # Already normalized

        # Shift masses
        for i in range(1, self.N + 1): # Do not include shell
            c[c_index] += (xs[i,0]**2 + xs[i,1]**2 + xs[i,2]**2) * masses[i] / self.asteroid.am**2

        return c


    def get_m(self, xs):
        mat = np.zeros(((self.asteroid.max_l + 1)**2 + 1, self.lms_dof), dtype=np.complex)
        lm_index = 0
        for l in range(0, self.asteroid.max_l + 1):
            for m in range(-l, l + 1):
                d_index = 0
                for i in range(1, self.N + 1):
                    for lp, mp in self.lms:
                        if lp > l or abs(m-mp) > l - lp:
                            mat[lm_index, d_index] = 0
                        else:
                            mat[lm_index, d_index] = rlm(l - lp, m - mp, xs[i,0], xs[i,1], xs[i,2]) * (-1)**(l-lp) / self.asteroid.am**(l - lp)
                        d_index += 1
                    mat[lm_index, d_index] = 0
                    d_index += 1
                lm_index += 1

        # a_m portion
        d_index = 0
        for i in range(1, self.N + 1):
            for lp, mp in self.lms:
                mat[lm_index, d_index] = 0
                d_index += 1
            mat[lm_index, d_index] = 1
            d_index += 1

        return mat
        