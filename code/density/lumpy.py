import numpy as np
from core import Method, rlm
from scipy.linalg import pinv
from scipy.optimize import minimize

class Lumpy(Method):
    def __init__(self, asteroid):
        print("Lumpy model")
        super().__init__(asteroid, False)
        self.N = 1
        self.d = 3
        if self.d == 1:
            self.lms = []
            lms_dof = 1
        elif self.d == 2:
            self.lms = []
            lms_dof = 4
            raise NotImplementedError
        elif self.d == 3:
            self.lms = [(2, -2), (2, -1), (2, 0), (2, 1), (2, 2)]
            lms_dof = 6
        else:
            raise Exception(f"For the lumpy model, d must be 1, 2, or 3. Not {self.d}")
        self.dof = (lms_dof + 1 + 3) * self.N - 3

        self.shell_lms = self.calc_shell()
        self.x0 = np.array([1, 0,0])# Get this from shell lms

        self.xs_result, self.masses_result = self.get_positions()


    def calc_shell(self):
        rlms = self.asteroid.moment_field()
        klms = np.zeros((self.asteroid.max_l + 1)**2 + 1, dtype=np.complex)
        i = 0
        for l in range(self.asteroid.max_l + 1):
            for m in range(-l, l+1):
                klms[i] = np.sum(rlms[i]) / self.asteroid.am**l * self.asteroid.division**3
                i += 1
        klms[i] = np.sum(rlms[i]) / self.asteroid.am**2 * self.asteroid.division**3
        return klms


    def get_positions(self):
        def minim(args):
            masses = np.append(args[-1], args[::4][:-1])
            poses = np.array([args[1::4], args[2::4], args[3::4]]).transpose().reshape(-1,3) * self.asteroid.am
            m0 = 1 - np.sum(masses)
            x1 = (m0 * self.x0 - np.dot(masses[1:], poses)) / masses[0]
            masses = np.append(m0, masses)
            poses = np.vstack((self.x0, x1, poses))
            m = self.get_m(poses)
            c = self._get_c(poses, masses)
            print(m, c)
            vec = (m @ pinv(m) - np.diag(np.ones_like(self.asteroid.data, dtype=float))) @ (self.asteroid.data - c)
            return np.sum(vec * vec.conj()).real

        print(minim([1]))
        print(minim([0.1]))

        #bounds = [(-1,1), (-2, 2), (-2, 2), (-2, 2)] * self.N
        #start_x = [1, 0, 0, 0] * self.N
        #result = minimize(minim, start_x[:-3], method="L-BFGS-B", bounds=bounds[:-3])
        #print(result)#

        masses = result.x[::4]
        poses = np.array([result.x[1::4], result.x[2::4], result.x[3::4]]).transpose() * self.asteroid.am
        return poses, masses
        


    def get_a(self):
        return pinv(self.get_m(self.xs_result))


    def get_b(self, x,y,z):
        raise NotImplementedError

    
    def get_c(self):
        return self._get_c(self.xs_result, self.masses_result)


    def _get_c(self, xs, masses):
        c = np.zeros((self.asteroid.max_l + 1)**2 + 1, dtype=np.complex)
        c_index = 0
        for l in range(0, self.asteroid.max_l + 1):
            for m in range(-l, l + 1):
                # Shift shell
                c[c_index] += self.shell_lms[c_index] # Shells are already shifted

                # Shift masses
                for i in range(1, self.N):
                    c[c_index] += rlm(l, m, xs[i,0], xs[i,1], xs[i,2]) * masses[i] * (-1)**(l) / self.asteroid.am**2
                c_index += 1

        # Shift shell
        c[c_index] += self.shell_lms[-1] # Already normalized

        # Shift masses
        for i in range(1, self.N): # Do not include shell
            c[c_index] += xs[i,0]**2 + xs[i,1]**2 + xs[i,2]**2 * masses[i] / self.asteroid.am**2

        return c


    def get_m(self, xs):
        mat = np.zeros(((self.asteroid.max_l + 1)**2 + 1, self.dof), dtype=np.complex)
        lm_index = 0
        for l in range(0, self.asteroid.max_l + 1):
            for m in range(-l, l + 1):
                d_index = 0
                for i in range(1, self.N):
                    for lp, mp in self.lms:
                        if lp > l:
                            mat[lm_index, d_index] = 0
                        else:
                            print(rlm(l - lp, m - mp, xs[i,0], xs[i,1], xs[i,2]) * (-1)**(l-lp) / self.asteroid.am**(l - lp))
                            mat[lm_index, d_index] = rlm(l - lp, m - mp, xs[i,0], xs[i,1], xs[i,2]) * (-1)**(l-lp) / self.asteroid.am**(l - lp)
                        d_index += 1
                    mat[lm_index, d_index] = 0
                    d_index += 1
                lm_index += 1

        # a_m portion
        d_index = 0
        for i in range(1, self.N):
            for lp, mp in self.lms:
                mat[lm_index, d_index] = 0
                d_index += 1
            mat[lm_index, d_index] = 1
            d_index += 1

        return mat
        