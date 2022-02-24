import numpy as np
from core import Method, rlm
from scipy.linalg import pinv


class Harmonic(Method):
    def __init__(self, asteroid):
        print("Harmonic model")
        super().__init__(asteroid)
        self.m = (self.asteroid.max_l + 1)**2

    def get_a(self):
        rlms = self.asteroid.moment_field()
        mat = np.zeros((len(rlms), len(rlms)-1), dtype=np.complex)
        for i in range(len(rlms)): # Primed
            for j in range(len(rlms) - 1): # Unprimed
                if i != len(rlms) - 1:
                    l = int(np.sqrt(i))
                else:
                    l = 2
                mat[i,j] = np.sum(rlms[i] * rlms[j].conj()) / self.asteroid.am**l * self.asteroid.division**3
        return pinv(mat)


    def get_b(self, x, y, z):
        b = np.zeros(self.m, dtype=np.complex)
        i = 0
        for l in range(0, self.asteroid.max_l + 1):
            for m in range(-l, l+1):
                b[i] = rlm(l, m, x, y, z).conj()
                i += 1
        return b
