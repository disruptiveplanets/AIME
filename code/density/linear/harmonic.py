import numpy as np
from core import Method, rlm
from scipy.linalg import pinv
from scipy.special import factorial

class Harmonic(Method):
    def __init__(self, asteroid):
        print("Harmonic model")
        super().__init__(asteroid, False)
        if self.asteroid.max_l is not None:
            self.m = (self.asteroid.max_l + 1)**2
    
    def get_method_name(self):
        return "harmonic"

    def get_a(self):
        rlms = self.asteroid.moment_field()
        mat = np.zeros((len(rlms), len(rlms)-1), dtype=np.complex)
        for i in range(len(rlms)): # Primed
            if i != len(rlms) - 1:
                lp = int(np.sqrt(i))
            else:
                lp = 2
            for j in range(len(rlms) - 1): # Unprimed
                l = int(np.sqrt(j))
                mat[i,j] = np.sum(rlms[i] * rlms[j].conj()) / self.asteroid.am**(lp + l) * self.asteroid.division**3

        return pinv(mat)


    def get_b(self, x, y, z):
        b = np.zeros(self.m, dtype=np.complex)
        i = 0
        for l in range(0, self.asteroid.max_l + 1):
            for m in range(-l, l+1):
                b[i] = rlm(l, m, x, y, z).conj() / self.asteroid.am**l
                i += 1
        return b
