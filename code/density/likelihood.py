import numpy as np
from core import Method
from scipy.linalg import pinv

class Likelihood(Method):
    def __init__(self, asteroid):
        print("Likelihood model")
        super().__init__(asteroid, True)

        self.grid_indices = np.cumsum(self.asteroid.indicator_map.reshape(-1)) # Indices in 1d row
        self.grid_indices = self.grid_indices.reshape(self.asteroid.indicator_map.shape) - 1# 3d row
        self.m = np.nanmax(self.grid_indices)


    def get_a(self):
        rlms = self.asteroid.moment_field()
        mat = []
        i = 0
        for l in range(0, self.asteroid.max_l+1):
            for m in range(-l, l+1):
                mat.append(self.asteroid.columnate(rlms[i]) * self.asteroid.division**3 / self.asteroid.am**l)
                i += 1
        mat.append(self.asteroid.columnate(rlms[i]) * self.asteroid.division**3 / self.asteroid.am**2)
        return pinv(np.array(mat))


    def get_b(self, x,y,z):
        ix = (x - self.asteroid.grid_line[0]) // self.asteroid.division
        iy = (y - self.asteroid.grid_line[0]) // self.asteroid.division
        iz = (z - self.asteroid.grid_line[0]) // self.asteroid.division
        return self.grid_indices[iy, ix, iz]
