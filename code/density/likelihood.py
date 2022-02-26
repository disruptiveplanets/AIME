import numpy as np
from core import Method
from scipy.linalg import pinv
import scipy.sparse

class Likelihood(Method):
    def __init__(self, asteroid):
        print("Likelihood model")
        super().__init__(asteroid)
        self.grid_indices = np.zeros((len(self.asteroid.grid_line), len(self.asteroid.grid_line), len(self.asteroid.grid_line)), dtype=int)
        self.m = 0
        for nx, x in enumerate(self.asteroid.grid_line):
            for ny, y in enumerate(self.asteroid.grid_line):
                for nz, z in enumerate(self.asteroid.grid_line):
                    if self.asteroid.indicator_map[nx,ny,nz]:
                        self.grid_indices[nx,ny,nz] = self.m
                        self.m += 1


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
        if len(x.shape) > 0:
            out = np.zeros((x.shape[0], x.shape[0], x.shape[0], self.m))
            out[ix, iy, iz, self.grid_indices[ix,iy,iz]] = 1
        else:
            out = np.zeros(self.m)
            out[self.grid_indices[ix, iy, iz]] = 1
        return scipy.sparse.coo_matrix(out)
