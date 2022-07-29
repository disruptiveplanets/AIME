import numpy as np
from core import Method
from scipy.linalg import inv

NUM_DRAWS = 10
SHIFT_FRACTION = 1

class FiniteElement(Method):
    def __init__(self, asteroid):
        print("Finite element linear model")
        super().__init__(asteroid, False)

        if self.asteroid.max_l is not None:
            self.num_grids = (self.asteroid.max_l + 1)**2 + 1
        self.points, self.grids = self.get_grids_centroid(self.num_grids)

    def get_grids_centroid(self, num_grids):
        # Shift each point to its cell's centroid and iterate
        points = self.get_rand_points(num_grids)
        x,y,z = np.meshgrid(self.asteroid.grid_line, self.asteroid.grid_line, self.asteroid.grid_line)
        grids = None
        for i in range(NUM_DRAWS):
            grids = self.get_grids_voronoi(points)
            for grid_index, grid in enumerate(grids):
                centroid = np.array([np.sum(grid * x), np.sum(grid * y), np.sum(grid * z)]) / np.sum(grid)
                points[grid_index] = points[grid_index] + (centroid - points[grid_index]) * SHIFT_FRACTION
            print(np.sort(np.sum(grids, axis=(1,2,3))))
            #print(np.sum(grids, axis=(1,2,3)))
        grids = self.get_grids_voronoi(points)
        return points, grids

    def get_grids_weights(self, num_grids):
        # Weight the points with small volumes stronger and iterate
        points = self.get_rand_points(num_grids)
        weights = np.ones(num_grids)
        grids = None
        equal_vol = np.sum(self.asteroid.indicator_map) / num_grids
        for i in range(NUM_DRAWS):
            grids = self.get_grids_weight(points, weights)
            vols = np.sum(grids, axis=(1,2,3)) / equal_vol
            weights = (1 / vols)**0.3
            print(vols)
        grids = self.get_grids_weight(points, weights)
        return points, grids

    def get_rand_points(self, num_grids):
        points = []
        while len(points) < num_grids:
            p = np.random.uniform(low=np.min(self.asteroid.grid_line), high=np.max(self.asteroid.grid_line), size=3)
            if self.asteroid.indicator(p[0], p[1], p[2]):
                points.append(p)
        return np.array(points)

    def get_grids_voronoi(self, points):
        x,y,z = np.meshgrid(self.asteroid.grid_line, self.asteroid.grid_line, self.asteroid.grid_line)
        dists = np.zeros((len(points), len(self.asteroid.grid_line), len(self.asteroid.grid_line), len(self.asteroid.grid_line)))
        for i in range(len(points)):
            dists[i] = (x - points[i,0]) ** 2 + (y - points[i,1]) ** 2 + (z - points[i,2]) ** 2
        grid_indices = np.argmin(dists, axis=0)

        grids = []
        for i in range(len(points)):
            grids.append((grid_indices == i) * self.asteroid.indicator_map)

        return np.array(grids)
    
    def get_grids_weight(self, points, weights):
        x,y,z = np.meshgrid(self.asteroid.grid_line, self.asteroid.grid_line, self.asteroid.grid_line)
        dists = np.zeros((len(points), len(self.asteroid.grid_line), len(self.asteroid.grid_line), len(self.asteroid.grid_line)))
        for i in range(len(points)):
            dists[i] = (x - points[i,0]) ** 2 + (y - points[i,1]) ** 2 + (z - points[i,2]) ** 2 / weights[i]**2
        grid_indices = np.argmin(dists, axis=0)

        grids = []
        for i in range(len(points)):
            grids.append((grid_indices == i) * self.asteroid.indicator_map)

        return np.array(grids)

    def get_method_name(self):
        return "fe"

    def get_a(self):
        rlms = self.asteroid.moment_field()
        mat = np.zeros((self.num_grids, self.num_grids), dtype=np.complex)
        klm_index = 0
        for l in range(self.asteroid.max_l + 1):
            for _ in range(-l, l + 1):
                for grid_index in range(self.num_grids):
                    mat[klm_index, grid_index] = np.sum(self.grids[grid_index] * rlms[klm_index] / self.asteroid.am**l)
                klm_index += 1
        for grid_index in range(self.num_grids):
            mat[klm_index, grid_index] = np.sum(self.grids[grid_index] * rlms[klm_index] / self.asteroid.am**2) # Radius term

        mat *= self.asteroid.division**3
        
        return inv(mat)

    def get_b(self, x,y,z):
        dist = [(x - p[0])**2 + (y - p[1])**2 + (z - p[2])**2 for p in self.points]
        indices = np.argmin(dist, axis=0)
        b = np.zeros([self.num_grids] + list(x.shape))
        b[indices] = 1
        return b
