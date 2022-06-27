from re import L
import numpy as np
from scipy.special import lpmv, factorial
from scipy.linalg import norm, eigh, pinv
from scipy.spatial.transform import Rotation
from display import make_gif, make_slices
import warnings, os

RLM_EPSILON = 1e-20
MAX_L = 3
NEGATIVE_VEC = None

def rlm(l,m,x,y,z):
    r = np.sqrt(np.maximum(RLM_EPSILON, x*x + y*y + z*z))
    return lpmv(m, l, z/r) / factorial(l + m) * r**l * np.exp(1j * m * np.arctan2(y, x))

def rlm_gen(l,m):
    return lambda x,y,z: rlm(l,m,x,y,z)

class Asteroid:
    def __init__(self, name, surface_am, division, max_radius, indicator, true_shape, true_moments=None):
        self.name = name
        self.surface_am = surface_am
        self.true_shape = true_shape
        self.true_moments = true_moments
        self.true_densities = None
        self.division = division
        self.grid_line = np.arange(-max_radius, max_radius, division)
        self.indicator = indicator
        self.indicator_map = self.get_indicator_map()
        self.moments = None

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


    def moment_field(self, surface_am):
        if self.moments is None:
            n = (MAX_L+1)**2 + 1
            self.moments = np.zeros((n, len(self.grid_line), len(self.grid_line), len(self.grid_line)), dtype=np.complex)
            i = 0
            for l in range(0, MAX_L+1):
                renorm = surface_am**(2 - l) if l > 0 else 1
                for m in range(-l, l+1):
                    self.moments[i] = self.map_np(rlm_gen(l,m)) * renorm
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

    def ell_3_shift(am, k22, k20, shift):
        a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22)
        b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22)
        c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20)
        return lambda x,y,z: (x - shift[0])**2/(a*a) + (y - shift[1])**2/(b*b) + (z - shift[2])**2 /(c*c) < 1

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

    def core(am, k22, k20, density_ratio, radius_ratio):
        a = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 + 12 * k22) * radius_ratio
        b = np.sqrt(5/3) * am * np.sqrt(1 - 2 * k20 - 12 * k22) * radius_ratio
        c = np.sqrt(5/3) * am * np.sqrt(1 + 4 * k20) * radius_ratio
        def func(x, y, z):
            dist = x*x / (a*a) + y*y / (b*b) + z*z / (c*c)
            return 1 * (dist > 1) + density_ratio * (dist <= 1)
        return func

    def two_core(density_ratio_1, radius_1, pos_1, density_ratio_2, radius_2, pos_2):
        def func(x, y, z):
            dist1 = (x - pos_1[0])**2 + (y-pos_1[1])**2 + (z - pos_1[2])**2
            dist2 = (x - pos_2[0])**2 + (y-pos_2[1])**2 + (z - pos_2[2])**2
            return 1 + (density_ratio_1 - 1) * (dist1 <= radius_1 * radius_1) + (density_ratio_2 - 1) * (dist2 <= radius_2 * radius_2)
        return func

    def core_sph(density_ratio, radius):
        def func(x, y, z):
            dist = x*x + y*y + z*z
            return 1 * (dist > radius * radius) + density_ratio * (dist <= radius * radius)
        return func

    def core_shift(density_ratio, radius, shift):
        def func(x, y, z):
            dist = x*x + (y-shift)*(y-shift) + z*z
            return 1 * (dist > radius * radius) + density_ratio * (dist <= radius * radius)
        return func

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

class TrueMoments:
    # Complex klms for the surface, bulk_klms
    def ell():
        k22, k20 = -0.05200629, -0.2021978
        return np.array([
            1.0, 0, 0, 0,
            k22, 0, k20, 0, k22,
            0, 0, 0, 0, 0, 0, 0
        ])
    def move_3():
        k22, k20 = -0.052298941182752544, -0.2023903027848194
        ik31, ik33 = -0.00801768900716791j, -0.0016420155228371945j
        return np.array([
            1.0,
            0.031607016470795626j, 0, 0.031607016470795626j,
            k22, 0, k20, 0, k22,
            ik33 * 1j, 0,ik31 * 1j, 0, ik31 * 1j, 0, ik33 * 1j
        ])
    def move_1_5():
        k22, k20 = -0.052024007380816216, -0.2022120206327561
        ik31, ik33 = -0.0020127008643331624j, -0.00041162951899840697j
        return np.array([
            1.0,
            0.007919641491764374j, 0, 0.007919641491764374j,
            k22, 0, k20, 0, k22,
            ik33 * 1j, 0,ik31 * 1j, 0, ik31 * 1j, 0, ik33 * 1j
        ])

class UncertaintyTracker:
    def __init__(self):
        self.num_maps = 0
        self.density_map_sum = None
        self.density_map_square_sum = None
    
    def update(self, density_map):
        if self.density_map_sum is None:
            self.density_map_sum = density_map
            self.density_map_square_sum = density_map**2
        else:
            self.density_map_sum += density_map
            self.density_map_square_sum += density_map**2
        self.num_maps += 1

    def generate(self):
        if self.num_maps == 0:
            return None, None
        mean_map = self.density_map_sum / self.num_maps
        return (mean_map, 
            np.sqrt(self.density_map_square_sum / self.num_maps - mean_map**2)
        )

    def __iadd__(self, other_tracker):
        if other_tracker.density_map_sum is None:
            return self
        if self.density_map_sum is None:
            self.num_maps = other_tracker.num_maps
            self.density_map_sum = other_tracker.density_map_sum
            self.density_map_square_sum = other_tracker.density_map_square_sum
        else:
            self.num_maps += other_tracker.num_maps
            self.density_map_sum += other_tracker.density_map_sum
            self.density_map_square_sum += other_tracker.density_map_square_sum
        return self

    def save(self, f):
        np.save(f, (
            self.density_map_sum, 
            self.density_map_square_sum,
            np.ones_like(self.density_map_sum) * self.num_maps)
        )

    def load(f):
        data = np.load(f)
        if len(data) == 3:
            map_sum, square_sum, nums = data
        else:
            print(f"Improperly formatted file: {data}")
            return UncertaintyTracker()
        unc_tracker = UncertaintyTracker()
        unc_tracker.density_map_sum = map_sum
        unc_tracker.density_map_square_sum = square_sum
        unc_tracker.num_maps = nums[0, 0, 0]
        return unc_tracker