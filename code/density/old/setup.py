import numpy as np
from scipy.special import lpmv, factorial
from scipy.linalg import norm

DIVISION = 10
MAX_RADIUS = 2001
pos_array = np.arange(-MAX_RADIUS, MAX_RADIUS, DIVISION)
TET_SHRINK = 1#np.sqrt(2 * 0.90233392 / 1.19533216 - 1)

TYPE = 5

if TYPE == 0:
    TAG = "asym-ell"
    KLMS = [1.0, 0, 0, 0, 0.05200629, 0, 0, 0, -0.2021978]
    A_M = 1000
    a = np.sqrt(5/3) * A_M * np.sqrt(1 - 2 * KLMS[8] + 12 * KLMS[4])
    b = np.sqrt(5/3) * A_M * np.sqrt(1 - 2 * KLMS[8] - 12 * KLMS[4])
    c = np.sqrt(5/3) * A_M * np.sqrt(1 + 4 * KLMS[8])
    def indicator(pos):
        return pos[0] * pos[0]/(a*a) + pos[1] * pos[1]/(b*b) + pos[2] * pos[2]/(c*c) < 1
elif TYPE == 1:
    TAG = "sym-ell"
    KLMS = [1.0, 0, 0, 0, 0, 0, 0, 0, -0.09766608]
    A_M = 1000
    a = np.sqrt(5/3) * A_M * np.sqrt(1 - 2 * KLMS[8] + 12 * KLMS[4])
    b = np.sqrt(5/3) * A_M * np.sqrt(1 - 2 * KLMS[8] - 12 * KLMS[4])
    c = np.sqrt(5/3) * A_M * np.sqrt(1 + 4 * KLMS[8])
    def indicator(pos):
        return pos[0] * pos[0]/(a*a) + pos[1] * pos[1]/(b*b) + pos[2] * pos[2]/(c*c) < 1
elif TYPE == 2:
    TAG = "asym-sph"
    KLMS = [1.0, 0, 0, 0, 0.05200629, 0, 0, 0, -0.2021978]
    A_M = 1000
    def indicator(pos):
        return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] < (A_M * A_M * 5/3)
elif TYPE == 3:
    TAG = "sym-sph"
    KLMS = [1.0, 0, 0, 0, 0, 0, 0, 0, -0.09766608]
    A_M = 1000
    def indicator(pos):
        return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] < (A_M * A_M * 5/3)
elif TYPE == 4:
    TAG = "tet"
    KLMS = [1.0, 0, 0, 0, 0, 0, 0, 0, 0]
    A_M = 1000
    tet_corners = 1.82688329031 * A_M * np.array([
        (1, 0, -1/np.sqrt(2)*TET_SHRINK),
        (-1, 0, -1/np.sqrt(2)*TET_SHRINK),
        (0, 1, 1/np.sqrt(2)*TET_SHRINK),
        (0, -1, 1/np.sqrt(2)*TET_SHRINK)])
    tet_norms = []
    for i in range(len(tet_corners)):
        v1 = tet_corners[i]
        v2 = tet_corners[(i+1)%4]
        v3 = tet_corners[(i+2)%4]
        normal = np.cross(v2-v1, v3-v1)
        tet_norms.append(normal / norm(normal))
    d = [np.dot(tet_norms[i], tet_corners[i]) for i in range(4)]
    def indicator(pos):
        return np.all([np.dot(pos, tet_norms[i]) / d[i] < 1 for i in range(4)])
elif TYPE == 5:
    TAG = "dumb-bell"
    db_rad = 1200
    KLMS = [1.0, 0, 0, 0, 25/608, 0, 0, 0, -25/304]
    A_M = np.sqrt(19/20) * db_rad
    def indicator(pos):
        r1 = np.sum((pos - np.array([db_rad/2, 0, 0]))**2)
        r2 = np.sum((pos + np.array([db_rad/2, 0, 0]))**2)
        return r1 < db_rad*db_rad or r2 < db_rad*db_rad
elif TYPE == 6:
    TAG = "high"# Asymmetric ellipsoid
    KLMS = [1.0, 0, 0, 0, 0, 0, 0, 0, 0,
           0.004444906415833378, 0.003477837670425714, 0.005378876304115907, 0.005783583661113503, -0.002305670990235569, -0.004974910518199871, 0.003351084537350686]
    A_M = 1000
    a = np.sqrt(5/3) * A_M * np.sqrt(1 - 2 * KLMS[8] + 12 * KLMS[4])
    b = np.sqrt(5/3) * A_M * np.sqrt(1 - 2 * KLMS[8] - 12 * KLMS[4])
    c = np.sqrt(5/3) * A_M * np.sqrt(1 + 4 * KLMS[8])
    def indicator(pos):
        return pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2] < (A_M * A_M * 5/3)
        #return pos[0] * pos[0]/(a*a) + pos[1] * pos[1]/(b*b) + pos[2] * pos[2]/(c*c) < 1
else:
    KLMS = None
    A_M = None
    def indicator(pos):
        return None

######################################################

def get_index(lp, mp):
    return lp **2 + mp + lp


def lm_from_index(index):
    l = int(np.sqrt(index))
    m = index - l**2 - l
    return (l, m)

MAX_L = int(np.sqrt(len(KLMS))) - 1

real_hlms = [k for k in KLMS]
real_klms = [k for k in KLMS]
complex_hlms = np.zeros_like(KLMS, dtype=np.cdouble)
complex_klms = np.zeros_like(KLMS, dtype=np.cdouble)
for l in range(MAX_L+1):
    for m in range(-l, l+1):
        if m == 0:
            complex_klms[get_index(l, m)] = KLMS[(l+1)**2 - 1]
        else:
            complex_klms[get_index(l, m)] = KLMS[(l)**2 + 2 * l - 2 * abs(m)] + 1.0j * KLMS[(l)**2 + 2 * l - 2 * abs(m) + 1]
        if m < 0:
            complex_klms[get_index(l, m)] = (-1)**m * complex_klms[get_index(l, m)].conj()
        real_hlms[get_index(l, m)] = real_klms[get_index(l, m)] * A_M**l
        complex_hlms[get_index(l, m)] = complex_klms[get_index(l, m)] * A_M**l

def rlm(l, m, pos):
    r = np.sqrt(max(0, pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]))
    return lpmv(m, l, pos[2] / r) / factorial(l + m) * r**l * np.exp(1j * m * np.arctan2(pos[1], pos[0]))

def get_hlms(l, m, densities):
    hlm = 0
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if not np.isnan(densities[nx,ny,nz]):
                    hlm += densities[nx,ny,nz] * rlm(l, m, [x,y,z])
    return hlm * DIVISION**3

def get_radius(densities):
    rad = 0
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if not np.isnan(densities[nx,ny,nz]):
                    rad += densities[nx,ny,nz] * (x**2 + y**2 + z**2)
    return np.sqrt(rad * DIVISION**3)

def get_mass(densities):
    return np.nansum(densities) * DIVISION**3

if __name__ == "__main__":
    densities = np.full((len(pos_array), len(pos_array), len(pos_array)), np.nan)
    for nx, x in enumerate(pos_array):
        for ny, y in enumerate(pos_array):
            for nz, z in enumerate(pos_array):
                if indicator([x, y, z]):
                    densities[nx, ny, nz] = -rlm(3,0,[x,y,z]).real / A_M**3
    densities /= np.nansum(densities) * DIVISION**3
    radius = get_radius(densities)
    print(radius, '\t', A_M)
    i = 0
    for l in range(4):
        print()
        for m in range(-l, l+1):
            print(get_hlms(l, m, densities) / radius**l, "\t", complex_hlms[i] / A_M**l)
            i += 1
