from scipy.special import lpmv
from math import factorial, atan2, sqrt, acos
from matplotlib import pyplot as plt
import numpy as np
import lmfit
from multiprocessing import Pool

DELTA = 0.01
EPS = 1e-5 # Error bars on fixed values
RADIUS = 1.2
SHAPE = [int(2 * RADIUS / DELTA)] * 3
xs = np.linspace(-RADIUS, RADIUS, SHAPE[0])
MAX_L = 2
RLMS = []
volume_element = DELTA**3

TRUE = np.asarray([
    1,
    0, 0, 0,
    0, 0, 0, 0, -0.1])
ERROR = np.asarray([
    EPS,
    EPS, EPS, EPS,
    0.01, EPS, EPS, EPS, 0.01])

def rlm(coord):
    (l, m) = coord
    this_rlm = np.zeros(SHAPE, dtype=np.complex)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            for k, z in enumerate(xs):
                r = sqrt(x**2 + y**2 + z**2)
                if r > 0:
                    costheta = z / r
                else:
                    costheta = 0
                this_rlm[i, j, k] = lpmv(m, l, costheta) / factorial(l + m) * r**l * np.exp(1j * m * atan2(y, x))
    return this_rlm

def make_r_squared():
    rsq = np.zeros(SHAPE)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            for k, z in enumerate(xs):
                rsq[i, j, k] = x**2 + y**2 + z**2
    return rsq

def make_rlms():
    lms = []
    for l in range(MAX_L+1):
        for m in range(0, l+1, 1):
            lms.append((l, m))
    with Pool() as pool:
        rlms = pool.map(rlm, lms)
    return rlms

def shape_mask(params):
    b = params["b"].value
    c = params["c"].value
    mask = np.zeros(SHAPE)
    for i, x in enumerate(xs):
        for j, y in enumerate(xs):
            for k, z in enumerate(xs):
                ellipse_scale = np.sqrt(x**2 + y**2 / b**2 + z**2 / c**2)
                if ellipse_scale <= 1 - DELTA / 2:
                    mask[i, j, k] = 1
                elif ellipse_scale >= 1 + DELTA / 2:
                    mask[i, j, k] = 0
                else:
                    mask[i, j, k] = (1 - ellipse_scale) / DELTA + 0.5
    return mask

def evaluate(params):
    mask = shape_mask(params)
    klmpred = [np.sum(r * mask) * volume_element for r in RLMS]
    scale_radius = np.sqrt(np.sum(mask * R_SQUARED))
    new_pred = []
    for l in range(MAX_L+1):
        for m in range(l, -1, -1):
            val = klmpred[(l * (l + 1)) // 2 + m]
            new_pred.append(val.real / scale_radius**l)
            if m != 0:
                new_pred.append(val.imag / scale_radius**l)


    # Convert integrals to the actual fit parameters
    new_pred /= new_pred[0]
    return new_pred

def chisq(params):
    return (TRUE - evaluate(params)) / ERROR

def fit():
    params = lmfit.Parameters()
    params.add("b", value=1, min=0, max=RADIUS)
    params.add("c", value=1, min=0, max=RADIUS)
    result = lmfit.minimize(chisq, params)
    print(lmfit.fit_report(result))
    return evaluate(result.params)

def ellipsoid_case():
    # Algebraically compute a, b, c, from moi
    k20 = TRUE[8]
    k22 = TRUE[4]
    Ixx = 2/3.0 * k20 - 4 * k22 + 2/3.0
    Iyy = 2/3.0 * k20 + 4 * k22 + 2/3.0
    Izz = -4/3.0 * k20 + 2/3.0
    asquared = 0.9 * (Iyy + Izz - Ixx)
    bsquared = 0.9 * (Ixx + Izz - Iyy) / asquared
    csquared = 0.9 * (Ixx + Iyy - Izz) / asquared
    return 1, np.sqrt(bsquared), np.sqrt(csquared)

if __name__ == "__main__":
    print("Analytical result")
    print(ellipsoid_case())
    RLMS = make_rlms()
    R_SQUARED = make_r_squared()
    print("Begin fit")
    print(fit())
