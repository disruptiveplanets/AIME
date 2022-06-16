import sys
sys.path.append("..")
sys.path.append("../..")
from mcmc_core import MCMCAsteroid
import fe
from core import Indicator, TrueShape
import numpy as np

DIVISION = 49
MAX_RADIUS = 2000
RUN_NAME = sys.argv[1]
AST_INDEX = sys.argv[2]
fe.grids.NUM_DRAWS = 0
k22a, k20a = -0.05200629, -0.2021978 # Surface
ELLIPSOID_AM = 1000

SURFACE_AMS = {
    "sph-3": 1000,
    "sph-1.5": 1000,
    "move-3": 1002.0081758422925,
    "move-1.5": 1000.1281468600504,
}

BULK_AMS = {
    "sph-3": 922.9234884822591,
    "sph-1.5": 978.4541044108308,
    "move-3": 933.1648422811957,
    "move-1.5": 980.8811439828254,
}

a = np.sqrt(5/3) * ELLIPSOID_AM * np.sqrt(1 - 2 * k20a + 12 * k22a)
b = np.sqrt(5/3) * ELLIPSOID_AM * np.sqrt(1 - 2 * k20a - 12 * k22a)
c = np.sqrt(5/3) * ELLIPSOID_AM * np.sqrt(1 + 4 * k20a)
core_displacement = 300
core_rad = 500
core_vol = np.pi * 4 / 3 * core_rad**3
ellipsoid_vol = np.pi * 4 / 3 * a * b * c
density_factor_low = 0.5
density_factor_high = 2
core_shift_low = core_displacement * (core_vol * density_factor_low) / ellipsoid_vol
core_shift_high = core_displacement * (core_vol * density_factor_high) / ellipsoid_vol

TRUE_SHAPES = {
    "sph-3": TrueShape.core_sph(3, 500),
    "sph-1.5": TrueShape.core_sph(1.5, 500),
    "move-3": TrueShape.core_shift(3, 500, core_displacement),
    "move-1.5": TrueShape.core_shift(1.5, 500, core_displacement),
}
INDICATORS = {
    "sph-3": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "sph-1.5": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "move-3": Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_high),
    "move-1.5": Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_low),
}



asteroid = MCMCAsteroid(f"den-core-{RUN_NAME}-{AST_INDEX}", f"../../samples/den-core-{RUN_NAME}-0-samples.npy", INDICATORS[RUN_NAME],
    TRUE_SHAPES[RUN_NAME], SURFACE_AMS[RUN_NAME], DIVISION, MAX_RADIUS, 9, BULK_AMS[RUN_NAME])

asteroid.pipeline(fe.FiniteElement, False, generate=True)