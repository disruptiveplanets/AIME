import sys
sys.path.append("..")
sys.path.append("../..")
from mcmc_core import MCMCAsteroid
import fe, lumpy
from core import Indicator, TrueShape, TrueMoments
import numpy as np

DIVISION = 49
MAX_RADIUS = 2000
RUN_NAME = sys.argv[1]
METHOD_NAME = sys.argv[2]
if METHOD_NAME == "lumpy":
    method_class = lumpy.Lumpy
    method_tag = "lump"
    dof = 2
elif METHOD_NAME == "fe":
    method_class = fe.FiniteElement
    method_tag = "fe"
    dof = 9
else:
    raise Exception(f"{METHOD_NAME} is not a valid method")
fe.grids.NUM_DRAWS = 0
k22a, k20a = -0.05200629, -0.2021978 # Surface
ELLIPSOID_AM = 1000

SURFACE_AMS = {
    "asym-ell": 1000,
    "sph-3": 1000,
    "sph-1.5": 1000,
    "move-3": 1002.0081758422925,
    "move-1.5": 1000.1281468600504,
}

BULK_AMS = {
    "asym-ell": 1000,
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
    "asym-ell": TrueShape.uniform(),
    "sph-3": TrueShape.core_sph(3, 500),
    "sph-1.5": TrueShape.core_sph(1.5, 500),
    "move-3": TrueShape.core_shift(3, 500, core_displacement),
    "move-1.5": TrueShape.core_shift(1.5, 500, core_displacement),
}
TRUE_MOMENTS = { # Moments of the known shape
    "asym-ell": TrueMoments.ell(),
    "sph-3": TrueMoments.ell(),
    "sph-1.5": TrueMoments.ell(),
    "move-3": TrueMoments.move_3(),
    "move-1.5": TrueMoments.move_1_5(),
}
INDICATORS = {
    "asym-ell": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "sph-3": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "sph-1.5": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "move-3": Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_high),
    "move-1.5": Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_low),
}
SAMPLE_NAMES = {
    "asym-ell": "den-asym",
    "sph-3": "den-core-sph-3",
    "sph-1.5": "den-core-sph-1.5",
    "move-3": "den-core-move-3",
    "move-1.5": "den-core-move-1.5",
}



asteroid = MCMCAsteroid(f"{RUN_NAME}-{method_tag}", f"../../samples/{SAMPLE_NAMES[RUN_NAME]}-0-samples.npy", INDICATORS[RUN_NAME],
    TRUE_SHAPES[RUN_NAME], SURFACE_AMS[RUN_NAME], DIVISION, MAX_RADIUS, dof, BULK_AMS[RUN_NAME], TRUE_MOMENTS[RUN_NAME])

asteroid.pipeline(method_class, True, generate=True)