import sys
sys.path.append("..")
sys.path.append("../..")
import mcmc_core
import fe, lumpy
from core import Indicator, TrueShape, TrueMoments
import numpy as np

RELOAD_UNC_TRACKER = False

if len(sys.argv) >= 3:
    generate = True if sys.argv[3] == "True" else False
else:
    generate = False
make_map = not generate


MAX_RADIUS = 2000
RUN_NAME = sys.argv[1]
METHOD_NAME = sys.argv[2]
if METHOD_NAME == "lumpy":
    cut_k2m = True
    method_class = lumpy.Lumpy
    method_tag = "lump"
    dof = 2
    if make_map:
        DIVISION = 9
    else:
        DIVISION = 99
    if RUN_NAME == "double":
        lumpy.MODEL_N = 2
        lumpy.SLOW_PRIOR = True
        dof = 7
        mcmc_core.NUM_SUCCESSES = 8
elif METHOD_NAME == "fe":
    cut_k2m = False
    DIVISION = 49
    method_class = fe.FiniteElement
    method_tag = "fe"
    dof = 5
else:
    raise Exception(f"{METHOD_NAME} is not a valid method")
fe.grids.NUM_DRAWS = 0
k22a, k20a = -0.05200629, -0.2021978
k22s, k20s = 0, -0.09766608

ELLIPSOID_AM = 1000

SURFACE_AMS = {
    "asym": 1000,
    "vel": 1000,
    "ori": 1000,
    "sym": 1000,
    "double": 1000,
    "sph-3": 1000,
    "sph-1.5": 1000,
    "move-3": 1002.0081758422925,
    "move-1.5": 1000.1281468600504,
}

BULK_AMS = {
    "asym": 1000,
    "vel": 1000,
    "ori": 1000,
    "sym": 1000,
    "double": 970.4652599064898,
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

blob_rad = 300

TRUE_SHAPES = {
    "asym": TrueShape.uniform(),
    "vel": TrueShape.uniform(),
    "ori": TrueShape.uniform(),
    "sym": TrueShape.uniform(),
    "sph-3": TrueShape.core_sph(3, 500),
    "sph-1.5": TrueShape.core_sph(1.5, 500),
    "move-3": TrueShape.core_shift(3, 500, core_displacement),
    "move-1.5": TrueShape.core_shift(1.5, 500, core_displacement),
    "double": TrueShape.two_core(3, blob_rad, [0, 500, 0], 3, blob_rad, [0, -500, 0]),
}
TRUE_MOMENTS = { # Moments of the known shape
    "asym": TrueMoments.ell(),
    "vel": TrueMoments.ell(),
    "ori": TrueMoments.ell(),
    "sym": TrueMoments.sph(),
    "sph-3": TrueMoments.ell(),
    "sph-1.5": TrueMoments.ell(),
    "move-3": TrueMoments.move_3(),
    "move-1.5": TrueMoments.move_1_5(),
    "double": TrueMoments.ell(),
}
INDICATORS = {
    "asym": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "vel": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "ori": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "sym": Indicator.ell(ELLIPSOID_AM, k22s, k20s),
    "sph-3": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "sph-1.5": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
    "move-3": Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_high),
    "move-1.5": Indicator.ell_y_shift(ELLIPSOID_AM, k22a, k20a, -core_shift_low),
    "double": Indicator.ell(ELLIPSOID_AM, k22a, k20a),
}
SAMPLE_NAMES = {
    "asym": "den-asym",
    "vel": "scale-vel",
    "ori": "ori-e-5",
    "sym": "den-sym",
    "sph-3": "den-core-sph-3",
    "sph-1.5": "den-core-sph-1.5",
    "move-3": "den-core-move-3",
    "move-1.5": "den-core-move-1.5",
    "double": "den-core-double",
}


if RELOAD_UNC_TRACKER:
    unc_tracker_file = f"{RUN_NAME}-{method_tag}-unc-tracker.npy"
else:
    unc_tracker_file = None


asteroid = mcmc_core.MCMCAsteroid(f"{RUN_NAME}-{method_tag}", f"../../samples/{SAMPLE_NAMES[RUN_NAME]}-0-samples.npy", INDICATORS[RUN_NAME],
    TRUE_SHAPES[RUN_NAME], SURFACE_AMS[RUN_NAME], DIVISION, MAX_RADIUS, dof, BULK_AMS[RUN_NAME], TRUE_MOMENTS[RUN_NAME])

result = None
while result is None:
    result = asteroid.pipeline(method_class, make_map, generate=generate, n_samples=1000, unc_tracker_file=unc_tracker_file, cut_k2m=cut_k2m)
    
with open(f"{RUN_NAME}-{method_tag}-unc-tracker.npy", 'wb') as f:
    result.save(f)