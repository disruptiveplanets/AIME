import sys
sys.path.append("..")
sys.path.append("../..")
import main
from core import Indicator

DIVISION = 49
MAX_RADIUS = 2000
RUN_NAME = sys.argv[1]
AST_INDEX = sys.argv[2]
main.grids.NUM_DRAWS = 0

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

k22, k20, surface_am = -0.05200629, -0.2021978, SURFACE_AMS[RUN_NAME] # For the shape

main.pipeline(f"den-core-{RUN_NAME}-{AST_INDEX}", f"../../samples/den-core-{RUN_NAME}-0-samples.npy", Indicator.ell(surface_am, k22, k20),
    surface_am, DIVISION, MAX_RADIUS, False, used_bulk_am=BULK_AMS[RUN_NAME])