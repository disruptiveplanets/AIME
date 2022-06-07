import sys
sys.path.append("..")
sys.path.append("../..")
import main
from core import Indicator

DIVISION = 49
MAX_RADIUS = 2000
AST_INDEX = sys.argv[1]
main.grids.NUM_DRAWS = 0

k22, k20, surface_am = -0.05200629, -0.2021978, 1000 # For the shape

main.pipeline(f"den-core-sph-{AST_INDEX}", "../../samples/den-core-sph-0-samples.npy", Indicator.ell(surface_am, k22, k20),
    surface_am, DIVISION, MAX_RADIUS, False, used_bulk_am=978.4541044108308)