import sys
from likelihood import Likelihood
from harmonic import Harmonic
from lumpy import Lumpy
from core import Asteroid, Indicator

DIVISION = 99
MAX_RADIUS = 2000
RELOAD = False

k22a, k20a = -0.05200629, -0.2021978
k22s, k20s = 0, -0.09766608


asteroids = {
    "sym-sph": ("samples/den-sym-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000)),
    "asym-sph": ("samples/den-asym-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000)),
    "sym-ell": ("samples/den-sym-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22s, k20s)),
    "asym-ell": ("samples/den-asym-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a)),
    "tet": ("samples/den-tet-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.tet(1000)),
    "db": ("samples/den-db-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.dumbbell(1000)),
    "high": ("samples/den-high-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000)),
    "in": ("samples/den-in-0-samples.npy", 1047.477436728389, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a)),
    "out": ("samples/den-out-0-samples.npy", 1050.660629058438, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a)),
    "blob": ("samples/den-blob-0-samples.npy", 894.3680254454393, DIVISION, MAX_RADIUS,
        Indicator.ell_x_shift(1000, k22a, k20a, -148.48101191304406)),
}

methods = {
    "likelihood": Likelihood,
    "harmonic": Harmonic, 
    "lumpy": Lumpy, 
}

if len(sys.argv) != 3:
    raise Exception("You must provide two arguments")

if sys.argv[1] not in asteroids:
    raise Exception("The first argument must be the asteroid type. One of {}".format(asteroids.keys()))

if sys.argv[2] not in methods:
    raise Exception("The second argument must be the method. One of {}".format(methods.keys()))

def make_asteroid(args):
    return Asteroid(args[0], args[1], args[2], args[3], args[4])

asteroid = make_asteroid(asteroids[sys.argv[1]])
method = methods[sys.argv[2]](asteroid)

if not RELOAD:
    print("Solving")
    method.solve()
    print("Getting densities")
    method.map_density()
    method.save_density(f"data/{sys.argv[1]}/{sys.argv[2]}-d.npy")
    print("Getting uncertainties")
    method.map_unc()
    method.save_unc(f"data/{sys.argv[1]}/{sys.argv[2]}-u.npy")
    method.check()
else:
    method.reload(f"data/{sys.argv[1]}/{sys.argv[2]}-d.npy", f"data/{sys.argv[1]}/{sys.argv[2]}-u.npy")

print("Plotting")
method.display(f"figs/{sys.argv[1]}/{sys.argv[2]}")