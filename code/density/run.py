import sys
from likelihood import Likelihood
from harmonic import Harmonic
from core import Asteroid, Indicator

DIVISION = 99
MAX_RADIUS = 2000
RELOAD = False

k22a, k20a = -0.05200629, -0.2021978
k22s, k20s = 0, -0.09766608


asteroids = {
    "sym-sph": Asteroid("samples/den-sym.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000)),
    "asym-sph": Asteroid("samples/den-asym.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000)),
    "sym-ell": Asteroid("samples/den-sym.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22s, k20s)),
    "asym-ell": Asteroid("samples/den-asym.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a)),
    "tet": Asteroid("samples/den-tet.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.tet(1000)),
    "db": Asteroid("samples/den-db.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.dumbbell(1000)),
    "high": Asteroid("samples/den-high.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000)),
    "in": Asteroid("samples/den-in.npy", 1047.477436728389, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a)),
    "out": Asteroid("samples/den-out.npy", 1050.660629058438, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a)),
    "blob": Asteroid("samples/den-blob.npy", , DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a)),
}

methods = {
    "likelihood": Likelihood,
    "harmonic": Harmonic, 
}

if len(sys.argv) != 3:
    raise Exception("You must provide two arguments")

if sys.argv[1] not in asteroids:
    raise Exception("The first argument must be the asteroid type. One of {}".format(asteroids.keys()))

if sys.argv[2] not in methods:
    raise Exception("The second argument must be the method. One of {}".format(methods.keys()))


method = methods[sys.argv[2]](asteroids[sys.argv[1]])

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