import sys
from likelihood import Likelihood
from harmonic import Harmonic
from core import Asteroid, Indicator

DIVISION = 99
MAX_RADIUS = 1500
RELOAD = False


asteroids = {
    "sym": Asteroid("samples/param-024-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000))
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
    method.save_density(f"data/{sys.argv[1]}-{sys.argv[2]}-d.npy")
    print("Getting uncertainties")
    method.map_unc()
    method.save_unc(f"data/{sys.argv[1]}-{sys.argv[2]}-u.npy")
    method.check()
else:
    method.reload(f"data/{sys.argv[1]}-{sys.argv[2]}-d.npy", f"data/{sys.argv[1]}-{sys.argv[2]}-u.npy")

print("Plotting")
method.display(f"figs/{sys.argv[1]}-{sys.argv[2]}")