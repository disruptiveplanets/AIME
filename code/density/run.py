import sys, os
from likelihood import Likelihood
from harmonic import Harmonic
from lumpy import Lumpy
from core import Asteroid, Indicator, TrueShape
import matplotlib.pyplot as plt

DIVISION = 9
MAX_RADIUS = 2000
RELOAD = False

k22a, k20a = -0.05200629, -0.2021978
k22s, k20s = 0, -0.09766608


asteroids = {
    "test": ("samples/den-in-0-samples.npy", 1047.477436728389, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a), TrueShape.in_(1000)),
    "sym-sph": ("samples/den-sym-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000), None),
    "asym-sph": ("samples/den-asym-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000), None),
    "sym-ell": ("samples/den-sym-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22s, k20s), TrueShape.uniform()),
    "asym-ell": ("samples/den-asym-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a), TrueShape.uniform()),
    "tet": ("samples/den-tet-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.tet(1000), TrueShape.uniform()),
    "db": ("samples/den-db-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.dumbbell(1000), TrueShape.uniform()),
    "high": ("samples/den-high-0-samples.npy", 1000, DIVISION, MAX_RADIUS,
        Indicator.sph(1000), TrueShape.uniform()),
    "in": ("samples/den-in-0-samples.npy", 809.2996416062074, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a), TrueShape.in_(1000)),
    "out": ("samples/den-out-0-samples.npy", 1264.373369575561, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a), TrueShape.out(1000)),
    "in-sph": ("samples/den-tet-0-samples.npy", 809.2996416062074, DIVISION, MAX_RADIUS,
        Indicator.sph(1000), TrueShape.in_sph(1000)),
    "out-sph": ("samples/den-tet-0-samples.npy", 1264.373369575561, DIVISION, MAX_RADIUS,
        Indicator.sph(1000), TrueShape.out_sph(1000)),
    "blob": ("samples/den-blob-0-samples.npy", 965.2268359486612, DIVISION, MAX_RADIUS,
        Indicator.ell_y_shift(1000, k22a, k20a, -57.02376624759285), TrueShape.blob(1000, k22a, k20a)),
    "rot-blob": ("samples/den-blob-0-samples.npy", 966.4144167196462, DIVISION, MAX_RADIUS,
        Indicator.ell_y_shift(1000, k22a, k20a, -57.02376624759285), TrueShape.rot_blob(1000, k22a, k20a)),
}

methods = {
    "likelihood": Likelihood,
    "harmonic": Harmonic, 
    "lumpy": Lumpy, 
}

def make_asteroid(args):
    return Asteroid(args[0], args[1], args[2], args[3], args[4], args[5], args[6])
    
if len(sys.argv) == 2:
    if sys.argv[1] != "plot":
        raise Exception("The only one-argument mode is 'plot'")

    for asteroid_name, args in asteroids.items():
        for method_name, method_class in methods.items():
            new_args = [a for a in args]
            #new_args[0] = ""

            if not os.path.exists(f"data/{asteroid_name}/{method_name}-d.npy"):
                print(f"Data was not present for asteroid {asteroid_name} method {method_name}")
                continue
            print(f"Plotting asteroid {asteroid_name} method {method_name}")

            asteroid = make_asteroid([asteroid_name]+new_args)
            method = method_class(asteroid)
            
            method.reload(f"data/{asteroid_name}/{method_name}-d.npy", f"data/{asteroid_name}/{method_name}-u.npy")
            method.display()

            plt.close("all")


elif len(sys.argv) == 3:    
    if sys.argv[1] not in asteroids:
        raise Exception("The first argument must be the asteroid type. One of {}".format(asteroids.keys()))

    if sys.argv[2] not in methods:
        raise Exception("The second argument must be the method. One of {}".format(methods.keys()))

    asteroid = make_asteroid([sys.argv[1]] + list(asteroids[sys.argv[1]]))
    method = methods[sys.argv[2]](asteroid)

    print("Solving")
    method.solve()
    print("Getting densities")
    method.map_density()
    method.save_density()
    print("Getting uncertainties")
    method.map_unc()
    method.save_unc()
    method.check()

    method.display()

else:
    raise Exception("One or two arguments must be provided.")