import sys, os
from likelihood import Likelihood
from harmonic import Harmonic
from lumpy import Lumpy
from core import Asteroid, Indicator, TrueShape
import matplotlib.pyplot as plt

DIVISION = 99
MAX_RADIUS = 2000
RELOAD = False

k22a, k20a = -0.05200629, -0.2021978
k22s, k20s = 0, -0.09766608


asteroids = {
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
    "in": ("samples/den-in-0-samples.npy", 1047.477436728389, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a), TrueShape.in_(1000, k22a, k20a)),
    "out": ("samples/den-out-0-samples.npy", 1050.660629058438, DIVISION, MAX_RADIUS,
        Indicator.ell(1000, k22a, k20a), TrueShape.out(1000, k22a, k20a)),
    "blob": ("samples/den-blob-0-samples.npy", 894.3680254454393, DIVISION, MAX_RADIUS,
        Indicator.ell_x_shift(1000, k22a, k20a, -148.48101191304406), TrueShape.blob(1000, k22a, k20a)),
}

methods = {
    "likelihood": Likelihood,
    "harmonic": Harmonic, 
    "lumpy": Lumpy, 
}

def make_asteroid(args):
    return Asteroid(args[0], args[1], args[2], args[3], args[4], args[5])
    
if len(sys.argv) == 2:
    if sys.argv[1] != "plot":
        raise Exception("The only one-argument mode is 'plot'")

    for asteroid_name, args in asteroids.items():
        for method_name, method_class in methods.items():
            new_args = [a for a in args]
            new_args[0] = ""

            if not os.path.exists(f"data/{asteroid_name}/{method_name}-d.npy"):
                print(f"Data was not present for asteroid {asteroid_name} method {method_name}")
                continue
            print(f"Plotting asteroid {asteroid_name} method {method_name}")

            asteroid = make_asteroid(new_args)
            method = method_class(asteroid)
            
            method.reload(f"data/{asteroid_name}/{method_name}-d.npy", f"data/{asteroid_name}/{method_name}-u.npy")
            method.display(f"figs/{asteroid_name}/{method_name}")

            plt.close("all")


elif len(sys.argv) == 3:    
    if sys.argv[1] not in asteroids:
        raise Exception("The first argument must be the asteroid type. One of {}".format(asteroids.keys()))

    if sys.argv[2] not in methods:
        raise Exception("The second argument must be the method. One of {}".format(methods.keys()))

    asteroid = make_asteroid(asteroids[sys.argv[1]])
    method = methods[sys.argv[2]](asteroid)

    print("Solving")
    method.solve()
    print("Getting densities")
    method.map_density()
    method.save_density(f"data/{sys.argv[1]}/{sys.argv[2]}-d.npy")
    print("Getting uncertainties")
    method.map_unc()
    method.save_unc(f"data/{sys.argv[1]}/{sys.argv[2]}-u.npy")
    method.check()

    #method.display(f"figs/{sys.argv[1]}/{sys.argv[2]}")

else:
    raise Exception("One or two arguments must be provided.")