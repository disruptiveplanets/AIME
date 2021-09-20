import os
import matplotlib.pyplot as plt
from shutil import copyfile
#copyfile("../code/fit_resolved/asteroids.cpython-38-x86_64-linux-gnu.so", "asteroids.cpython-38-x86_64-linux-gnu.so")
import display

PATH = "minimizer"

for fname in os.listdir("../staged/"):
    fname = os.path.splitext(fname)[0]
    if not os.path.isdir("{1}/{0}".format(fname, PATH)):
        os.mkdir("{1}/{0}".format(fname, PATH))
        i = 0
        while True:
            out = os.system("scp jdinsmore@txe1-login.mit.edu:~/asteroid-tidal-torque/code/fit_resolved/{0}-{2}.h5 {1}/{0}/{0}-{2}.h5".format(fname, PATH, i))
            if out != "":
                print(out)
                break
            i += 1

    os.rename("../staged/{0}.dat".format(fname), "{1}/{0}/{0}.dat".format(fname, PATH))
    print(fname)
    i = 0
    while True:
        try:
            disp = display.Display("{1}/{0}/{0}".format(fname, PATH), "{1}/{0}/{0}-{2}".format(fname, PATH, i))
        except Exception as e:
            if i == 0:
                raise e
            break
        disp.show_redchi()
        disp.show_params()
        disp.show_corner()
        disp.show_compare()
        disp.show_results()
        plt.show()
        print()
        i += 1
