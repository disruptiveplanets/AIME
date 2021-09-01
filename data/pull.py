import os
import matplotlib.pyplot as plt
from shutil import copyfile
copyfile("../code/fit_resolved/asteroids.cpython-38-x86_64-linux-gnu.so", "asteroids.cpython-38-x86_64-linux-gnu.so")
import display

PATH = "minimizer"

for fname in os.listdir("../staged/"):
    fname = os.path.splitext(fname)[0]
    if not os.path.isdir("{1}/{0}".format(fname, PATH)):
        os.mkdir("{1}/{0}".format(fname, PATH))
        print(os.system("scp jdinsmore@txe1-login.mit.edu:~/search-first/code/fit_resolved/{0}.h5 {1}/{0}/{0}.h5".format(fname, PATH)))
    os.rename("../staged/{0}.dat".format(fname), "{1}/{0}/{0}.dat".format(fname, PATH))
    print(fname)
    disp = display.Display("{1}/{0}/{0}".format(fname, PATH))
    disp.show_redchi()
    disp.show_params()
    disp.show_corner()
    disp.show_compare()
    disp.show_results()
    plt.show()
    print()
