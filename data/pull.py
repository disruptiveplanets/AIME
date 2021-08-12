import os
import display
import matplotlib.pyplot as plt

PATH = "converge"

for fname in os.listdir("../staged/"):
    fname = os.path.splitext(fname)[0]
    try:
        os.mkdir("{1}/{0}".format(fname, PATH))
        print(os.system("scp jdinsmore@txe1-login.mit.edu:~/iterative-fit/code/fit_resolved/{0}.h5 {1}/{0}/{0}.h5".format(fname, PATH)))
    except:
        pass
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
