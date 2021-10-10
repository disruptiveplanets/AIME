import os
from shutil import copyfile

PATH = "minimizer"

I_LIMIT = -1

for fname in os.listdir("../staged/"):
    fname = os.path.splitext(fname)[0]
    if not os.path.isdir("{1}/{0}".format(fname, PATH)):
        os.mkdir("{1}/{0}".format(fname, PATH))
        i = 0
        while True:
            if i >= I_LIMIT and I_LIMIT > 0: break
            out = os.system("scp jdinsmore@txe1-login.mit.edu:~/asteroid-tidal-torque/code/fit_resolved/{0}-{2}.h5 {1}/{0}/{0}-{2}-all.png".format(fname, PATH, i))
            if out != 0:
                print(out)
                break
            i += 1

    #os.rename("../staged/{0}.dat".format(fname), "{1}/{0}/{0}.dat".format(fname, PATH))
    print(fname)
