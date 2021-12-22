import os
from shutil import copyfile

PATH = "probe-space-1"
FORBIDDEN_NAMES = ["ob"]

I_LIMIT = -1

for fname in os.listdir("../staged/"):
    fname = os.path.splitext(fname)[0]

    forbidden = False
    for n in FORBIDDEN_NAMES:
        if n in fname:
            forbidden = True
            break
    
    if forbidden:
        continue

    if not os.path.isdir("{1}/{0}".format(fname, PATH)):
        os.mkdir("{1}/{0}".format(fname, PATH))
    i = 0
    while True:
        if i >= I_LIMIT and I_LIMIT > 0: break
        out = os.system("scp jdinsmore@txe1-login.mit.edu:~/asteroid-tidal-torque/code/fit_resolved/{0}-{2}-samples.npy {1}/{0}/{0}-{2}-samples.npy".format(fname, PATH, i))
        if out != 0:
            break
        i += 1

    os.rename("../staged/{0}.txt".format(fname), "{1}/{0}/{0}.txt".format(fname, PATH))
    print(fname)
