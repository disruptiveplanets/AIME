import os, sys
from shutil import copyfile

PATH = "cad-speed-contour"
FORBIDDEN_NAMES = []

I_LIMIT = -1

MOVE = False
if len(sys.argv) == 2:
    if sys.argv[-1] == 'move':
        MOVE = True
    else:
        raise Exception("Second argument must be `move`")

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
        out = os.system("scp jdinsmore@txe1-login.mit.edu:~/asteroid-tidal-torque/code/fit_resolved/{0}-{2}-all.png {1}/{0}/{0}-{2}-all.png".format(fname, PATH, i))
        if out != 0:
            break
        i += 1

    if MOVE:
        os.rename("../staged/{0}.txt".format(fname), "{1}/{0}/{0}.txt".format(fname, PATH))
    print(fname)
