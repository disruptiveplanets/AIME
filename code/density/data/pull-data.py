import os

names = [
    "sym-sph", 
    "asym-sph", 
    "sym-ell", 
    "asym-ell", 
    "tet", 
    "db", 
    "high", 
    "in", 
    "out", 
    "blob"
]


for fname in names:
    if not os.path.isdir(fname):
        os.mkdir(fname)
    out = os.system("scp jdinsmore@txe1-login.mit.edu:~/asteroid-tidal-torque/code/density/data/{0}/* {0}".format(fname))
    print(fname)
