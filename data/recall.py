import os
from shutil import copyfile

PATH = "probe"

for fname in os.listdir(PATH):
    if not os.path.isdir("{}/{}".format(PATH, fname)): continue
    try:
        os.rename("{0}/{1}/{1}.txt".format(PATH, fname), "../staged/{0}.txt".format(fname))
        print("Success")
    except:
        pass
