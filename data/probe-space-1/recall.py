import os
from shutil import copyfile

for fname in os.listdir():
    if not 'ob' in fname:
        continue
    if not os.path.isdir(fname): continue
    try:
        os.rename("{0}/{0}.txt".format(fname), "../../staged/{0}.txt".format(fname))
        print("Success")
    except:
        pass
