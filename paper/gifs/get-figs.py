import os, shutil

ROOT = "/home/jtdinsmo/Dropbox (MIT)/12.420 project/"

f = open(ROOT + "paper/gifs/fig-sources.txt", 'r')
for line in f.readlines():
    if line == '':
        continue
    if line.startswith('#'):
        continue
    args = line[:-1].split(' ')
    if len(args) != 2:
        continue
    name, path = args
    if not os.path.isfile(ROOT + path):
        print("Could not find " + path)
        continue
    if path[-4:] not in [".eps", ".pdf"]:
        print("WARNING: .eps or .pdf file formats are preferred")
    shutil.copyfile(ROOT+path, ROOT + "paper/gifs/"+name)

f.close()
