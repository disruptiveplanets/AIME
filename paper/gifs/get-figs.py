import os, shutil

ROOT = "/Users/jtd/Documents/research/12.420 project/"

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
    if path[-4:] not in [".gif"]:
        print("WARNING: .gif file formats are preferred")
    shutil.copyfile(ROOT+path, ROOT + "paper/gifs/"+name)

    bare = ROOT + "paper/gifs/"+name[:-4]
    print(bare)
    os.system(f"ffmpeg -f gif -i \"{bare}.gif\" \"{bare}.mp4\"")
    os.remove(f"{bare}.gif")

f.close()
