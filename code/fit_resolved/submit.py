import os, sys

tf = open("template-submit", 'r')
template_text = tf.read()[:-1].split('\n')
tf.close()

RELOAD = False

if os.path.exists("errors.dat"):
    os.remove("errors.dat")

def submit(fname):
    if os.path.exists(fname+"-*.h5") and not RELOAD:
        res = ''
        while res not in ['y', 'n']:
            res = input("The path {} exists. Are you sure you want to overwrite it? (y/n)".format(fname+"-*.h5"))
        if res == 'n':
            return
    f = open("{}.sh".format(fname), 'w')
    reload_text = ' reload' if RELOAD else ''
    write_text = '\n'.join(template_text[:3])+ '\n' + template_text[3] + \
        fname + '.log\n' + '\n'.join(template_text[4:]) + ' ' + fname + reload_text

    f.write(write_text)
    f.close()
    os.system("sbatch {}.sh".format(fname))

if __name__ == "__main__":
    if len(sys.argv) == 1:
        for fname in os.listdir("../../staged/"):
            fname = os.path.splitext(fname)[0]
            submit(fname)

    elif len(sys.argv) >= 2:
        for fname in sys.argv[1:]:
            submit(fname)
    else:
        raise Exception("You must pass 0 or 1 arguments.")