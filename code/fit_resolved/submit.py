import os, sys

RELOAD = False
if len(sys.argv) == 2:
    if sys.argv[1] == "reload":
        RELOAD = True
    else:
        raise Exception("You passed an invalid second argument.")

tf = open("template-submit", 'r')
template_text = tf.read()[:-1].split('\n')
tf.close()

if RELOAD:
    reload = " reload"
else:
    reload = ''

for fname in os.listdir("../../staged/"):
    fname = os.path.splitext(fname)[0]
    if os.path.exists(fname+".h5") and not RELOAD:
        res = ''
        while res not in ['y', 'n']:
            res = input("The path {} exists. Are you sure you want to overwrite it? (y/n)".format(fname+".h5"))
        if res == 'n':
            continue
    f = open("submit-{}.sh".format(fname), 'w')
    write_text = '\n'.join(template_text[:3])+ '\n' + template_text[3] + \
        fname + '.log\n' + '\n'.join(template_text[4:]) + ' ' + fname + reload

    f.write(write_text)
    f.close()
    os.system("sbatch submit-{}.sh".format(fname))
