import os

RELOAD = False
tf = open("template-submit.sh", 'r')
template_text = tf.read()[:-1]
tf.close()

if RELOAD:
    reload = " reload"
else:
    reload = ''

for fname in os.listdir("../../staged/"):
    fname = os.path.splitext(fname)[0]
    f = open("submit-{}.sh".format(fname), 'w')
    f.write(template_text + fname + reload)
    f.close()
    os.system("sbatch submit-{}.sh".format(fname))
