import os

RELOAD = False
tf = open("template-submit.sh", 'r')
template_text = tf.read()[:-1].split('\n')
tf.close()

if RELOAD:
    reload = " reload"
else:
    reload = ''

for fname in os.listdir("../../staged/"):
    fname = os.path.splitext(fname)[0]
    f = open("submit-{}.sh".format(fname), 'w')
    write_text = '\n'.join(template_text[:3])+ '\n' + template_text[3] + \
        fname + '.log\n' + '\n'.join(template_text[4:]) + ' ' + fname + reload

    f.write(write_text)
    f.close()
    os.system("sbatch submit-{}.sh".format(fname))
