import os

RELOAD = False
tf = open("template-submit.sh", 'r')
template_text = tf.read()[:-1]
tf.close()
os.mkdir("fits")

if RELOAD:
    reload = " reload"
else:
    reload = ''

for fname in os.listdir("../../staged/"):
    f = open("fits/submit-{}.sh".format(fname), 'w')
    f.write(template_text + fname + reload)
    print(os.system("cd fits; sbatch submit-{}.sh").format(fname))
