import matplotlib.pyplot as plt
import numpy as np
import os

alphas = np.linspace(-np.pi, np.pi, 50)
spins = 2 * np.pi / (9 * 3600) * np.array([np.cos(alphas), np.sin(alphas), 0*alphas]).transpose()

with open('2-params-r.dat','r') as f:
    record_text = f.readlines()
with open('2-params-rp.dat','r') as f:
    record_p_text = f.readlines()

def replace(text, spin):
    lines = [l for l in text]
    lines[6] = ', '.join([str(s) for s in spin])+"\n"
    return ''.join(lines)

def get_spins(filename):
    sx = []
    sy = []
    sz = []
    with open(filename, 'r') as f:
        for line in f.readlines()[:-1]:
            sx_, sy_, sz_ = line.split(' ')
            sx.append(float(sx_))
            sy.append(float(sy_))
            sz.append(float(sz_))
    return np.array([sx, sy, sz])

resids = []
for spin in spins:
    with open('2-params-r.dat','w') as f:
        f.write(replace(record_text, spin))
    with open('2-params-rp.dat','w') as f:
        f.write(replace(record_p_text, spin))
    
    os.system("./bin/generate 2-params-r.dat")
    os.system("./bin/generate 2-params-rp.dat")
    normal=get_spins("2-params-r-resolved.dat")
    perfect=get_spins("2-params-rp-resolved.dat")
    resid_norm = np.sum((normal - perfect)**2)
    resids.append(resid_norm)

plt.plot(alphas, resids)
plt.xlabel("Initial alpha")
plt.ylabel("Residual squared sum")
plt.savefig("resids.png")
plt.show()