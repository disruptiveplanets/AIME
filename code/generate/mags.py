import matplotlib.pyplot as plt

mags = []
phis = []
with open("torque-mags.dat", 'r') as f:
    for line in f.readlines()[:-2]:
        if line == '':
            continue
        smag, sphi = line.split(',')
        mags.append(float(smag))
        phis.append(float(sphi))

plt.plot(phis, mags)
plt.xlabel("phi")
plt.ylabel("Torque magnitude")
plt.show()
