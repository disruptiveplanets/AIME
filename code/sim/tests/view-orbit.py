from matplotlib import pyplot as plt

xs = []
ys = []
zs = []
f = open("orbit-path.txt", 'r')
start=False
for line in f.readlines():
    if line == '\n':
        start=True
        continue
    if line == '':
        continue
    if line[0] == 'S':# End of text (Simulation ended message)
        continue
    if start:
        x, y, z = line.split(' ')
        xs.append(float(x))
        ys.append(float(y))
        zs.append(float(z))

plt.figure(figsize=(6, 6))
plt.axis('equal')
plt.xlabel("y (m)")
plt.ylabel("z (m)")
plt.scatter(ys, zs, color='k', marker='.')
plt.plot(ys, zs, color='k')
plt.scatter(0, 0, color='C1')
plt.savefig("img-orbit-path.png")
plt.show()
