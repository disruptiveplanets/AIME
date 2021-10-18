from matplotlib import pyplot as plt
import numpy as np

t0 = np.array([0, 0])
t1 = np.array([-0.125, -0.25])
t2 = np.array([0.125, -0.25])

SPACE_DIM = 20

def fn(x, y):
    if abs(x) - 1e-3 > -y / 0.25 * 0.125:
        return np.nan
    return x * y

X = []
Y = []
data = []
for i, y in enumerate(np.linspace(-0.25, 0, SPACE_DIM)):
    data_line = []
    X_line = []
    Y_line = []
    delta = 0
    if i % 2 == 1:
        delta = 0.25 / SPACE_DIM / 2
    for x in np.linspace(-0.125+delta, 0.125+delta, SPACE_DIM):
        data_line.append(fn(x, y))
        X_line.append(x)
        Y_line.append(y)
    data.append(data_line)
    X.append(X_line)
    Y.append(Y_line)
data = np.array(data)
X = np.array(X)
Y = np.array(Y)

c = plt.contourf(X, Y, data)
plt.colorbar(c)

plt.xlim(-0.125, 0.125)
plt.ylim(-0.25, 0)

plt.scatter(X, Y)

plt.show()
