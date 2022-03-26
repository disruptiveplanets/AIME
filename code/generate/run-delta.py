import matplotlib.pyplot as plt
import numpy as np
import os, sys
sys.path.append("../fit_resolved")
import random_vector

'''Test how diferent the spin data is when you vary k22 as a function of the direction of initial velocity'''

SIGMA = [0.01, 1e-5]
D_K = 0.00001

os.system("g++ -DASTEROIDS_MAX_J=0 -DASTEROIDS_MAX_K=2 -Wall -std=c++17 ../sim/*.cpp main.cpp -o \"bin/generate\" -O3")

def get_texts(theta, phi, dk=D_K):
    spin_mag = 0.00019392547
    spin=(spin_mag * np.cos(phi) * np.sin(theta), spin_mag * np.sin(phi) * np.sin(theta), spin_mag * np.cos(theta))
    return f"""0 2 # Max j, Max k
1.0 # J00 (must be 1)
1.0 # K00 (must be 1)
0 0 0 # K1 (must be 0)
0.05200629, 0, 0, 0, -0.2029789 # K2 (must be x 0 0 0 y)
1000.0 # radius (remember it's not really the radius)
{spin[0]}, {spin[1]}, {spin[2]}
0.39269908169 # Initial roll
31850000 # Perigee
6000 # Speed
""", f"""0 2 # Max j, Max k
1.0 # J00 (must be 1)
1.0 # K00 (must be 1)
0 0 0 # K1 (must be 0)
{0.05200629}, 0, 0, 0, {-0.2029789-dk} # K2 (must be x 0 0 0 y)
1000.0 # radius (remember it's not really the radius)
{spin[0]}, {spin[1]}, {spin[2]}
{0.39269908169} # Initial roll
31850000 # Perigee
6000 # Speed
"""

def read_spin(f):
    spins = []
    for l in f.readlines():
        if l == '':
            continue
        x,y,z = l.split(' ')
        spins.append((float(x), float(y), float(z)))
    return spins

def get_spins(text):
    with open("2-params.dat",'w') as f:
        f.write(text)
    os.system("./bin/generate 2-params.dat")
    with open("2-params-resolved.dat", 'r') as f:
        return np.array(read_spin(f))

def get_shift(theta, phi):
    t0, t1 = get_texts(theta, phi)
    s0, s1 = get_spins(t0), get_spins(t1)
    like_zero = random_vector.log_likelihood(random_vector.TILT_UNIFORM_TRUE, s0, s0, len(s0), SIGMA)
    like_off = random_vector.log_likelihood(random_vector.TILT_UNIFORM_TRUE, s0, s1, len(s0), SIGMA)
    #diff = get_spins(t0) - get_spins(t1)
    #norms = np.sqrt(diff[:,0]**2 + diff[:,1]**2 + diff[:,2]**2)
    return np.abs(like_off-like_zero)/D_K

thetas = np.linspace(0, np.pi, 20)
phis = np.linspace(0, np.pi, 20)
data = []
for theta in thetas:
    line = []
    for phi in phis:
        line.append(get_shift(theta, phi))
    data.append(line)

print(get_shift(0, 0))
print(get_shift(np.pi, 0))

plt.pcolormesh(phis, thetas, data)
plt.xlabel("Phi")
plt.ylabel("Theta")
plt.colorbar()
plt.savefig("deltas.png")
plt.show()
