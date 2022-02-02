# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np
import os
from scipy.spatial import Delaunay
from scipy.linalg import norm

TRUE_THETA = np.array([0.39269908169, 0.05200629, -0.2021978, 0, 0, 0, 0, 0, 0, 0])
DIST_DIMENSION = 1

plt.style.use("jcap")

file_names = os.listdir()
file_names.sort()

pf = open("percentiles.dat", 'w')

for bare in file_names:
    if not os.path.isdir(bare):
        continue

    min_dist = None
    min_text = None
    for file in os.listdir(bare):
        if not file[-11:] == "samples.npy": continue
        with open("{}/{}".format(bare, file), 'rb') as f:
            array = np.load(f)
        
        flat_samples = array.reshape(-1, array.shape[-1])

        # Get means
        means = np.mean(flat_samples, axis=0)
        dist = (means - TRUE_THETA)[DIST_DIMENSION]
        if min_dist is None or dist < min_dist:
            data = []
            for i in range(flat_samples.shape[1]):
                data.append(", ".join([
                    str(np.mean(flat_samples[:,i])),
                    str(np.percentile(flat_samples[:,i], 50 + 95.45 / 2)),
                    str(np.percentile(flat_samples[:,i], 50 + 68.27 / 2)),
                    str(np.percentile(flat_samples[:,i], 50)),
                    str(np.percentile(flat_samples[:,i], 50 - 68.27 / 2)),
                    str(np.percentile(flat_samples[:,i], 50 - 95.45 / 2)),
                ]))

            min_text = file + ": " +  ": ".join(data)+"\n"
            min_dist = dist

    #if min_text is None:
    #    continue
    pf.write(min_text)
    print(bare, min_dist)

pf.close()


lons = np.linspace(-np.pi, np.pi, 180)
lats = np.linspace(-np.pi/2, np.pi/2, 90)
xyzs = []

for dir_name in file_names:
    if not os.path.isdir(dir_name):
        continue
    with open(f"{dir_name}/{dir_name}.txt", 'r') as f:
        max_j, max_l = f.readline().split(", ")
        max_j, max_l = (int(max_j), int(max_l))
        cadence = int(f.readline())
        perigee = int(f.readline())
        radius = float(f.readline())
        speed = float(f.readline())
        spin = [float(x) for x in f.readline().split(',')]
        jlms = [float(x) for x in f.readline().split(',')]
        theta_true = [float(x) for x in f.readline().split(',')]
        theta_high = [float(x) for x in f.readline().split(',')]
        theta_low = [float(x) for x in f.readline().split(',')]
        sigma = [float(d) for d in f.readline().split(',')]
    xyzs.append(np.array(spin) / norm(spin))
xyzs = np.array(xyzs)
tris = Delaunay(xyzs)

def get_tri(p):
    num_admitted = 0
    min_ts = None
    for ts in tris.convex_hull:
        pts = [xyzs[t] for t in ts]
        normals = [np.cross(pt, p) for pt in pts]
        dots = [[np.dot(normals[i], pts[j]) for i in range(3)] for j in range(3)]
        if dots[0][1] <= dots[0][0] <= dots[0][2] or dots[0][2] <= dots[0][0] <= dots[0][1]:
            if dots[1][0] <= dots[1][1] <= dots[1][2] or dots[1][2] <= dots[1][1] <= dots[1][1]:
                if dots[2][0] <= dots[2][2] <= dots[2][1] or dots[2][1] <= dots[2][2] <= dots[2][0]:
                    big_normal = np.cross(pts[1] - pts[0], pts[2] - pts[0])
                    my_dot = np.dot(big_normal, p)
                    plane_dot = np.dot(big_normal, pts[0])
                    if not ((my_dot < 0) ^ (plane_dot < 0)):
                        num_admitted += 1
                        min_ts = ts
    return min_ts

adjusted_cart_array_data = []
for lat in lats:
    print(lat)
    cart_line = []
    for lon in lons:
        shrink = 0.999
        p = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
        tri = get_tri(p)
        normal = np.cross(xyzs[tri[1]] - xyzs[tri[0]], xyzs[tri[2]] - xyzs[tri[0]])
        p *= np.dot(normal, xyzs[tri[0]]) / np.dot(normal, p) * shrink
        cart_line.append(p)
    adjusted_cart_array_data.append(cart_line)

with open("projecteds.dat", 'wb') as f:
    np.save(f, adjusted_cart_array_data)