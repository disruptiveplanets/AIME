import numpy as np
import matplotlib.pyplot as plt
import corner
import os, sys

PULL = True

if PULL:
    os.system("scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/orientation/fit_resolved/ori-0-samples.npy .")
    os.system("scp jdinsmore@txe1-login.mit.edu:asteroid-tidal-torque/code/fit_resolved/scan-cadence/cad-00-0-samples.npy .")

with open("cad-00-0-samples.npy", 'rb') as f:
    vel_moments = np.load(f)

with open("ori-0-samples.npy", 'rb') as f:
    o_moments = np.load(f)

print(vel_moments.shape)
print(o_moments.shape)
sys.exit()

vel_resids = vel_moments - np.mean(vel_moments, axis=0)
o_resids = o_moments - np.mean(o_moments, axis=0)
ratio = np.mean(np.std(vel_resids, axis=0)) / np.mean(np.std(o_resids, axis=0))
print(f"Ratio {ratio}")
vel_moments = vel_resids / ratio + np.mean(vel_moments, axis=0)

fig = corner.corner(vel_moments)
corner.corner(o_resids, fig=fig)
plt.show()