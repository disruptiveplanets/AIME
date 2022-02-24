import sys
import numpy as np
sys.path.append("..")
from setup import *

MAX_L_COMPUTE = 3
print(TAG)

mass = 0
radius = 0
klms = np.zeros((MAX_L_COMPUTE + 1)**2, dtype=np.complex)
for nx, x in enumerate(pos_array):
    print(f"{nx} / {len(pos_array)}")
    for ny, y in enumerate(pos_array):
        for nz, z in enumerate(pos_array):
            if indicator([x, y, z]):
                mass += 1
                radius += x*x + y*y + z*z
                lmindex = 0
                for l in range(0, MAX_L_COMPUTE + 1):
                    for m in range(-l, l+1):
                        klms[lmindex] += rlm(l, m, [x, y, z])
                        lmindex += 1

mass *= DIVISION**3
radius *= DIVISION**3 / mass
radius = np.sqrt(radius)
lmindex = 0
for l in range(0, MAX_L_COMPUTE + 1):
    for m in range(-l, l+1):
        klms[lmindex] *= DIVISION**3 / mass / radius**l
        lmindex += 1

lmindex = 0
out_text = ""
for l in range(0, MAX_L_COMPUTE + 1):
    for m in range(-l, l+1):
        out_text += f"K{l}{m}:\t{klms[lmindex]}\n"
        lmindex += 1

print(out_text)