# Goal: spawn a bunch of runs that cover the parameter space.
from matplotlib import pyplot as plt
import numpy as np

EARTH_RADIUS = 6378000
EARTH_MASS = 3.986004418e14

JUPITER_RADIUS = 71492000
JUPITER_MASS = 1.26686534e17

ORIGINAL_PERIGEE = 5 * EARTH_RADIUS
ORIGINAL_VEL = 6000


original_a = -EARTH_MASS / ORIGINAL_VEL**2
eccentricity = -(ORIGINAL_PERIGEE + np.sqrt(ORIGINAL_PERIGEE**2 - 4 * original_a * (ORIGINAL_PERIGEE - original_a))) / (2 * original_a)

new_perigee = ORIGINAL_PERIGEE * JUPITER_RADIUS / EARTH_RADIUS
new_vel = np.sqrt(JUPITER_MASS / new_perigee * (eccentricity - 1))

print(f"Old perigee {ORIGINAL_PERIGEE}\t old velocity {ORIGINAL_VEL}")
print(f"Old a {original_a}\t old e {eccentricity}")
print(f"New perigee {new_perigee}\t new velocity {new_vel}")
print(f"New a {-JUPITER_MASS / new_vel**2}\t new e {eccentricity}")


def get_text(vel):
    return f"""0, 3
120
{new_perigee / EARTH_RADIUS}
1000
{vel}
0.00006464182, 0.00012928364, -0.00012928364
1.0
0.39269908169, 0.05200629, -0.2021978, 0, 0, 0, 0, 0, 0, 0
0.78539816339, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
-0.78539816339, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5"""


f = open("../../staged/jupiter-same-orbit.txt", 'w')
f.write(get_text(new_vel))
f.close()

f = open("../../staged/jupiter-same-vel.txt", 'w')
f.write(get_text(ORIGINAL_VEL))
f.close()
