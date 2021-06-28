import asteroids, time

CADENCE = 3600.0

L = 1
n = 1
m = 1
clms = [1, 1, 2, 3]
densities = [1, 2, 3, 4, 5, 6]
spinx = 0.0001
spiny = 0.00005
spinz = 0.00004
impact_parameter = 63700000
speed = 4000
central_mass = 5.972e24

start = time.time()
resolved_data = asteroids.simulate(CADENCE, L, n, m, clms, densities, spinx, spiny, spinz,
    impact_parameter, speed, central_mass)
end = time.time()
print("Time taken: {} s".format(end - start))

f = open("resolved.dat", 'w')
for i, dat in enumerate(resolved_data):
    f.write(str(dat) + ' ')
    if i %3 == 2:
        f.write('\n')
f.close()
