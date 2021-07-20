import asteroids, time

CADENCE = 3600.0

jlms = [5.972e24, 5.972e22, -5.972e22, 4.972e22]
mlms = [
    1, 2, 3, 4, 5,#m2
    0, 0, 0, 0, 0, 0, 0 #m3
    ]
spin = 0.00012
impact_parameter = 31850000
speed = 4000

start = time.time()
resolved_data = asteroids.simulate(CADENCE, jlms, mlms, spin,
    impact_parameter, speed)
end = time.time()
print("Time taken: {} s".format(end - start))

f = open("resolved.dat", 'w')
for i, dat in enumerate(resolved_data):
    f.write(str(dat) + ' ')
    if i %3 == 2:
        f.write('\n')
f.close()
