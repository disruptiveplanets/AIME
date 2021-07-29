import asteroids, time, sys

CADENCE = 3600.0

jlms = [5.972e24, 0, 0, 0]
klms = [
    1e6, 1e5, -5e5,
    0, 0, 0, 0, 0, 0, 0 #m3
    ]
spin = [0.00012, 0.00012, 0.00012]
initial_roll = 0
impact_parameter = 31850000
speed = 4000

start = time.time()
try:
    resolved_data = asteroids.simulate(CADENCE, jlms, klms,
        spin[0], spin[1], spin[2], initial_roll, impact_parameter, speed)
except RuntimeError as err:
    print(err)
    sys.exit()
end = time.time()
print("Time taken: {} s".format(end - start))

f = open("resolved.dat", 'w')
for i, dat in enumerate(resolved_data):
    f.write(str(dat) + ' ')
    if i %3 == 2:
        f.write('\n')
f.close()
