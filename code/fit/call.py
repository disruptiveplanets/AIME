import os

EARTH_MASS = 5.972e24

class Asteroid:
    def __init__(self):
        self.L = 1
        self.n = 1
        self.m = 1
        self.clms = [[1], [1, 2, 3]]
        self.rhos = [1, 1, 1, 1, 1, ]
        self.spins = [0.0001, 0.0002, 0.0003]
        self.b = 63700000
        self.v = 4000
        self.m = EARTH_MASS

    def write(self, num):
        # THrow if file exists.
        f = open(str(num) + ".dat", 'r')
        f.write("""{} {} {}
{}
{}
{} {} {}
{}
{}
{}
""".format(self.L, self.n, self.m, ' '.join(self.clms), ' '.join(self.rhos), self.spins[0], self.spins[1], self.spins[2], self.b, self.v, self.m))
        f.close()


def call(asteroid):
    num =
    # Write the asteroid data
    asteroid.write(num)

    # Call the c++ code
    os.system

    # Read c++ code
    f =

    # Delete old file
