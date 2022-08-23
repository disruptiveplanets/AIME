import numpy as np

TRUE_CENTER_OF_MASS = [0, 15.839935068775791, 0]

TRUE_MOMENTS = [
    -0.05200629, -0.2021978,
    0, 0,
    0, 0,
    0, 0,
    0
]

UNIFORM_MOMENTS = [
    -0.05195945615840046, -0.20085472167788565,
    0, -0.0003845813086399978,
    0, 0,
    0, -0.0019368958923091067,
] # Uncorrected
ELLIPSOID_AM = 1000
BULK_AM = 1000

samples = np.load("../density/samples/den-cpre-move-1.5-0-samples.npy")

## Compute chi squared agreement between true moments and data and between uniform moments and data and compute the delta chisquared 