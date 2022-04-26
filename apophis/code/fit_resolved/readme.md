## April 24

Made this folder and modified it for a tumbling asteroid. Confirmed that the simulation built, but didn't test anything. Especially in danger are the initializing quaternions (both for the axis conversions from north to ecliptic to inertial) and the rotational periods.

## April 25

Performed a thorough check of the tumbling math, and checked that the true Apophis periods are reproduced. Also checked the ecliptic coordinate system quaternion.

Free parameters: K3m, initial roll & precession, radius, cadence, period gap, observational precision

## April 26

Tested the fit code with the test backend to make sure it runs