## April 24

Made this folder and modified it for a tumbling asteroid. Confirmed that the simulation built, but didn't test anything. Especially in danger are the initializing quaternions (both for the axis conversions from north to ecliptic to inertial) and the rotational periods.

## April 25

Performed a thorough check of the tumbling math, and checked that the true Apophis periods are reproduced. Also checked the ecliptic coordinate system quaternion.

Free parameters: K3m, initial roll & precession, radius, cadence, period gap, observational precision

## April 26

Tested the fit code with the test backend to make sure it runs

## April 28

Fixed a bug with the initial orientation (it's now z-y-z and confirmed with as long of a fourier transform as I can make with 4 GB Memory)..

## April 30

To test the degeneracy, I edited main.cpp so that the initial roll is always pi/4. The redchi and log-like are adjusted to make the first parameter still meaningful; that way the fit still converges. This requires the use of fit-less and recompiling main.cpp