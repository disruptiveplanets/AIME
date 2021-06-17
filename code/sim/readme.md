# Simulation

This code simulates the asteroid over its flyby. A description of the asteroid should be put in `params.dat`, which should have the following format:
```
L n m
C00
C1-1 C10 C11
<...other Clms...>
rho00 rho01 rho02 rho03 rho04 rho05
<...other densities...>
spinx, spiny, spinz
impact parameter
velocity
central mass
```
All units are given in SI units unless otherwise specified. Please see *summaries/21-06* for a description of what each variable means, and the math going into the internal code.

Results of the simulated are saved in `params-unresolved.dat`, which assumes that the the asteroid is unresolved (and only light curve data can be captured) and `params-resolved.dat`, which assumes that the asteroid is resolved (and full rotation data can be captured)

## Tests
Done
- The integrals calculated in `math/polyhedron.nb` are accurately reproduced in the code
- The mass, center of mass, and moment of inertia of a cube are correct.
- For large `n` and only Clm00 is nonzero, the mass and moment of inertia reduce to that of a sphere. (Done for `n=20`)
- The orbits follow hyperbolae
- Neglecting torque and with zero initial angular momentum, the object does not change orientation
- Orbit is not sensitive to changes in integration time step, or the way that the integration time step is lengthened at the edges of the orbit
- Neglecting torque and with nonzero initial angular momentum, the object spins at a constant rate
- Neglecting torque, angular momentum is conserved
- Neglecting torque, the object does not change orientation

Not done
- The object rotates the proper number of times per second
- Torque is not generated on a spherically symmetric object
- Nutation is stable
- Results are not sensitive to changes in integration time step, or the way that the integration time step is lengthened at the edges of the orbit
- What is the best dt(r, omega) profile that leads to the most stable
    - orbital path?
    - rotational angular momentum conservation?
    - orientation conservation in the case of no initial spin?

## Features
Implemented
- The ability to read asteroid data from a text file and generate an asteroid model
- Code to automatically calculate the asteroid's mass, center of mass, moment of inertia matrix, and tidal torque
- A program to plot the shape of the asteroid as projected onto various planes with Painter's algorithm and `matplotlib`

Not yet implemented
- A simulation of the asteroid's flight through the sphere of influence of a planet
- An MCMC to regress for asteroid shape.
