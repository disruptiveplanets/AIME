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

Results of the simulated are saved in `params-unresolved.dat`, which assumes that the the asteroid is unresolved (and only light curve data can be captured) and `params-resolved.dat`, which assumes that the asteroid is resolved (and full rotation data can be captured) For the resolved case, the results are the spin vector (x, y, z in global coordinates) at periodic intervals. By default, the period (`CADENCE`, defined in main.cpp), is one hour.

## Use
To generate data, run

`./target/sim <ASTEROID PARAMETERS>`

`ASTEROID PARAMETERS` should be a `.dat` file containing asteroid parameters in the format given above. If no file is provided, `params.dat` will be assumed.

To draw the asteroid, run

`./target/sim <ASTEROID PARAMETERS> draw <DIRECTION>`

`ASTEROID PARAMETERS` fulfills the same role as in the first paragraph, except that it is required here. `DIRECTION` should be either `x`, `y`, or `z`. It represents the axis that the asteroid is projected onto. Non-axial directions would be easy to implement, but are not yet simply because there has not yet been a use for them. If `DIRECTION` is not provided, `z` is assumed.

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
- The object rotates the proper number of times per second
- Nutation is stable
- Torque is not generated on a spherically symmetric object
- Ang mom increases for steady torque.
- Results are not sensitive to changes in integration time step, or the way that the integration time step is lengthened at the edges of the orbit

## Features
Implemented
- The ability to read asteroid data from a text file and generate an asteroid model
- Code to automatically calculate the asteroid's mass, center of mass, moment of inertia matrix, and tidal torque
- A program to plot the shape of the asteroid as projected onto various planes with Painter's algorithm and `matplotlib`
- A simulation of the asteroid's flight through the sphere of influence of a planet. It runs to my satisfaction in **16.5 seconds**

## To do:
- Fix the fact that my build command is wrong
- Fix up the drawing code: it doesn't free asteroid and I might not have finished porting it.
- Test to make sure I haven't caused any bugs with the new port
    - I.e., test that the output is the same.
- Test the speed: make sure it's comparable
- Implement emcee! 
