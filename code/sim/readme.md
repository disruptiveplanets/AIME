# Simulation

This code simulates the asteroid over its flyby. It cannot be compiled by itself, but should instead be compiled from the other folders in the `code` directory. The simplest is `generate`.

In this and all other code, units are given in SI units unless otherwise specified. Please see *summaries/21-06* for a description of what each variable means, and the math going into the internal code.

Results of the simulation are written to the `resolved_data` and `unresolved_data` vectors passed to `Asteroid::simulate`. They assume that the the asteroid is unresolved respectively; for the resolved case, the entire spin data is recorded, and for the unresolved case, light curve data is recorded. Data is recorded at time intervals set by the `cadence` argument.

This naming is slightly misleading; a resolved data has both spin data _and_ light curve data available, so that both the resolved and unresolved data vectors can be used in that case, not just the resolved one.

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
- The Euler angles are correct.
- The D matrix calculation is correct
- Spherical harmonics are rotated correctly

Not done:
* How does mass come into play in the MOI definition?
* MOI rotates like a matrix
* All the numbers I assumed are real are actually real
* New torque matches old torque
