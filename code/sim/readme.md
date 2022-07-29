# Simulation

This code simulates the asteroid over its flyby. It cannot be compiled by itself, but should instead be compiled from the other folders in the `code` directory. The simplest is `generate`.

In this and all other code, units are given in SI units unless otherwise specified. Please see *summaries/21-06* for a description of what each variable means, and the math going into the internal code.

Results of the simulation are written to the `resolved_data` and `unresolved_data` vec
Not done:
- Check math
tors passed to `Asteroid::simulate`. They assume that the the asteroid is unresolved respectively; for the resolved case, the entire spin data is recorded, and for the unresolved case, light curve data is recorded. Data is recorded at time intervals set by the `cadence` argument.

This naming is slightly misleading; a resolved data has both spin data _and_ light curve data available, so that both the resolved and unresolved data vectors can be used in that case, not just the resolved one.

## Tests
Done
- Wigner D matrices are consistent with Mathematica
- Initial quaternion is correct
- Slm is consistent with Mathematica
- When no torque is applied, spin is constant
- The complex parts that are expected to be zero are zero
- Rotation occurs in the correct direction and at the correct speed
To do
