# Generate

This code generates data for a simulated asteroid over its flyby. A description of the asteroid should be put in `params.dat`, which should have the following format:
```
Jlmax, Mlmax
J00
J11Re J11Im J00Re
<...other Jlms...>
M22Re M22Im M21Re M21Im M20Im
<...other Mlms...>
spin_mag
impact parameter
velocity
```
All units are given in SI units unless otherwise specified.
Results of the simulated are saved in `params-unresolved.dat` and `params-resolved.dat`. See the readme under the `/code/sim` directory to see what what these files record.

## Use
To generate data, run

`./target/sim <ASTEROID PARAMETERS>`

`ASTEROID PARAMETERS` should be a `.dat` file containing asteroid parameters in the format given above. If no file is provided, `params.dat` will be assumed.
