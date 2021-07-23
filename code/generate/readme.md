# Generate

This code generates data for a simulated asteroid over its flyby. A description of the asteroid should be put in `params.dat`, which should have the following format:
```
Jlmax, Mlmax
J00
J11Re J11Im J10Re
<...other Jlms...>
M00
0 0 0
M22Re 0 0 0 M20Re
M33Re M33Im M32Re M32Im M31Re M31Im M30Re
<...other Mlms...>
spinx spiny spinz
initial roll
impact parameter
velocity
```
All units are given in SI units unless otherwise specified.
Results of the simulated are saved in `params-unresolved.dat` and `params-resolved.dat`. See the readme under the `/code/sim` directory to see what what these files record.

## Use
To generate data, run

`./target/sim <ASTEROID PARAMETERS>`

`ASTEROID PARAMETERS` should be a `.dat` file containing asteroid parameters in the format given above. If no file is provided, `params.dat` will be assumed.
