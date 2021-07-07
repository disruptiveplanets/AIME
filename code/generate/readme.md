# Generate

This code generates data for a simulated asteroid over its flyby. A description of the asteroid should be put in `params.dat`, which should have the following format:
```
L n m
C00
C1-1 C10 C11
<...other Clms...>
rho00 rho01 rho02 rho03 rho04 rho05
<...other densities...>
spin
impact parameter
velocity
central mass
```
All units are given in SI units unless otherwise specified.
Results of the simulated are saved in `params-unresolved.dat` and `params-resolved.dat`. See the readme under the `/code/sim` directory to see what what these files record.

## Use
To generate data, run

`./target/sim <ASTEROID PARAMETERS>`

`ASTEROID PARAMETERS` should be a `.dat` file containing asteroid parameters in the format given above. If no file is provided, `params.dat` will be assumed.
