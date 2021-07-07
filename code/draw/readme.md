# Draw

This code draws an asteroid's shape and density distribution as projected onto a given axis. A description of the asteroid should be put in `params.dat`, which should have the following format:
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
The output is given in `params_<axis>.pdf`, where `<axis>` is either `1`, `2`, or `3`, depending on the user's choice. 1 indicates the principal axis with the largest eigenvalue, 3 with the smallest.

## Use
To draw the asteroid, run

`./target/sim <ASTEROID PARAMETERS> <DIRECTION>`

`ASTEROID PARAMETERS` is the location of the parameters file (include the name of the file itself). It is a required argument. `DIRECTION` should be either `1`, `2`, or `3`. It represents the principal axis that the asteroid is projected onto. Non-axial directions would be easy to implement, but are not yet simply because there has not yet been a use for them. If `DIRECTION` is not provided, `1` is assumed.
