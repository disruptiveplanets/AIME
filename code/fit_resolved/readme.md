# Fit

This code generates data for a specific asteroid parameterization, then fits to it with an MCMC. It works by wrapping the C++ simulation code found in `../sim` with `main.cpp` in this directory, then compiles it with `fast-build` to the `.so` found here. The python module `pybind` allows the `.so` to be loaded and called by python, so that we use the `emcee` module to do the fit.

All units are given in SI units unless otherwise specified.

The mcmc checkpoints are saved in `asteroids.h5`, which is excluded from the github. But after running `fit.py`, corner plots and other data can be extracted from `asteroids.h5`.

On the first day of implementing the MCMC, I realized that torque can be written in terms of moment of inertia, which means that the time to run the simulation should be based only on the matrix multiplication and number of integration steps. No dependence on L, n,  should be found; if there is any, it's an implementation mistake. It's a good test to check this.

## Use
Simply run `python fit.py`. The true parameterization to be used is written into the python file.

## To do
* Try for simple asteroid: L, n, m = (1, 1, 1)
* Run on supercomputer
* Try for bigger model

## Timing data
I'm recording (L, n, m), as well as the time on my laptop and the supercomputer time, so that I can track how long this really takes.

|Parameterization id | Local time | Supercomputer time |
|--------------------|------------|--------------------|
| 1 (1, 1, 1)        | 58:45:--   | ?                  |
