# Fit

This code generates data for a specific asteroid parameterization, then fits to it with an MCMC. It works by wrapping the C++ simulation code found in `../sim` with `main.cpp` in this directory, then compiles it with `fast-build` to the `.so` found here. The python module `pybind` allows the `.so` to be loaded and called by python, so that we use the `emcee` module to do the fit.

All units are given in SI units unless otherwise specified.

## Use
Simply run `python fit.py <params>`, where `<params>` is your params file. If you want to run multiple fits, put your params in the `staged` directory in the root of your project and run `python submit.py` on Supercloud.
