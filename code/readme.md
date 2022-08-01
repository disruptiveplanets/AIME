# Code

* The `sim` directory contains code to simulate the asteroid's flight. It is not executable
* The `generate` directory generates data for a set parameterization of an asteroid, which can be changed by the user. This is mostly for debugging purposes
* The `fit_resolved` generates sample data for a specific asteroid (described in a file in `../staged`) and extracts moments from it
* The `density` directory extracts a density distribution from these moments
* The `thresholds` directory iterates over data that has been collected, generates density distributions for each, and determines at what point those distributions become too uncertain.
* The `other-figs` directory contains other figures that might be useful.

# Note

Please note that, due to changes in the names and definitions of variables as the project was developed, the variables in the code no longer reflect those defined in the paper. For example, the definitions of klm and jlm used in the `sim` directory do not use $a_A$ (defined in the paper); they use $a_m$ which is defined such that $I_A = \mu_A a_m^2$. In equation 1 for $K_{lm}$, $a_m$ is used in place of $a_A$. In some places, $a_A$ is referred to as `a_surface` and $a_m$ is referred to as `a_bulk`.

For uniform asteroids, $a_A = a_m$ and there is no inconsistency. But for non uniform asteroids, the density moments obtained by a fit with this definition of $a_m$ must be multiplied by factors of ($a_m / a_A$) to obtain density moments consistent with the definitions in the paper. This is done in the code so that **the final data products accurately reflect the definitions used in the paper.** But intermediate data products may be inaccurate.