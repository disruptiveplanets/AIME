# Code

* The `sim` directory contains code to simulate the asteroid's flight. It is not executable
* The `generate` directory generates data for a set parameterization of an asteroid, which can be changed by the user. This is mostly for debugging purposes
* The `fit_resolved` generates sample data for a specific asteroid (described in a file in `../staged`) and extracts moments from it
* The `density` directory extracts a density distribution from these moments
* The `thresholds` directory iterates over data that has been collected, generates density distributions for each, and determines at what point those distributions become too uncertain.
* The `other-figs` directory contains other figures that might be useful.
