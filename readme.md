# Asteroid tidal torque research project

## Description

This repo was created by Jack Dinsmore in June of 2021 to hold code, data, and progress reports for his research project on the effect of tidal torque on asteroids. The project was done with Professor Julien de Wit of MIT. All data referred to in the paper is either contained in this repo or can be generated via the scripts in this repo.

My goal is to be able to determine asteroid shape and density distribution based on observations of a close flyby to a planet. We will assume that light curve data has revealed the rotational velocity of the asteroid as a function of time.

The corresponding paper, _Constraining the Interiors of Asteroids Through Close Encounters_, was submitted to MNRAS in late July.

## File description
- **paper/** contains the paper associated with this work.
- **code/** contains both the C++ asteroid simulation and the fitting software.
- **data/** contains all the data, processing code, and figures computed on a supercomputer. Many of the large files are not on GitHub.
- **orientation/** contains a copy of the `code` directory, but specialized for a data set of asteroid orientation.
- **apophis/** contains some code copied from this directory specifically dedicated to the 2029 encounter of Apophis.
- **math/** contains Mathematica files used to do some of the analysis.
- **summaries/** contains progress reports and summaries of mathematical calculations.
