# Asteroid tidal torque research project

## Description

This repo was created by Jack Dinsmore in June of 2021 to hold code, data, and progress reports for his research project on the effect of tidal torque on asteroids. The project was done with Professor Julien de Wit of MIT.

My goal is to be able to determine asteroid shape and density distribution based on observations of a close flyby to a planet. We will consider both the resolved and unresolved cases. In the resolved case, full knowledge of the asteroid's rotation data is assumed as a function of time. In the unresolved case, light curve data is assumed.

In both the resolved and unresolved cases, we will build a simulation to generate observables for a given parameterization of the asteroid. We will use the simulation to generate sample data, then fit to that data to test the strength of the fit. We will do this for multiple asteroid shapes for completeness.

## File description
- **summaries/** contains progress reports and summaries of mathematical calculations.
- **math/** contains Mathematica files used to do some of the analysis.

## Project timeline
1. Create a robust and general model for asteroid shape and density profile (Done)
2. Write code to simulate an asteroid flyby for certain parameters

### Resolved case

3. Generate a sample dataset for asteroid rotation data with the simulation in step 2
4. Use an MCMC to fit a shape and density distribution to the data generated in step 3

### Unresolved case
5. Extend step 2 to generate light curves for the asteroid
6. Generate a sample dataset for asteroid rotation data with the simulation in step 5
7. Use an MCMC to fit a shape and density distribution to the data generated in step 6

### Actual data (Optional step)
8. Apply the fitting software to an actual flyby.

The intent is to start writing the paper for this project by the end of the summer of 2021, the summer between Jack's junior and senior years.