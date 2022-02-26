## Feb 24

Started this new density finder. Implemented the likelihood and harmonic model and tested them.

Next steps are
1. To run some sample fits to get correct uncertainties
2. To run the likelihood and harmonic models on all existing shape models
3. Design and implement the lump model

I also wrote some code in get-exotic-moments which computes the density moments of standard density distributions for use in making a fit with the same distributions. I computed them for an ellipsoid and a sphere to confirm that the program worked, then for the rest. Then I pushed all the sample asteroids to the supercomputer.

## Feb 26

I pulled the fits that were finished (asym, blob, in, out) and started experimenting with them. Wrote the interface to the supercomputer. But the code is a little broken; there are still some bugs to be worked out. The fit doesn't always converge.