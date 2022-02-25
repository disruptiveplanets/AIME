## Feb 24

Started this new density finder. Implemented the likelihood and harmonic model and tested them.

Next steps are
1. To run some sample fits to get correct uncertainties
2. To run the likelihood and harmonic models on all existing shape models
3. Design and implement the lump model

I also wrote some code in get-exotic-moments which computes the density moments of standard density distributions for use in making a fit with the same distributions. I computed them for an ellipsoid and a sphere to confirm that the program worked, then for the rest. Then I pushed all the sample asteroids to the supercomputer.

I'm a little worried about the blob indicator function. It might need to be numpy compatible, and it isn't.

(1+0j)
(0.004923771870408657+0.00397950275280673j)
(0.00548358106504681+0j)
(-0.004923771870408657+0.00397950275280673j)
(-0.05075584790467253+0.0009834059113698675j)
(0.0013408794786330584+0.0010914539656453523j)
(-0.19944944236136655+0j)
(-0.0013408794786330584+0.0010914539656453523j)
(-0.05075584790467253-0.0009834059113698675j)
(-4.779784624535005e-05+9.47940413744119e-05j)
(5.5145967143651265e-05+0.00026566118251419846j)
(-7.561156170951632e-05-5.882721025690665e-05j)
(-0.0004899253900947085+0j)
(7.561156170951632e-05-5.882721025690665e-05j)
(5.5145967143651265e-05-0.00026566118251419846j)
(4.779784624535005e-05+9.47940413744119e-05j)
(995.2700542582672+0j)