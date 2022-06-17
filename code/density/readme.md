## Feb 24

Started this new density finder. Implemented the likelihood and harmonic model and tested them.

Next steps are
1. To run some sample fits to get correct uncertainties
2. To run the likelihood and harmonic models on all existing shape models
3. Design and implement the lump model

I also wrote some code in get-exotic-moments which computes the density moments of standard density distributions for use in making a fit with the same distributions. I computed them for an ellipsoid and a sphere to confirm that the program worked, then for the rest. Then I pushed all the sample asteroids to the supercomputer.

## Feb 26

I pulled the fits that were finished (asym, blob, in, out) and started experimenting with them. Wrote the interface to the supercomputer. But the code is a little broken; there are still some bugs to be worked out. The fit doesn't always converge.

## Feb 27

No more of the fits were finished, but the density extractions were done. They had a couple bugs left in them which I fixed today. The density distribution has been thoroughly checked and seems very correct; the uncertainty also seems correct and I can't find any problems

The data for the asymmetric ellipsoid was made for the wrong k22; I need to redo it. I pushed all density extractions to supercloud and that k22 calc.

## Feb 28

I pulled the results for mapping the asteroid today. They look only slightly screwy (not uniform enough density), but that's probably expected. These models are designed to be hard to fit to. I'm looking forward to the asymmetric fit. The story there may be that the method is easily distracted by low-certainty K3m components. A very well-supported result to put in a paper, but a disappointing one.

I also have been working on the lumpy model. I've realized one symmetric rock is too degenerate. The mass and lever arm are degenerate, so you can choose one arbitrarily, and then the radius is determined by the radius of the asteroid. There are no other free parameters. I can add more rocks, but I feel that this is a pretty bad problem. I'd like to talk to Julien about it.

In the mean time, maybe it's best to take a break on the lumpy model. I pushed the asymmetric model to the supercomputer.

## March 5

I spent yesterday and today pushing all the fits to the plotters and generating figures. That worked today, but I sent the figs back to be re-plotted to speed up the gifs and remove a typo.

## March 7

Today I updated the lumpy model. It works. Here are the untested components though:
- Once an ellipsoid lump is solved for, is it inserted into the asteroid with correct orientation to get the density map?
- Implementation bugs for N > 1 (I know these bugs exist)

I will not test them until I have a use for them; I'd like to talk briefly with Julien about exactly which examples he thinks are interesting, then write and test just enough code to implement those examples.

Knowing that everything works, I re-pushed the blob calculation to the supercomputer, both because I moved the blob and because I wanted to check that the K3m components, which are quite important, were precise. I also ran lumpy-1-1 for all models except blob.

In case I end up testing the rotated blob case, I generated an example called rot-blob. But to find moments for it, we'd need to rotate it into the principal axis frame, which sounds deeply annoying. So I'm not going to yet.

I also learned of a mistake in my in and out density distributions, so I recalculated those and pushed the fits to the supercomputer again. Now blob, in, and out are running.

## March 9

Pulled data from supercomputer. Learned that the ellipsoid radii are imaginary for negative mass (like out example) but I'll have to fix it later because I'm in Caltech and sleepy.

## March 12

I fixed the problem with asteroid dimensions. Rerunning out.

## March 13

I pulled the figure and now neither in nor out has a lump, but I think that's OK because the COM is not affected by the presense of a lump at the origin. If it were an axial lump (d=3) then the ratios a/c, b/c would be constrained but still not c because the COM is not affected?

## March 28

Today I learned that the reason why the radius of the lump for the in and out models is unconstrained is because the shell has Klm proportional to the true Klm, by virtue of the density being foliated. The shell's mass therefore is set via this proportionality constant and the radius of the lump is set to zero to avoid error in the am estimate.

Therefore, I'm changing the in and out definition from a foliated distribution to a spherically symmetric one. Hopefully this will yield better results.

Submitted the simulation for the new in and out moments

## March 29

The new in and out fits finished, so now I'm executing all three models on them. I also pulled the data locally and executed with a very sparse grid; you do get the correct lump. So it should work.

I pulled in and out data today. Looks correct.

I realized that the blob data may be for an old distribution, since the text file is out of date, so I'm rerunning.

I confirmed by directly plugging in the Klms that the density distribution is indeed correct for lumpy blob if there are no high order K3ms, and happens to be true for exactly accurate K3ms too. Simulating the run on the sc then running the solver for blob.

I also learned that the uncertainty equation was slightly wrong, so I'm rerunning all the figures so that it can be redone.

Uncertainty can be very low sometimes, when a certain point depends only on Klm components which are not uncertain. The equatorial plane seems to be an example, but I don't understand why; they depend on both inertial axes, right?

Maybe this doesn't matter. I'm just going to run the figs tomorrow.


## March 30

I fixed a bug with square rooting a complex number and redrew the shapes in the Sc.

## April 19

Fixed the uncertainties bug that was leading to negative uncertainty.

## April 22

Made new figs. Everything, including the uncertainty, should be correct.

## June 2

Found moments for a core asteroid. They are:
    Elliptical core
        K22 = -0.05200645295935643 (old)
        K20 = -0.20219868297920368 (old)
        aA = 891.777
    Spherical core
        K22 = -0.03973412605878553
        K20 = -0.1754365400999416
        aA = 851.5964018739621
Pushed them to supercomputer to run extraction.

The core densities are three times the shell densities (ceres; Julie C. Castillo-Rogez and Thomas B. McCord, 2010) and the radii were chosen such that the core and shell have equal mass.

## June 3

Sent a lot of things to Julien. Changed the prior density bounds so that only asym-ell, den-core-ell, and den-core-sph are valid.

I also realized that the core of Spherical core might extend beyond the asteroid. So I manually set its radius (500 m) such that it didn't. New moments are 
    Spherical core
        K22 = -0.05040959164683727
        K20 = -0.19599016007639866
        aA = 922.9234884822591

I set the new spherical shell running in the supercomputer.

## June 6

Pulled the spherical core run from the supercomputer and ran lots of fits to it. Mostly still running on the supercomputer, but their priors are 0.5-9 and they run full number of elements.

## June 7

The average appears to be going well. It had a core density of 3x the shell and prior 0.5-9. I'm a little confused as to why the core only appears to have density 2x the shell but I think it could be the fault of the prior or not taking enough averages. Rather than continue down this path, I moved to a core of 1.5x the shell and 0.25-3 prior. Currently extracting moments. I'm also upping the DIVISION from 99 to 49. The new moments are
    Spherical core
        K22 = -0.0515978274784774
        K20 = -0.20060996601498285
        aA = 978.4541044108308

I also ran many density extractions on the new data.

## June 8

I finished averaging the new data, but I don't like it. The K3m moments are too high and they move the core off-center. I think it's because so few walkers converge; the image shows that the true values are outside the contours, and only 1 walker converged. I would like to achieve a better fit to these moments, so I reran the extraction. We will see what happens.

## June 13

I'm also doing the average for the 3x density moments. These are stored in avg-core-1.5 and avg-core-3 for cores of 1.5 or 3x the surrounding density. Furthermore, I'm doing one for a moved core of 1.5x surrounding density, 500 m in radius (same core) but moved 300 m. 

    Spherical core x3
        K22 = -0.05040959164683728
        K20 = -0.19599016007639866
        aA = 922.9234884822591
    Spherical core x1.5
        K22 = -0.05159782747847741
        K20 = -0.20060996601498282
        aA = 978.4541044108308
    Moved core x3
        K22 = -0.05203775196773803
        K20 = -0.19716764233198797
        IK33 = -0.0015577797254511872
        IK31 = -0.007833505702120666
        aA = 933.1648422811957
        surface_aA = 1002.0081758422925
    Moved core x1.5
        K22 = -0.05195945615840046
        K20 = -0.20085472167788565
        IK33 = -0.0003845813086399978
        IK31 = -0.0019368958923091067
        aA = 980.8811439828254
        surface_aA = 1000.1281468600504

I had to adjust my get_exotic_moments file a bit, and main, in order to jive with the new moments. I'm going to do recompute all moments just to make things more certain. I also recomputed the blob moments to ensure that they match the old blob moments.

## June 14

The moments code was running out of memory, unbeknownst to me. I got it working, computed the moments, confirmed that the blob K2 and K3 moments are reproduced, then ran the moment extraction for all four core types.

## June 16

I rewrote mcmc so that it's flexible and can be executed for multiple models. I also pulled the core extractions but they had a bug.

## June 17

I fixed a couple flaws in the code that were failing for non-9 degrees of freedom. I also tried to figure out what was causing the move-3 and move-1.5 core extractions to fail. Their moments were different from the data moments.