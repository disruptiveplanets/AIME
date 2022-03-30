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