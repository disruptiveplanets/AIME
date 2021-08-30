# Converge data files

My goal is to get the fit to converge.

## August 4

I want to try setting a large cloud, so I have run 0. Then I implement Julien's advice of setting three parameters to the true values and varying one (runs 1.0-1.3)

## August 5

The fits look good! They're much closer to converging, but I need to run them again. I'm adding run 2 fits which have smaller spread.

## August 9

I just realized that I still have degeneracy in my model. So I'm removing it and rerunning some fits. The run 2 fits convinced me that big spread is better, or at least faster I think. So run 0 tests two fits for all parameters, with two different spread values. Run 1 tests fits only to one parameter, with the lower spread value.

## August 16

I'm going to rerun August 9, this time with reloading possible. Runs 0.* fit all three parameters at the same time, with different sigma and different degrees of initial separation. Fits 1.* fit one parameter at a time with different sigma. In this case, smaller sigma are better.

## August 18

It seems as though sigma controls the speed of convergence. I'm rerunning 0.1 and 0.2 which are the only fits that didn't converge and which I think have a chance of doing so. I'm also fitting runs 2.*, which fix one parameter and fit two.

## August 20

It seems hardest to fit to theta 2, possibly because it's the smallest compared to its sigma. So I'm shortening sigma and running a very long fit in preparation for vacation as 3.*. It fixes zero, one, and two parameters. All the fits converge except for three parameters and two parameters with theta 1 fixed; I'm running these for another 100,000 iterations.
