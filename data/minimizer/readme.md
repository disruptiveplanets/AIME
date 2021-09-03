# Converge data files

My goal is to get the fit to converge.

## Sept 1

I ran my first fit by minimizing first and it converged! The only task now is to make the initial minimization better. In the converge branch, I learned that l=3 parameters are very hard to fit to. Later, I checked this and confirmed that you can fit l=3 parameters, but only with a real honking big asteroid (for an Earth flyby). Order of km.

## Sept 2

Today I fixed the minimizer. The fit now runs exactly as it should, finding minima with the minimizer and then fitting on the lowest few points found.

Note: After I pull from the supercomputer, I need to un-comment-out the pull.py line that takes the .so from the fit_resolved directory. 
