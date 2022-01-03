# Nov 14

I started this, as a repeat of probe-sigma-space 0, to get more data.

# Dec 12

After a couple days break to work on grad school and finish the GCE project, I fixed bugs in fit.py and will have to rerun. Used to be 24 cores per task, now 8, to see if speed increases. Runtime for 1 run is ~24 hrs local.

# Dec 20

I ran the previous test and it succeeded. I'm working on my compiling software now. I found a few bugs, but most importantly, I've changed the thin to exactly 10 in all cases so that I can automatically set the burn-in locally. Then I can get cleaner numbers.

# Dec 22

Having changed the minimization process so that the walkers more closely match the likelihood, I'm re-running the parameter search and oblateness test, with the hope that the minimization procedure will work. It works on the test and on some real parameters I tested locally.

Next step is to fix any more minimization problems, update the paper, and work on the intro.

# Jan 2

The minimization process failed, partially because of a segfault in the sim, and partially because the minimization technique sometimes failed. I added another tier of minimization, switched to Nelder-Mead, and fixed the segfault. Of those fits that succeeded, the contours are better but still not perfect.

But now that I look at the data properly, the fits that succeeded actually did quite well! The contours are OK!

I also changed the plotting software to save the original data and plot against that. Residuals are now more reliable.