# Nov 14

I started this, as a repeat of probe-sigma-space 0, to get more data.

# Dec 12

After a couple days break to work on grad school and finish the GCE project, I fixed bugs in fit.py and will have to rerun. Used to be 24 cores per task, now 8, to see if speed increases. Runtime for 1 run is ~24 hrs local.

# Dec 20

I ran the previous test and it succeeded. I'm working on my compiling software now. I found a few bugs, but most importantly, I've changed the thin to exactly 10 in all cases so that I can automatically set the burn-in locally. Then I can get cleaner numbers.

# Dec 22

Having changed the minimization process so that the walkers more closely match the likelihood, I'm re-running the parameter search and oblateness test, with the hope that the minimization procedure will work.