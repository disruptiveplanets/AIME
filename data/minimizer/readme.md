# Converge data files

My goal is to get the fit to converge.

## Sept 1

I ran my first fit by minimizing first and it converged! The only task now is to make the initial minimization better. In the converge branch, I learned that l=3 parameters are very hard to fit to. Later, I checked this and confirmed that you can fit l=3 parameters, but only with a real honking big asteroid (for an Earth flyby). Order of km.

## Sept 2

Today I fixed the minimizer. The fit now runs exactly as it should, finding minima with the minimizer and then fitting on the lowest few points found.

Note: After I pull from the supercomputer, I need to un-comment-out the pull.py line that takes the .so from the fit_resolved directory.

## Sept 19

Today a run was completed wherein I scanned for l=3 parameters (run 4). However, none of them can resolve the high order parameters. If sigma is too low, the fit fails. But sigma = 0.01 works well. The lack of response to l=3 is probably because they have little effect on the dynamics (check with generate). I could force an effect by increasing the radius or lowering the perigee, but the radius is 100 m and the perigee is 5 earth radii: as low as it can be.

I'll run the fit again with radii of 1000 m, (run 5.0, 5.1, 5.2) and 10 km (5.3, 5.4, 5.5). Beyond that, the only thing to do is change J.

## Sept 23

The other fits failed to resolve a difference; some of them even failed to converge! On closer inspection, so many failed to converge that I fear a minimizer issue. I am therefore trying a tiered minimizer technique, where you minimize for small l, then fix those parameters and increment l. Run 6 is a rerun of 5 using this technique.

## Sept 27
However, I'm still running into problems with minima finding. Run 2 returns to l=2 order fits, but I don't have data for it because it's only useful for checking the minimizing parameters. I've learned that stopping at r = 2rp is worse than 5rp which is worse than 10 rp, but stopping at the end point of the simulation is bad. So, I'm running run-7, which is run-6 except minimizing up to rp=10. I will also run it for twice as long.

## Oct 2
I ran some fits again, and they are still running, but I have already determined that some are not converging. I'm moving over to l=2 to test my model again, with 100k iterations. It's run 8.

## Oct 5
Well, the l=2 run goes exactly as expected. I must now determine why the l=3 fit fails, despite running on the same code. I could try to half-ass this, or do it this weekend instead.

## Oct 10
It is now the weekend. I spent a long time on this code, building a test model which can be turned off and on with the TEST variable in display and fit.py. I used this model to make changes to the matrix calculation technique, the initial point picking, and extra log files. I got successful fits for the toy model. Now I want to see if I get good fits for the true model.

I'm using run 9 as my test, which contains both l=2 and l=3. I've also adjusted the fit and pull code so that the corner plotting happens on the cluster. That way, I don't get memory problems, and I don't have to store GB of data locally.

I had to introduce a fix that took the abs of the eigenvalues, because I wasn't getting pos def Hessians. Fix this in the future.
