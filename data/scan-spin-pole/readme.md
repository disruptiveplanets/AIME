## Jan 28

Created today. Wrote scatter code and untested compile code.

## Jan 29

Pushed to supercomputer while gap 2.0 and peri are still running

## Jan 31

Pulled data today with two processes still running (16, 21). Code works. Waiting for two processes to finish.

## Feb 1 

Pulled last two data points, made display prettier, and wrote text. However, I hit a snag when I realized that I could not explain the white spots at \pm y. So I made a covariance invariant formulation which still has white spots, although not as bright. I learned that you always have zero torque when the phi coordinate is equal to alpha and beta is pi/2. Maybe torque engages at the semi latus rectum point? Can't figure out why spin pointing along y would give low precision, when it gives torque along z. However, spin along x gives zero torque and is observable.

## Feb 4

In generate, I did some bug testing. I believe that the bald spot at +- Y is a balancing effect over the orbit, choosing an alpha that minimizes total torque. Looking at one point in the orbit is not valid. Wrote up the paper and I'm counting this issue as finally finished. I also made the symmetric plot for the display of dependence on initial conditions.

## March 31

I finally resolved the problem of why the data is more uncertain at Y. It's because here tau propto z, whereas at x it is zero and third order and non-perigee effects dominate. Propto z is bad for certainty because some klm are poorly constrained and because it induces no tumbling. The section is now done.