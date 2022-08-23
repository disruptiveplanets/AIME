## Aug 23

This readme is specifically for `delta-chisq.py`. Today, I invented the "smooth" distribution, which is uniform except it sweeps any center of mass offset into a linear trend in density. Moments need to be calculated for it. For move-1.5, the difference in density distribution between the low and high-density ends of the asteroid is 10%.

The smoothed moments are:
    K22 = -0.051987377852116975
    K20 = -0.20218802922577933
    IK33 = -0.00011727668026358638
    IK31 = -0.0005743325909623937

I used the code to calculate disagreement between these moments and the true moments.

This led to moment delta chisqs of 
    Move 1.5:   2409766477.3123336 (267751830.81248152)
    Sph 3:      11323332523.87726 (1258148058.2085845)
    Double:     35714773304.59333 (3968308144.9548144)
    (Confirmed that all true chi squares are good)

The density chisqs for lumpy were
    Move 1.5:   44085061042967.4 (6481657.902872656)
    Sph 3:      -103779770.21643382 (-15.258338177810144)
    Double:     -429852111.02974606 (-63.19949314648655)
    (Chi squares are not good for sph-3 and double)

And for finite element:
    Move 1.5:   -4971.586730161602 (-0.11798905283277011)
    Sph 3:      -35938.16997105251 (-0.8519384119820906)
    Double:     -44384.311664549234 (-1.0521598630890678)
    (Confirmed that all redchis with true are reasonable)

Negative means that uniform is more consistent. But remember, the uncertainties don't take into account correlations across the data set. Mention this in the paper, and the fact that the finite element model is so close to uniform.

The reason why the double chi squared is so bad is not because the model is inconsistent with the truth (it is) it's because there are so many minima that the MCMC can get stuck in the wrong one and not represent the true prob distro.

Next, make figures, put everything in the text, and justify a new threshold. I set it to be the lumpy uncertainty at the perigee value that gives fe uncertainty of 30%.