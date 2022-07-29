# Uncertainty threshold

There are two methods here to determine the maximum Klm uncertainty permissible to produce useable density distros. They are 

1. `max-klms`: a file to compute the maximum K_lm value attainable for uniform asteroids which have non-negative density. You would then use this as the upper uncertainty limit.

2. `scan-uncs` and `compile`: a file to compute the uncertainty in density moments which will produce density models with uncertainty equal to the density.

## April 29

Made this directory to contain estimates for how to get the sigma threshold. There are two approaches. One (the better one) is to determine how large sigma needs to typically be in order to have sigma_rho / rho = 1. Maybe =1 on average? In order to do this, I need a scan of uncertainty in sigma all over the parameter space up to K3m, and I only have K2m. That's in `scan-uncs`.

So I also am using another method, which computes how large the largest K3m is when you have a density distribution whose lowest point is rescaled to zero. That's in `max-klms`. It gives roughly the right answer.

Both of these use shape models consistent with the uniform ellipsoid with those k2ms.

In order to generate uncertainties for K3m, I'm running scan-space with K3m instead.