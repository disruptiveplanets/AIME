## April 29

Made this directory to contain estimates for how to get the sigma threshold. There are two approaches. One (the better one) is to determine how large sigma needs to typically be in order to have sigma_rho / rho = 1. Maybe =1 on average? In order to do this, I need a scan of uncertainty in sigma all over the parameter space up to K3m, and I only have K2m. That's in `scan-uncs`.

So I also am using another method, which computes how large the largest K3m is when you have a density distribution whose lowest point is rescaled to zero. That's in `max-klms`. It gives roughly the right answer.

Both of these use shape models consistent with the uniform ellipsoid with those k2ms.

In order to generate uncertainties for K3m, I'm running scan-space with K3m instead.