## Jan 5

I restarted this folder pursuant to the param search completing successfully. The first task is to adapt the code so that it works for l=3 and for sigma other than 0.01.
- I fixed bugs that prevented the tiered minimization working
- I need to find the scaling pattern between data sigma and result sigma by running the test code and finding how the Hessian eigenvalues scale with data sigma. (I confirmed that the data sigmas do appear to be similar between the test and the true function)
- Things take forever to converge. Should I be worried?
    - No, I don't think so. This convergence takes forever with sigma=0.01, and I don't expect good convergence here because the parameters are so badly known.
- After finding the scaling, I need to confirm it with the true likelihood function
- Then can I just run the fit, I think. Let's start with the smallest sigma.