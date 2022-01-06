## Jan 5

I restarted this folder pursuant to the param search completing successfully. The first task is to adapt the code so that it works for l=3 and for sigma other than 0.01.
- I fixed bugs that prevented the tiered minimization working
- I need to find the scaling pattern between data sigma and result sigma by running the test code and finding how the Hessian eigenvalues scale with data sigma.
    - Sigma=e-2: 6.62540362e-04 2.43756485e-06 4.33391169e-07 6.75305954e-01 6.70141655e-01 5.58199672e-01 5.45043393e-01 4.60189682e-01 4.07491528e-01 3.37137117e-01
    - Sigma=e-3: 3.59251208e-05 2.43721152e-08 4.33398769e-09 7.07493544e-02 6.59762513e-02 4.92739649e-02 4.58450613e-02 3.90775900e-02 3.30948550e-02 3.22413912e-02
    - Sigma=e-4: 3.31099491e-06 2.36902195e-10 4.21283783e-11 1.07493565e-02 7.55600674e-03 7.28915976e-03 4.16987694e-03 3.57737289e-03 3.24460038e-03 2.75694140e-03
    - Sigma=e-5: 3.38815197e-07 2.98488230e-12 5.30800769e-13 1.70288302e-03 1.55018236e-03 7.85707165e-04 5.86544680e-04 5.17829492e-04 4.18368070e-04 3.69664785e-04
- Now let's try the true likelihood function:
    - Sigma=e-2: 
    - Sigma=e-3: 
    - Sigma=e-4: 
    - Sigma=e-5: 
- Things take forever to converge. Should I be worried?
    - No, I don't think so. This convergence takes forever with sigma=0.01, and I don't expect good convergence here because the parameters are so badly known.
- Then can I just run the fit, I think. Let's start with the smallest sigma.