## Jan 5 & 6

I restarted this folder pursuant to the param search completing successfully. The first task is to adapt the code so that it works for l=3 and for sigma other than 0.01.
- I fixed bugs that prevented the tiered minimization working
- I need to find the scaling pattern between data sigma and result sigma by running the test code and finding how the Hessian eigenvalues scale with data sigma.
    - Sigma=e-2: 6.62540362e-04 2.43756485e-06 4.33391169e-07 6.75305954e-01 6.70141655e-01 5.58199672e-01 5.45043393e-01 4.60189682e-01 4.07491528e-01 3.37137117e-01
    - Sigma=e-3: 3.59251208e-05 2.43721152e-08 4.33398769e-09 7.07493544e-02 6.59762513e-02 4.92739649e-02 4.58450613e-02 3.90775900e-02 3.30948550e-02 3.22413912e-02
    - Sigma=e-4: 3.31099491e-06 2.36902195e-10 4.21283783e-11 1.07493565e-02 7.55600674e-03 7.28915976e-03 4.16987694e-03 3.57737289e-03 3.24460038e-03 2.75694140e-03
    - Sigma=e-5: 3.38815197e-07 2.98488230e-12 5.30800769e-13 1.70288302e-03 1.55018236e-03 7.85707165e-04 5.86544680e-04 5.17829492e-04 4.18368070e-04 3.69664785e-04
- Now let's try the true likelihood function:
    - Sigma=e-2: 1.22114038e-03 3.87518446e-05 1.47211470e-06 4.99365230e+00 1.44707464e+00 5.21729781e-01 3.14525469e-01 1.77536109e-01 5.08089245e-02 9.07788160e-03
    - Sigma=e-5: 9.68010976e-09 9.76839687e-11 1.80019698e-12 2.20580924e-05 8.05107737e-06 1.88920324e-06 3.86877265e-07 2.25057974e-07 6.22742431e-08 1.10804406e-08
- Things take forever to converge. Should I be worried?
    - No, I don't think so. This convergence takes forever with sigma=0.01, and I don't expect good convergence here because the parameters are so badly known.
- Then can I just run the fit, I think. Let's start with the smallest sigma.

I learned that in order to minimize for very low data sigma over l=3, I need to add in another minimization tier. After some experimenting, I found that a tier at drc ___ was sufficient. This was done via a supercomputer scan (Jan 6)
    - -2: Worked for sigma-5. However, nearly all sigmas didn't find their target. The second minimization (-2 to end) is hard for them, not so much the beginning. (Try -2 and something else? 4)
    - -3: Also worked
    - -4: Maybe?
    - -5: Way too far back for sigma-2, but doesn't matter because the minimization can be done in run round

The above information was used to set the maximum and minimum eigenvalues and the way they scale with sigma, as well as the minimization pattern for l=3

# Jan 13
I'm starting the run again with sigma theta to rho ratio of 1e-5. There may be some bugs still left over from the changes I made (corrected orientation update code and uncertainty model). I'm doing the asymmetric model, due to the irritating features of the symmetric model and its uncertainty with gamma zero