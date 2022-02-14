# Density

## Jan 20

I basically wrote all the code in this folder today. Some of the plotting code was already there. These are the four models I have for making density distributions of asteroids. Each outputs a dat file of density distributions defined on a lattice, a png of a cross section, and a gif of a scan through the asteroids. Parameters can be tweaked in setup.py.

Next, I need to run some checks (How to make agreement better for harmonic?) (Is my rlm / ylm correct?). Then I need to think about some possible new models that could be added (polygonal?) Then I need to think of some good test cases, make data, and put this in the paper.

## Jan 21

I finished my checks and started making paper figures. I also made the slice plot. The paper figures are 
    * Symmetric and asymmetric models with correct ellipsoid forms
    * Symmetric and asymmetric models with spheres
Need to review the surface model I think, and the layer figure. Then put these in the paper.

## Jan 26 (after wisdom teeth hiatus)

I finished up the above tests and put them in the paper. I also thought of modelling a symmetric tetrahedral asteroid and perhaps a contact binary or an asteroid with a bumps, and implemented the tetrahedron and c.b. I also fixed the issue with the surface model.

## Jan 27

I refined the parameters to get better estimates and pushed everything to the supercomputer. It appears that this will take a while.

## Jan 28

Logged in today and everything was finished. I ran all the displays.

All the algorithms executed accurately except for surface, which has some slight issues. The ensemble method repeated seeds across threads, so I have fixed that bug and fixed it.

## Jan 29

I pulled the ensembles and ran the displays. Everything looks good

## Feb 13

I added some l=3 components and didn't get sensible results (big density blowups). At first I thought this was a bug, but in fact it's because the high order klms are just inconsistent with the shape. I wrote likelihood-pt-2 to debug this. Now I'm using a sphere with no 2 components to test the 3 components, and you still get some issues.

## Feb 14

I made the test directory, which creates a bunch of random density distributions and gets their Klms. It's difficult to get Klm much over 1e-2. So I generated random Klms with abs < 0.01 and reran in the supercomputer.