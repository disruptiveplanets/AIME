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

I finished up the above tests and put them in the paper. I also thought of modelling a symmetric tetrahedral asteroid and perhaps a contact binary or an asteroid with a bumps, and implemented the tetrahedron and c.b.

Next I need to figure out algorithm runtime and consider asteroid with bumps? Maybe redundant with other two.