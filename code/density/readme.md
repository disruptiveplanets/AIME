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

# Feb 15

I drew the computation from the supercomputer today and it worked. Wrote this up in the summary. Then computed the higher order moments of each surface in the test directory: compute-moments. The results were:

Ellipsoid:
    1
    9.872208985875585e-06   -5.988990207635196e-06  1.4768533426670653e-06
    0.0520078479466095      5.012975622254911e-06   9.229553849303065e-07   5.1694428460585835e-08  -0.20219088782517713
    2.6165713950624314e-07  -9.245164121042991e-08  -1.4112544753842155e-06 -7.558819324529767e-08  -1.9172846960120366e-06 1.1389056361745457e-06  1.547794652201397e-06

Sphere: 
    1
    5.31317868864823e-06    -5.31317868864823e-06   1.062635737729647e-05
    -3.2290956144954023e-24 1.253894317039141e-06   -2.5077886340783027e-06 2.5077886340782692e-06  -3.2290956144954025e-23
    4.858140521680214e-07   4.858140566039582e-07   1.0611029770443995e-21  2.313772663612163e-06   -1.4574421698215781e-06 1.4574421698215781e-06  3.886512452857546e-06

Tetrahedron:
    1
    -2.13001803656672e-06   1.2738036364613414e-06  8.713672218542483e-05
    2.13001803656672e-06    1.2738036364613414e-06  -8.537198871839977e-06  2.4798978685652396e-07  -1.8206721115635396e-06
    2.6853926316436226e-07  -1.332727705764822e-07  -0.035859509908525926   -1.846508022712845e-09  -5.188517778675532e-07  -3.678734768900552e-07  -6.347944163140826e-07

Dumbbell: 
    1
    2.2046275124745951e-07  -1.5066996314417975e-06 3.013399262883594e-06
    0.041119225728060795    -1.3969817668597828e-06 2.7939635337195765e-06  -3.3476919773238653e-06 -0.08223845145612159
    8.41226036772929e-08    1.1146775335242793e-07  -6.504760357466113e-08  -7.683153619586698e-07  -2.523678032297567e-07  2.043080528937105e-07   7.61646819632193e-07
