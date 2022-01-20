# Density

## Description
This directory explores a density distribution finding method I came up with. It develops some initial code to test whether the method works and tries it out on a few scenarios. It is only included in the `rotating` github branch, but is not specific to that branch.

The task is as follows. Suppose we have an accurate shape model of an asteroid, obtained from lightcurve analysis, and suppose we have Klms for that asteroid as well, obtained from a flyby. Then we want to impose a density distribution on the already-known shape to satisfy both the Klm and the shape constraints.

### Files
* **density.py**: the main module. Contains all the required code
* **call.py**: an example of how to use **density.py**.
* **example.gif**: the result of **call.py**.

### Use
The shape model of the asteroid is determined by a python function `is_inside`, which is passed as an argument to **density.py**. `is_inside` takes a point and returns `True` if the point is inside the asteroid and `False` otherwise.

Once this function is written and fit language is written like in **call.py**, the whole thing is called using

```mpirun -n {NUMBER OF CORES} python3 call.py```

The number of cores argument can be omitted to use the default value. The program can also be run without MPI, simply with

```python3 call.py```

## Method
1. Break the asteroid volume into *N* Voronoi cells, where *N* is the number of *Klm*s available.
    - You can do this by distributing *N* points uniformly randomly around a box in which the asteroid is contained, and redrawing all the points that are outside the asteroid until all points are inside. The construct cells around the points.
2. Set up a matrix *Klmn* which contains the integral of *Rlm* over the volume of Voronoi cell *n*, where each instance of *n* is a new column, and each *l, m* combination is a new row.
    - Note that this matrix is square, since *N* is equal to the number of *Klm*s.
3. Write the *Klm*s you want to achieve as a column vector, and the densities of each Voronoi cell as another column vector. Since the sum of the *Klmn* entry of every Voronoi cell equals the total *Klm* entry, we may write the matrix equation *Klm = (Klmn) (density)*.
4. Invert the matrix *Klmn* and solve for density.
5. Restart from step 1 to get another density map. Average these two density maps to get a third density map which still conforms to the shape model and still reproduces the *Klm*s, but no longer has the hard edges of the Voronoi cells.

As implemented in this code, some things such as a map of *Rlm* over all points in the asteroid are precomputed to decrease run time, but the overall method is not affected. This procedure is greatly parallelizable because every new instance of Voronoi cell choice can occur independently, with the average carried out only at the end. So it is implemented using MPI.
