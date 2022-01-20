import density
from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

TRIALS_PER_TASK = 10

np.random.seed(3096518683**(rank+1)%3279872989)

klms = [
    1.2e6,
    0, 0, 0,
    2e3, 0, 0, 0, -1e4,
    ]

def is_inside(x, y, z):
    return x*x + y*y + z*z < 1

distro = density.Density(10, 1, klms, is_inside)
print("Populating")
distro.populate(TRIALS_PER_TASK, discard_negative=True)

if rank != 0:
    comm.Send([distro.densities, MPI.DOUBLE], dest=0, tag=77)
    comm.send(distro.trials, dest=0, tag=47)

if rank == 0:
    print("Waiting for other processes")
    for i in range(1, nprocs):
        print("{}/{} received".format(i, nprocs))
        new_densities = np.empty_like(distro.densities)
        comm.Recv([new_densities, MPI.DOUBLE], source=i, tag=77)
        newtrials = comm.recv(source=i, tag=47)
        distro.add_densities(new_densities, newtrials)

    print("Saving gif")
    for i, k in enumerate(klms):
        print(distro.get_klm(i))
    distro.save_densities("example.txt")
    print("Done")
