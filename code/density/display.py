from density import Density

distro = Density.load_densities("example.txt")
distro.save_gif("example.gif", 2.0)
