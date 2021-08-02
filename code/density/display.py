from density import Density

distro = Density.load_densities("example.txt")
distro.save_gif("example2.gif", 2.0)
