import numpy as np
from core import Method
from scipy.linalg import pinv
import scipy.sparse

class Lumpy(Method):
    def __init__(self, asteroid):
        print("Lumpy model")
        super().__init__(asteroid)


    def get_a(self):
        raise NotImplementedError


    def get_b(self, x,y,z):
        raise NotImplementedError
