from ctypes import *

dll = CDLL('./main.so')
dll.foo.restype = py_object
a = dll.foo('foo')
print(a)
