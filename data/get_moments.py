import os, sys
import numpy as np
from scipy.stats import moment
from scipy.special import factorial

if len(sys.argv) != 2:
    print("Please pass a samples file to find the moments of")

with open(sys.argv[1], 'rb') as f:
    array = np.load(f)
    array = array.reshape(-1, array.shape[-1])

std = np.std(array, axis=0)
print("Mean:", np.mean(array, axis=0))
print("Std:", np.std(array, axis=0))
for i in range(3, 11):
    m = moment(array, moment=i, axis=0) / std**i
    g = 2 / np.sqrt(2)**i * factorial(i-1) / factorial(i/2 - 1) if i % 2 == 0 else 0
    f = m - g / g if g != 0 else "--"
    print(f"Moment {i}: {m}.\tGaussian value: {g}.\tFraction difference: {f}")
