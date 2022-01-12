# https://ssd-api.jpl.nasa.gov/doc/cad.html

import pandas as pd
import numpy as np
import json

DIAM_LIMIT = 100
DIST_LIMIT = 1# Lunar distances
YEAR_LIM = 2032

data = pd.read_csv("cneos_closeapproach_data.csv")
print(data.columns)
diameter_mask = np.full(len(data), False)
print(len(data))
for i in range(len(data)):
    if str(data["Diameter"][i]) == "nan":
        continue
    space_index = data["Diameter"][i].index(' ')
    unit = data["Diameter"][i][space_index+1:space_index+3]
    try:
        space_index = min(space_index, data["Diameter"][i].index('±'))
    except Exception:
        pass
    diam = float(data["Diameter"][i][:space_index])
    if unit == "km":
        diam *= 1000
    diameter_mask[i] = diam > DIAM_LIMIT

dist_mask = np.full(len(data), False)
for i in range(len(data)):
    line = data["CA Distance Minimum (LD | au)"][i]
    ld, au = line.split(" | ")
    dist_mask[i] = float(ld) < DIST_LIMIT

year_mask = np.full(len(data), False)
for i in range(len(data)):
    line = data["Close-Approach (CA) Date"][i]
    hyphen_index = line.index('-')
    year = line[:hyphen_index]
    year_mask[i] = float(year) <= YEAR_LIM


        
    
print(np.sum(diameter_mask))
print(np.sum(dist_mask))
print(np.sum(year_mask))
print(np.sum(diameter_mask & dist_mask & year_mask))
print(data[diameter_mask & dist_mask & year_mask])