import numpy as np
import scipy as sp
from scipy import interpolate
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.stats import mad_std
from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u

import os
import linecache
import string
import math
import multiprocessing
from multiprocessing import Pool
import sys



# Calculates the distance btween two given points

def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)



#file = fits.open('SimulatedPSFSUBComb14-r3.fit')
file = fits.open('2175SUBComb30-r.fit')

data = file[0].data.astype(float)

signal = []
for i in range(len(data[0])):
	for j in range(len(data[0])):
		if distance(i, j, 25, 25) <= 1.4 / 0.396 * 1.5:
			signal.append(data[i, j])

#print(signal)
#print(22.5 - 2.5 * np.log10((np.sqrt(np.mean(np.array(signal - np.mean(signal))**2)))))
print(22.5 - 2.5 * np.log10(sum(signal) + 3 * np.sqrt(np.mean(signal)**2)))