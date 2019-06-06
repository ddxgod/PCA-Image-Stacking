from astropy.io import fits
import numpy as np
import scipy as sp
import os
from astropy.stats import sigma_clipped_stats

mgdirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/')
refdirs = os.listdir('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/')

print("MG")
for i in range(len(mgdirs)):
    filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + mgdirs[i]
    hdulist = fits.open(filename)
    scidata = hdulist[0].data.astype(float)
    mean, median, stddev = sigma_clipped_stats(scidata, sigma=3.0, iters=5)

    if mean < 1:
        print(mgdirs[i])

        
print("REF")
for i in range(len(refdirs)):
    filename = '/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + refdirs[i]
    hdulist = fits.open(filename)
    scidata = hdulist[0].data.astype(float)
    mean, median, stddev = sigma_clipped_stats(scidata, sigma=3.0, iters=5)

    if mean < 1:
        print(refdirs[i])
    
