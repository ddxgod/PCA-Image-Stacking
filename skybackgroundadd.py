# Calculate the full image sky background to be added back to the original image Regular bilinear interpolation is sufficiently accurate




import numpy as np
import scipy as sp
from scipy import interpolate
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.stats import mad_std
from astropy.table import Table
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
import photutils
from photutils import DAOStarFinder
from photutils.centroids import centroid_2dg
import os
import linecache
import string
import math
import multiprocessing
from multiprocessing import Pool
import sys
#import ezgal
sys.setrecursionlimit(15000000)



f = fits.open('95-887-r.fit')
field_data = f[0].data
field_data /= f[0].header['NMGY']
sky_binn = f[2].data['ALLSKY']
spline = interpolate.interp2d(np.arange(256), np.arange(192), sky_binn)
print(f[2].data['XINTERP'])
sky_data = spline(f[2].data['XINTERP'].ravel(), f[2].data['YINTERP'].ravel())
fits.writeto('sky_test.fit', sky_data, header=f[0].header, clobber=True)