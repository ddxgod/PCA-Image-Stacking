from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import scipy as sp
import os
import linecache
import math
from scipy import interpolate
from scipy.stats import powerlaw
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from decimal import Decimal
import powerlaw

#pseudo code for SB profile and error calculation here:
def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)



def photoncount(scidata, radius1, radius2):
	flux = 0
	length = 0
	#print(np.shape(scidata))

	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if distance(i, j, 50, 50) <= radius1 and distance(i, j, 50, 50) >= radius2:
				flux += scidata[i][j] / 2000 / 10**8
				length += 1

	return flux / length


#load the data here:
hdulist = fits.open('TotComb50-g.fit')
scidata = hdulist[0].data.astype(float)

scale = cosmo.kpc_proper_per_arcmin(0.48) * u.arcmin / u.kiloparsec * 0.396 / 60 #kpc/pixel


#get the size of image:
width = scidata.shape[0]
max_boundary = width/scale*0.5


boundaries = [3, 10, 13, 16, 19, 27, 37, 51, 67, 100, 140] #kpc


#calcluate the sky error here:
sky_count = photoncount(scidata, boundaries[-1] / scale,  max_boundary/ scale)
sky_flux = sky_count * 2000 * 10**8

SBarray = []
error_array = []
outter = []
for j in range(len(boundaries) - 1):

	f_count = photoncount(scidata, boundaries[j + 1] / scale, boundaries[j] / scale)
	print("%f, %f" % (boundaries[j + 1], boundaries[j]))
	print(f_count)

	#outter.append(6 * (25 / 3)**((1.0 / 8)*j))
	outter.append(boundaries[j + 1] / scale * 0.396) # the unit is arcsec
	f_flux = f_count * 2000 * 10**8 #covert the count

	error = np.sqrt(f_flux)
	error / np.sqrt(2000 * 10**8)
	f_flux / 2000 / 10**8

	#flux per pxiel**2 to surface_brightness per kpc**2
	mag = -2.5 * np.log10(f_count/scale**2)

	# start to calculate the error by assuming is the gaussain error.
	mean_f = f_flux
	sigma_f = error
	f_array = np.random.normal(mean_f, sigma_f, 1000)
	mag_array = -2.5 * np.log10(f_array/scale**2)
	error_mag = astropy.sigma_clipped_stats(mag_array)[3]
	SBarray.append(mag)
	error_array.append(error_mag)
