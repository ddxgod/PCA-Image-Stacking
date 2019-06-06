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
import itertools


# Create the individual PSF frames

def psf_frame():	

	'''
	add sky background and noise here
	'''
	final_image = np.ones((101, 101))

	sky = 170 * 0.0052
	final_image *= sky
	final_image /= 0.0052
	#plt.imshow(final_image)
	#plt.show()
	for k, j in itertools.product(range(final_image.shape[0]), range(final_image.shape[1])):		
		arr = np.random.normal(final_image[k,j], final_image[k,j]**0.5, 1000)

		index_need = np.random.randint(arr.size, size=1)
		final_image[k, j] = arr[index_need]


	#plt.imshow(final_image)
	final_image *= 0.0052
	#print('before sky sub: ' + str(sum(sum(final_image))))

	final_image -= sky
	#final_image -= sum(sum(final_image)) / 10201
	print('after sky sub: ' + str(sum(sum(final_image))))
	#print("%f, %f" % (np.min(final_image), np.max(final_image)))
	return final_image / math.sqrt(1582)



file = fits.open('SimulatedPSFSUBComb14-r.fit')
scidata = file[0].data.astype(float)

scidata += psf_frame()
#for i in range(len(scidata)):
#	scidata[i] += psf_frame()

fits.writeto('SimulatedPSFSUBComb14-r2.fit', scidata, clobber=True)