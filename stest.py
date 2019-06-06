# Test the simulation of one frame





from astropy.io import fits
from astropy.table import Column
from astropy.table import Table
import scipy as sp
from scipy.ndimage.filters import gaussian_filter
import scipy.optimize as opt
from scipy import interpolate
import numpy as np
np.set_printoptions(threshold=np.inf)
from numpy import *
import functools
import string
import decimal
import matplotlib.pyplot as plt
import fileinput
import math
import os
import photutils
from photutils import DAOStarFinder
from photutils import IRAFStarFinder
from photutils import subtract_psf
from photutils import centroid_2dg
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
import multiprocessing
from multiprocessing import Pool
import linecache
from sklearn.decomposition import NMF, PCA, IncrementalPCA
import astropy.units as u
import astropy.coordinates as coord
import sys
import itertools
sys.setrecursionlimit(15000)





def gaussian_2d(x, y, x0, y0, xsig, ysig):
	return np.exp(-0.5*(((x-x0) / xsig)**2 + ((y-y0) / ysig)**2))





# Calculates the distance btween two given points

def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)




# Method that calculates the total photon count within 3 sigma of the quasar centroid

def photonCount(xc, yc, sigma, data):
	count = 0
	for i in range(len(data)):
		for j in range(len(data)):
			if distance(i, j, xc, yc) <= sigma:
				count += data[i][j]
	return count




# Create the individual PSF frames

def psf_frame(comb_flux):
	FWHM = 1.4 / 0.396
	plate_scale = 1
	
	x_qso = 50
	y_qso = 50
	x = np.arange(0, 101, plate_scale)
	y = np.arange(0, 101, plate_scale)
	X, Y = np.meshgrid(x, y)
	QSO = gaussian_2d(X, Y, x_qso, y_qso, FWHM/2.355, FWHM/2.355)
	QSO *= comb_flux / sum(sum(QSO))
	print(sum(sum(QSO)))
	#plt.imshow(QSO + make_noise_image((101, 101), type='poisson', mean=0.001))
	

	'''
	add sky background and noise here
	'''
	sky = 0#170 * 0.0052
	final_image = QSO + sky
	final_image /= 0.0052
	#plt.imshow(final_image)
	#plt.show()
	for k, j in itertools.product(range(QSO.shape[0]), range(QSO.shape[1])):		
		arr = np.random.normal(final_image[k,j], final_image[k,j]**0.5, 1000)

		index_need = np.random.randint(arr.size, size=1)
		final_image[k, j] = arr[index_need]


	#plt.imshow(final_image)
	final_image *= 0.0052
	#print('before sky sub: ' + str(sum(sum(final_image))))

	final_image -= sky
	print('after sky sub: ' + str(sum(sum(final_image))))
	#print("%f, %f" % (np.min(final_image), np.max(final_image)))
	return final_image






psf = psf_frame(10**(-1.2/2.5))
plt.imshow(psf)
plt.show()
#print(10**(1/2.5))
print(sum(sum(psf)) / 10**(-1.2/2.5))