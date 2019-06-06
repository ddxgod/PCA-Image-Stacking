# Final simulation procedure, where each reference QSO is simulated through a Gaussian PSF, 10 PSF stars are created, and PCA PSF Subtraction is performed




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
				count += data[i, j]
	return count







def begin(index):
	#i = int(index)
	i = int(index.split('-')[0])
	#mgi = int(index.split('-')[1])
	color = index.split('-')[1]
	
	#try:
	print(index)
	
	try:
		qsocut = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + index + '_MG.fit')
		qsodata = qsocut[0].data.astype(float)

		psfcut = fits.open('/data/marvels/billzhu/MG II PSF Cut/0.37 - 0.55/' + color + '/' + index+ '_PSF.fit')
		psfdata = psfcut[0].data.astype(float)

		psfcount = photonCount(50, 50, 1, psfdata)
		psfdata /= psfcount
		qsocount = photonCount(50, 50, 1, qsodata)
		psfdata *= qsocount
		print("%f, %f" % (psfcount, qsocount))
		residue = qsodata - psfdata

		fits.writeto('/data/marvels/billzhu/MG II PSF Subtract 2/0.37 - 0.55/' + color + '/' + index + '_SUB.fit', residue, qsocut[0].header, overwrite = True)
	except:
		print('ERROR')



# Code that opens up a maximum of 8 processes for concurrent execution
	
if __name__ == '__main__':
	#multiprocessing.set_start_method('spawn')

	#gdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/g/')
	#rdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/r/')
	#idirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/z/')
	#dirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/')


	#gdirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/g/')
	#rdirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/r/')
	#idirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/z/')


	gdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/g/')
	rdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/r/')
	idirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/i/')
	zdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/z/')
	#udirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/u/')
	
	
	"""
	gdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/g/')
	rdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/r/')
	idirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/i/')
	zdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/z/')
	"""

	check_dirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/r/')
	rangelist = []
	#rangelist.append('99795-r')
	#rangelist.append('89-230830-g')

	#rangelist.append('100728-g')

	#rangelist.append('96034-g')

	#rangelist.append('95823-g')
	#rangelist.append(gdirs[0].split('_')[0])
	

	"""
	for d in gdirs:
		index = d.split('_')[0]
		if str(index) + '_SUB.fit' in check_dirs:
			rangelist.append(index)
	

	"""
	for d in rdirs:
		index = d.split('_')[0]
		if str(index) + '_SUB.fit' in check_dirs:
			rangelist.append(index)
	"""
	
	
	for d in idirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	

	
	for d in zdirs:
		index = d.split('_')[0]
		if str(index) + '_SUB.fit' not in check_dirs:
			rangelist.append(index)
	"""
	
	
	print(len(rangelist))
	#try:
	pool = Pool(multiprocessing.cpu_count())
	pool.map(begin, rangelist)
	#except:
	#print("Error: Unable to process file")


