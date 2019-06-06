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




# Create the Absorber frame

def absorber_frame(flux, impact_parameter):
	FWHM = 0.6 / 0.396
	plate_scale = 1
	
	x_qso = 50
	y_qso = 50
	x_host = x_qso + impact_parameter
	y_host = 50
	x = np.arange(0, 101, plate_scale)
	y = np.arange(0, 101, plate_scale)
	X, Y = np.meshgrid(x, y)
	QSO = gaussian_2d(X, Y, x_qso, y_qso, FWHM/2.355, FWHM/2.355)
	QSO *= flux / sum(sum(QSO))
	host = gaussian_2d(X, Y, x_host, y_host, FWHM/2.355, FWHM/2.355)
	host /= sum(sum(host))
	host *= 10**(-1.2/2.5)

	'''
	add sky background and noise here
	'''
	sky = 170 * 0.0052
	final_image = QSO + host + sky
	final_image /= 0.0052
	#plt.imshow(final_image)
	#plt.show()
	for k, j in itertools.product(range(QSO.shape[0]), range(QSO.shape[1])):
		arr = np.random.normal(final_image[k,j], final_image[k,j]**0.5, 1000)
		index_need = np.random.randint(arr.size, size=1)
		final_image[k, j] = arr[index_need]


	final_image *= 0.0052
	#final_image -= sigma_clipped_stats(final_image, sigma=3.0, iters=5)[1]
	final_image -= sky
	return final_image
	



# Create the individual PSF frames

def psf_frame(flux, comb_flux):
	FWHM = 0.6 / 0.396
	plate_scale = 1
	
	x_qso = 50
	y_qso = 50
	x = np.arange(0, 101, plate_scale)
	y = np.arange(0, 101, plate_scale)
	X, Y = np.meshgrid(x, y)
	QSO = gaussian_2d(X, Y, x_qso, y_qso, FWHM/2.355, FWHM/2.355)
	QSO *= comb_flux / sum(sum(QSO))
	#plt.imshow(QSO + make_noise_image((101, 101), type='poisson', mean=0.001))
	

	'''
	add sky background and noise here
	'''
	sky = 170 * 0.0052
	final_image = QSO + sky
	final_image /= 0.0052
	#plt.imshow(final_image)
	#plt.show()
	for k, j in itertools.product(range(QSO.shape[0]), range(QSO.shape[1])):
		arr = np.random.normal(final_image[k,j], final_image[k,j]**0.5, 1000)
		index_need = np.random.randint(arr.size, size=1)
		final_image[k, j] = arr[index_need]


	final_image *= 0.0052
	#final_image -= sigma_clipped_stats(final_image, sigma=3.0, iters=5)[1]
	final_image -= sky
	return final_image





# Method that creates the simulaed final frames

def begin(index):
	i = int(index.split('-')[0])
	mgi = int(index.split('-')[1])
	color = index.split('-')[2]

	prop_data = TableDR12[mgi]

	qso_flux = prop_data['PSFFLUX'][2]
	largearr = []

	absorber_data = []
	for i in range(8):
		adi = absorber_frame(qso_flux, i)
		absorber_data.append(adi)



	# Calculate the residues with offsets from 1 to 8 pixels

	pixel = np.arange(0, 8)

	for offset in pixel:
		largearr = []
		adi = absorber_data[offset]
		comb_flux = photonCount(50, 50, 0.6 / 0.396, adi)

		for i in range(10):
			psf = psf_frame(qso_flux, comb_flux)
			largearr.append(np.reshape(psf, 10201))


		largearr = np.array(largearr)

		print(np.shape(largearr))



		# Set number of components in PCA, use Incremental PCA (IPCA) due to high efficiency and speed  
		numcomp = len(largearr)


		print(numcomp)
		mean_vector = []
	 
		#print(np.shape(largearr))

		try:
			for j in range(0, 10201):
				mean_vector.append(np.mean(largearr[:, j]))
		except:
			print("NO SOURCE FOUND")
			return

		largearr -= mean_vector
			
		ipca = IncrementalPCA(n_components=numcomp)
		ipca.fit(largearr)
		ipca_comp = ipca.components_
		#print(np.shape(ipca_comp))


		chunk_size = 50


		# Only use the components of the central portion of the quasar, since there may be overfitting due to strength of ipca
		
		new_comp = []
		for j in range(len(largearr)):
			temp = np.reshape(ipca_comp[j, :], (2 * chunk_size + 1, 2 * chunk_size + 1))
			new_comp.append(np.reshape(temp[47 : 54, 47 : 54], 49))

		new_comp = np.array(new_comp)
		new_comp = new_comp.T
		print(np.shape(new_comp))

		ipca_comp = ipca_comp.T

		
		#print(ipca_comp)
		take_final = numcomp    
		
		# Final fitting of the first n components, as determined by take_final, into the quasar to build a PSF fit
		#print(np.shape(ipca_comp))
		adi = np.reshape(adi, (2 * chunk_size + 1)**2)
		adi -= mean_vector
		adi = np.reshape(adi, (2 * chunk_size + 1, 2 * chunk_size + 1))
		coeff = np.dot(np.reshape(adi[47 : 54, 47 : 54], 49), new_comp)
		
		#coeff = np.dot(qsodata, ipca_comp)

		final_fit = np.dot(ipca_comp[:, 0:take_final], coeff[0:take_final])
		final_fit += mean_vector
		final_fit = np.reshape(final_fit, (2 * chunk_size + 1, 2 * chunk_size + 1))
		#final_fit /= len(largearr)

		adi = np.reshape(adi, (2 * chunk_size + 1)**2)
		adi += mean_vector
		adi = np.reshape(adi, (2 * chunk_size + 1, 2 * chunk_size + 1))



		mean, median, stddev = sigma_clipped_stats(final_fit, sigma=3.0, iters=5)
		final_fit -= median

		#final_fit *= adi[50][50] / final_fit[50][50]
		final_fit /= photonCount(50, 50, 0.6 / 0.396, final_fit)
		final_fit *= photonCount(50, 50, 0.6 / 0.396, adi)
		#mean, median, stddev = sigma_clipped_stats(final_fit)
		#final_fit -= median
		#final_fit *= np.max(adi) / np.max(final_fit)

		print("%f, %f" % (np.max(adi), np.max(final_fit)))




		# Final residue from subtraction of PSF from QSO

		residue = adi - final_fit
		#mean, median, stddev = sigma_clipped_stats(residue, sigma=3.0, iters=5)
		#residue -= median
		

		#fits.writeto('/data/marvels/billzhu/Sim.fit', np.array(absorber_data), clobber=True)

		try:
			fits.writeto('/data/marvels/billzhu/Simulated PSF Subtract/' + color + '06/' + str(index) + '-' + str(offset) + '_SUB.fit', residue, clobber = True)
			fits.writeto('/data/marvels/billzhu/Simulated PSF Cut/' + color + '06/' + str(index) + '-' + str(offset) + '_PSF.fit', final_fit, clobber = True)

			print('\n')

			print("DONE TO BOTTOM")
		except:
			print('HEADER IS CORRUPT')





# Code that opens up a maximum of 8 processes for concurrent execution
	
if __name__ == '__main__':
	#multiprocessing.set_start_method('spawn')
	TableDR12 = Table.read('DR12Q.fits')
	TableDUST = Table.read('final_catalog_full.fit')

	#gdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/g/')
	rdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/r/')
	#idirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/z/')
	#dirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/')

	#gdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/g/')
	#rdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/r/')
	#idirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/z/')
	#udirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/u/')


	#gdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/i/')
	#rdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/r/')
	#idirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/z/')
	
	#check_dirs = os.listdir('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/i/')
	rangelist = []

	#rangelist.append(gdirs[0].split('_')[0])
	#for d in gdirs:
	#	index = d.split('_')[0]
	#	#if str(index) + '_SUB.fit' not in check_dirs:
	#	rangelist.append(index)


	
	for d in rdirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	

	"""
	for d in idirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	
	for d in zdirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	"""
	
	print(len(rangelist))
	#try:
	pool = Pool(multiprocessing.cpu_count())
	pool.map(begin, rangelist)
	#except:
	#print("Error: Unable to process file")



