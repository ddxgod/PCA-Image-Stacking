# First iteration of experiment simulation. using real reference quasar data and host galaxy imaging data



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
sys.setrecursionlimit(15000)



# Checks if the given coordinates are within the bounds of the image

def inbounds(x, y):
	if x > 0 and x < 2048 and y > 0 and y < 1489:
		return True
	return False



# Calculates the distance btween two given points

def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)



# Method that checks if a particular residue is noise or not by drawing a ring of similar distance as the "bright" point around the centroid

def checkNoise(x, y, qX, qY, data):
	halo = []
	for i in range(42):
		for j in range(42):
			if abs(distance(i, j, qX, qY) - distance(x, y, qX, qY)) <= 2 and distance(i, j, x, y) >= 2:
				halo.append(data[j, i])
	mean, median, std = sigma_clipped_stats(halo, sigma=3.0, iters=5)

	if data[y, x] > mean + 3 * std:
		return True
	return False
				

	
# Method that calculates the total photon count within 3 sigma of the quasar centroid

def photonCount(xc, yc, sigma, data):
	count = 0
	for i in range(len(data)):
		for j in range(len(data)):
			if distance(i, j, xc, yc) <= sigma:
				count += data[i][j]
	return count



# Method that checks if a source > 1 sigma is outside of a certain radius, and if so, masks it by putting mean value

def checkOutter(data, mean, std, visited):
	count = 0
	for i in range(len(data)):
		for j in range(len(data[0])):
			if data[i][j] > mean + 3 * std and visited[i][j] == False:
				data[i][j] = mean
	return data



# Checks if a sources > 2 sigma is inside of a certain radius, and if so, masks the source by calling floodfill()

def checkInner(data, sources, qX, qY, mean, stddev, pointer):
	count = 0
	for i in range(len(sources)):


		#yes = False
		#if abs(sources['COLC'][i][pointer] - qX - 3) < 3 and abs(sources['ROWC'][i][pointer] - qY - 46) < 3:
		#    print(True)
		#    yes = True

		if distance(sources['COLC'][i][pointer] - qX + 50, sources['ROWC'][i][pointer] - qY + 50, 50, 50) > 8:
			visited = np.zeros((len(data), len(data)), dtype=bool)


			# Radial masking with the isophotal radii given is inefficient

			#if sources['iso_a'][i][pointer] > 0:
			# 	data = radiusmask(data, int(sources['COLC'][i][pointer] - qX + 50), int(sources['ROWC'][i][pointer] - qY + 50), sources['iso_a'][i][pointer], mean, visited)

			
			#print("%f, %f" % (sources['COLC'][i][pointer] - qX + 50, sources['ROWC'][i][pointer] - qY + 50))
			#else:
			data = floodfill(data, int(sources['COLC'][i][pointer] - qX + 50), int(sources['ROWC'][i][pointer] - qY + 50), mean, mean + stddev, visited)

			#print(True)
			#print("%f, %f" % (sources['COLC'][i][pointer] - qX + 50, sources['ROWC'][i][pointer] - qY + 50))


			# If no isophotal radii are given, then use radius = 6 pixels as an initial estimate

			"""
			for j in range(-6, 6):
				for k in range(-6, 6):
					if distance(sources['COLC'][i][pointer], sources['ROWC'][i][pointer], sources['COLC'][i][pointer] + j, sources['ROWC'][i][pointer] + k) < 6 and int(sources['ROWC'][i][pointer] - qY + 50 + k) >= 0 and int(sources['ROWC'][i][pointer] - qY + 50 + k) < 101 and int(sources['COLC'][i][pointer] - qX + 50 + j) >= 0 and int(sources['COLC'][i][pointer] - qX + 50 + j) < 101 and data[int(sources['ROWC'][i][pointer] - qY + 50 + k)][int(sources['COLC'][i][pointer] - qX + 50 + j)] > mean + 3 * stddev:
						data[int(sources['ROWC'][i][pointer] - qY + 50 + k)][int(sources['COLC'][i][pointer] - qX + 50 + j)] = mean
						#print(data[int(sources['COLC'][i][pointer] - qX + 50 + j)][int(sources['ROWC'][i][pointer] - qY + 50 + k)])
			"""

	return data




# If the isophotal major axis is given, then use radiusmask() for accurate masking

def radiusmask(data, xc, yc, isophotal_radius, mean, visited):
	if xc < 0 or xc >= len(data) or yc < 0 or yc >= len(data):
		return data

	for r in range(len(data)):
		for c in range(len(data)):
			if distance(r, c, xc, yc) < 1.2 * isophotal_radius / 2 + 3:
				data[c][r] = mean

	return data




# Main method for masking by recursively supressing all connecting bright pixels to the mean value
# Specfically used for comet tails and other non-circular sources

def floodfill(data, x, y, mean, threshold, visited):
	if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False:
		data[y][x] = mean
		visited[y][x] = True

	else:
		return data

	floodfill(data, x - 1, y, mean, threshold, visited)
	floodfill(data, x + 1, y, mean, threshold, visited)
	floodfill(data, x, y - 1, mean, threshold, visited)
	floodfill(data, x, y + 1, mean, threshold, visited)

	return data




# Specifically checks for cut-off bright sources on the boundaries of the image i.e. centroid not on cutout

def perimeter(data, mean, stddev):
	visited = np.zeros((len(data), len(data)), dtype=bool)
	for i in range(101):
		if data[len(data) - 1][i] > mean + 3 * stddev:
			data = floodfill(data, i, len(data) - 1, mean, mean + stddev, visited)
		if data[i][len(data) - 1] > mean + 3 * stddev:
			data = floodfill(data, len(data) - 1, i, mean, mean + stddev, visited)
		if data[0][i] > mean + 3 * stddev:
			data = floodfill(data, i, 0, mean, mean + stddev, visited)
		if data[i][0] > mean + 3 * stddev:
			data = floodfill(data, 0, i, mean, mean + stddev, visited)
			
	return data




# Normalizes a PSF by recursively checking all connecting bright points greater than a specified and treshold (2 sigma) and scaling them

def normalize(data, x, y, p1, p2, threshold, visited):
	if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False:
		data[y][x] /= p1
		data[y][x] *= p2
		visited[y][x] = True
	else:
		return data

	normalize(data, x - 1, y, p1, p2, threshold, visited)
	normalize(data, x + 1, y, p1, p2, threshold, visited)
	normalize(data, x, y - 1, p1, p2, threshold, visited)
	normalize(data, x, y + 1, p1, p2, threshold, visited)

	return data



# Stand in for Atlas files, determines which pixels are associated with each source by checking if greater than threshold (2 sigma)
# Returns a boolean array

def connected(data, x, y, threshold, visited, radius):
	if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False and distance(x, y, 50, 50) < radius:
		visited[y][x] = True
	else:
		return visited

	connected(data, x - 1, y, threshold, visited, radius)
	connected(data, x + 1, y, threshold, visited, radius)
	connected(data, x, y - 1, threshold, visited, radius)
	connected(data, x, y + 1, threshold, visited, radius)

	return visited



	
# Disect the objid to gain useful information like id

def objid_extract(obj_id, full=False):
	masks={'sky_version':0x7800000000000000,
		   'rerun':0x07FF000000000000,
		   'run':0x0000FFFF00000000,
		   'camcol':0x00000000E0000000,
		   'first_field':0x0000000010000000,
		   'field':0x000000000FFF0000,
		   'id':0x000000000000FFFF}

	run=(obj_id & masks['run']) >> 32
	rerun=(obj_id & masks['rerun']) >> 48
	camcol=(obj_id & masks['camcol']) >> 29
	field=(obj_id & masks['field']) >> 16
	id=(obj_id & masks['id']) >> 0
	sky_version=(obj_id & masks['sky_version']) >> 59
	first_field=(obj_id & masks['first_field']) >> 28

	return {'run':run,
			'rerun':rerun,
			'camcol':camcol,
			'field':field,
			'id':id,
			'first_field':first_field,
			'sky_version':sky_version}





# Determine all flags associated with an object

def detflags(total):
	flags = np.zeros(32, dtype=bool)

	if total < 0:
		total += 2<<31
		flags[31] = True

	for i in range(1, 32):
		if total > (2 << (31 - i)):
			total -= 2 << (31 - i)
			flags[31 - i] = True

	return flags





# Class that assists in choosing the optimal PSF sources by storing index and Z value

class tuplet:
	def __init__(self, i, z):
		self.i = i
		self.z = z
		
	def getSelf(z):
		return self

	def getIndex():
		return i

	def getZ():
		return z




# Main method for executing PSF Subtraction
# Only searches in a 600x600 box around the quasar for a PSF source
# If none are found then the file is passed on and not used

def begin(index):
	#i = int(index)
	i = int(index.split('-')[0])
	mgi = int(index.split('-')[1])
	color = index.split('-')[2]
	
	#try:
	print(index)
	#filename = 'Test Data Extract/' + str(i) + '.fit'
	#filename = str(i) + '-g.fit'


	try:
		# Open the Experiment Absorber's data
		qlist = fits.open('/data/marvels/billzhu/Experiment Absorbers/' + color + '/' + str(index) + '_EA.fit')
		abs_data = qlist[0].data.astype(float)


		# Open a reference QSO's data
		qsofile = fits.open('/data/marvels/billzhu/2175 Reference Quasar Cut/' + color + '/' + str(index) + '_REF.fit')
		qsodata = qsofile[0].data.astype(float)

		# Open the reference QSO's field image
		scifile = fits.open('/data/marvels/billzhu/2175 Reference Dataset/' + color + '/' + str(index) + '.fit')
		scidata = scifile[0].data.astype(float)

		# Open the object binary table
		obj_table = Table.read('/data/marvels/billzhu/2175 Reference Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1)
	except:
		return

	
	pointer = 0
	if color == 'g':
		pointer = 1
	if color == 'r':
		pointer = 2
	if color == 'i':
		pointer = 3
	if color == 'z':
		pointer = 4
	if color == 'u':
		pointer = 0


	
	#line_data = linecache.getline('Full Data.txt', i).split()
	#line_data = linecache.getline('DR12 QSO.txt', i).split()
	#print(len(line_data))
	#obj_id = int(line_data[52])



	# Recalculate the sigma stats after background subtraction
	#mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
	#print("%f, %f, %f" % (mean, median, std))

	largearr = []
	stars = []
	chunk_size = 50


	qx = qsofile[0].header['XCOORD']
	qy = qsofile[0].header['YCOORD']
	obj_id = qsofile[0].header['ID']

	for j in range(len(obj_table)):
		sx = obj_table['COLC'][j][pointer]
		sy = obj_table['ROWC'][j][pointer]
		flags1 = detflags(obj_table['OBJC_FLAGS'][j])
		flags2 = detflags(obj_table['OBJC_FLAGS2'][j])

		#try:
		if obj_table['OBJC_TYPE'][j] == 6 and flags1[12] == False and flags1[17] == False and flags1[18] == False and flags2[27] == False and distance(sx, sy, qx, qy) > 5 and inbounds(sx + chunk_size + 6, sy + chunk_size + 6) and inbounds(sx - chunk_size - 5, sy - chunk_size - 5) and obj_table['PSFMAG'][j][pointer] < 18 and obj_table['M_RR_CC'][j][pointer] > 0 and abs(obj_table['M_RR_CC'][j][pointer] - obj_table['M_RR_CC'][obj_id - 1][pointer]) < 0.1 * obj_table['M_RR_CC'][obj_id - 1][pointer]:
			
			#try:
			"""
			preshift = scidata[int(sy - 10) : int(sy + 11), int(sx - 10) : int(sx + 11)]
			xc, yc = centroid_2dg(preshift, mask = None)
			xc += quasar['COLC'][pointer] - 10
			yc += quasar['ROWC'][pointer] - 10
			"""

			xc = obj_table['COLC'][j][pointer]
			yc = obj_table['ROWC'][j][pointer]
			preshift = scidata[int(yc - chunk_size - 5) : int(yc + chunk_size + 6), int(xc - chunk_size - 5) : int(xc + chunk_size + 6)]
			spline = interpolate.interp2d(np.arange(int(xc - chunk_size - 5), int(xc + chunk_size + 6)),
										  np.arange(int(yc - chunk_size - 5), int(yc + chunk_size + 6)),
										  preshift)

			xrang = np.arange(xc - chunk_size, xc + chunk_size + 1)
			yrang = np.arange(yc - chunk_size, yc + chunk_size + 1)

			if len(xrang) > 2 * chunk_size + 1:
				xrang = xrang[:-1]
			if len(yrang) > 2 * chunk_size + 1:
				yrang = yrang[:-1]

			shifted1 = spline(xrang, yrang)


			# Subtract off accurate local sky background
			shifted1 -= sigma_clipped_stats(shifted1, sigma=3.0, iters=5)[1]


			"""
			background = []

			for r in range(len(shifted1)):
				for c in range(len(shifted1)):
					if distance(r, c, 150, 150) > 150:
						background.append(shifted1[r, c])

			true_bkg = sigma_clipped_stats(background, sigma=3.0, iters=5)[1]
			shifted1 -= true_bkg
			"""

			# Take out the PSF star's image, and any absorbers associated with it if very close

			mean, median, std = sigma_clipped_stats(shifted1, sigma=3.0, iters=5)

			
			visited = connected(shifted1, chunk_size, chunk_size, mean + std/4, np.zeros((2 * chunk_size + 1, 2 * chunk_size + 1), dtype=bool), obj_table['PETROTHETA'][j][pointer] / 0.396 * 1.2 + 3)
			"""
			poststamp = np.zeros((2 * chunk_size + 1, 2 * chunk_size + 1))
			for r in range(len(shifted1)):
				for c in range(len(shifted1)):
					if visited[c][r] == True:
						poststamp[c][r] = shifted1[c][r]
		 
			shifted1 -= poststamp
			"""



			"""
			mag22 = 0
			if color == 'g':
				mag22 = 2500
			if color == 'r':
				mag22 = 1900
			if color == 'i':
				mag22 = 1400
			if color == 'z':
				mag22 = 300
			"""

			#print("%f, %f" % (xc, yc))
			shifted1 = checkInner(shifted1, obj_table, xc, yc, mean, std, pointer)
			shifted1 = perimeter(shifted1, mean, std)
			#shifted1 += poststamp

			#gauss1 = photutils.fit_2dgaussian(shifted1[40 : 61, 40 : 61], mask = None)
			#fit_fwhm_x = 2*np.sqrt(2*np.log(2))*np.sqrt(abs(gauss1.x_stddev.value))
			#fit_fwhm_y = 2*np.sqrt(2*np.log(2))*np.sqrt(abs(gauss1.y_stddev.value))
			

			#print("%f, %f, %f, %f" % (fit_fwhm_x, fit_fwhm_y, fwhm_x, fwhm_y))


			#if abs(fit_fwhm_x - fwhm_x) < 0.25 and abs(fit_fwhm_y - fwhm_y) < 0.25:
			#print(shifted1[50][50])



			if shifted1[50][50] > 0.05:

				#shifted1 = normalize(shifted1, 50, 50, np.max(shifted1), np.max(qsodata), mean1 + 5 * std1, np.zeros((101, 101), dtype=bool))

				#visited = connected(shifted1, 50, 50, mean + 1 * std, np.zeros((101, 101), dtype=bool))
				#shifted1 = gaussian_filter(shifted1, sigma=std)
				#smax = np.max(shifted1)
				
				for r in range(len(visited)):
					for c in range(len(visited)):
						#if distance(r, c, 50, 50) < 1.2 * obj_table['iso_a'][obj_id - 1][pointer]/2 + 3:
						shifted1[c][r] /= obj_table['PSFFLUX'][j][pointer]
						shifted1[c][r] *= obj_table['PSFFLUX'][obj_id - 1][pointer]
				
				#shifted1 /= obj_table['PSFFLUX'][j][pointer]
				#shifted1 *= obj_table['PSFFLUX'][obj_id - 1][pointer]

				
				largearr.append(np.reshape(shifted1, 10201))
				stars.append(j)
				#print(True)

				#stars.append(shifted1)
		#except:
			#continue

	largearr = np.array(largearr)
	#largearr *= 1 + 10**((obj_table['PSFMAG'][obj_id - 1][pointer] - 24.5) / 2.5)
	#largearr /= photonCount(50, 50, obj_table['PETROTHETA'][obj_id - 1][pointer] / 0.396 , qsodata)
	#largearr *= photonCount(50, 50, obj_table['PETROTHETA'][obj_id - 1][pointer] / 0.396 , qsodata_copy)
	print(np.shape(largearr))



	# Set number of components in PCA, use Incremental PCA (IPCA) due to high efficiency and speed	
	numcomp = len(largearr)


	# Need a healthy number of sources to make the PSF fitting in order to decrease noise, setting at 5% threshold

	#if len(largearr) < 10:
	#    print('No Sources')
	#    return

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


	# Calculate the residues with offsets from 1 to 8 pixels

	pixel = np.arange(0, 8)

	for offset in pixel:	
		qsodata_copy = np.array(qsodata)
		abs_data_copy = np.array(abs_data)

		print(qlist[0].header['PSFMAG'])
		factor = 10 ** (qlist[0].header['PSFMAG'] / -2.5 - 19.5 / -2.5)
		abs_data_copy /= (factor * 10)
		qsodata_copy[:, offset : 101] += abs_data_copy[:, 0 : 101 - offset]




		# Final fitting of the first n components, as determined by take_final, into the quasar to build a PSF fit
		#print(np.shape(ipca_comp))
		qsodata_copy = np.reshape(qsodata_copy, (2 * chunk_size + 1)**2)
		qsodata_copy -= mean_vector
		qsodata_copy = np.reshape(qsodata_copy, (2 * chunk_size + 1, 2 * chunk_size + 1))
		coeff = np.dot(np.reshape(qsodata_copy[47 : 54, 47 : 54], 49), new_comp)
		
		#coeff = np.dot(qsodata, ipca_comp)

		final_fit = np.dot(ipca_comp[:, 0:take_final], coeff[0:take_final])
		final_fit += mean_vector
		final_fit = np.reshape(final_fit, (2 * chunk_size + 1, 2 * chunk_size + 1))
		#final_fit /= len(largearr)

		qsodata_copy = np.reshape(qsodata_copy, (2 * chunk_size + 1)**2)
		qsodata_copy += mean_vector
		qsodata_copy = np.reshape(qsodata_copy, (2 * chunk_size + 1, 2 * chunk_size + 1))


		"""
		background = []

		for r in range(len(final_fit)):
			for c in range(len(final_fit)):
				if distance(r, c, 50, 50) > 50:
					background.append(final_fit[r, c])
		"""

		mean, median, stddev = sigma_clipped_stats(final_fit)
		final_fit -= median


		#final_fit /= final_fit[50, 50]
		#final_fit *= qsodata_copy[50, 50]
		final_fit /= photonCount(50, 50, 2 * obj_table['PSF_FWHM'][obj_id - 1][pointer], final_fit)
		final_fit *= photonCount(50, 50, 2 * obj_table['PSF_FWHM'][obj_id - 1][pointer], qsodata_copy)

		print("%f, %f" % (np.max(qsodata_copy), np.max(final_fit)))




		# Final residue from subtraction of PSF from QSO

		residue = qsodata_copy - final_fit
		mean, median, stddev = sigma_clipped_stats(residue, sigma=3.0, iters=5)
		residue -= median
		

		try:
			fits.writeto('/data/marvels/billzhu/Experiment PSF Subtract/' + color + '/' + str(index) + '-' + str(offset) + '_SUB.fit', residue, qlist[0].header, clobber = True)
			fits.writeto('/data/marvels/billzhu/Experiment PSF Cut/' + color + '/' + str(index) + '-' + str(offset) + '_PSF.fit', final_fit, qlist[0].header, clobber = True)

			print('\n')

			print("DONE TO BOTTOM")
		except:
			print('HEADER IS CORRUPT')


	
	
# Code that opens up a maximum of 8 processes for concurrent execution
	
if __name__ == '__main__':
	#multiprocessing.set_start_method('spawn')

	#gdirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/g/')
	#dirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/')




	gdirs = os.listdir('/data/marvels/billzhu/Experiment Absorbers/g/')
	#rdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/r/')
	#idirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/z/')
	#udirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/u/')


	#gdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/i/')
	#rdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/r/')
	#idirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/z/')
	
	check_dirs = os.listdir('/data/marvels/billzhu/Experiment PSF Subtract/g/')
	rangelist = []


	#rangelist.append(gdirs[0].split('_')[0])
	
	for d in gdirs:
		index = d.split('_')[0]
		#if str(index) + '_1' + '_SUB.fit' not in check_dirs:
		rangelist.append(index)


	
	#for d in rdirs:
	#	index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
	#	rangelist.append(index)
	
	#for d in idirs:
	#	index = d.split('_')[0]
	#	if str(index) + '_SUB.fit' not in check_dirs:
	#		rangelist.append(index)
	
	#for d in zdirs:
	#	index = d.split('_')[0]
	#	#if str(index) + '_SUB.fit' not in check_dirs:
	#	rangelist.append(index)
	
	
	print(len(rangelist))
	#begin(rangelist[0])
	#try:
	pool = Pool(os.cpu_count())
	pool.map(begin, rangelist)
	#except:
	#print("Error: Unable to process file")





