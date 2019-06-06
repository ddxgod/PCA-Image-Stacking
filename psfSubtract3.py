#Last updated: 1/14/2019
#Description: Final version of PCA PSF subtraction that uses the fpObj files given, gauranteeing optimal normalization, background + sky subtractions, only bright stars are used, etc. 

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
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
sys.setrecursionlimit(1500000)



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
	mean, median, std = sigma_clipped_stats(halo, sigma=3.0, maxmaxiters=5)

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




# Calculate the background by taking the sigma clipped median value between the radii specified

def calc_background(data, x, y, radius1, radius2):
	bkg_array = []
	for i in range(len(data[0])):
		for j in range(len(data)):
			if abs(x - i) < radius2 and abs(y - j) < radius2 and inbounds(i, j) and (i - x)**2 + (j - y)**2 >= radius1**2 and (i - x)**2 + (j - y)**2 <= radius2**2:
				bkg_array.append(data[j, i])

	true_bkg = sigma_clipped_stats(bkg_array, sigma=3.0, maxiters=5)[0]
	print(true_bkg)
	return true_bkg




# Checks if a sources > 2 sigma is inside of a certain radius, and if so, masks the source by calling floodfill()

def checkInner(data, sources, qX, qY, qX1, qY1, mean, stddev, visited, pointer):
	count = 0
	for i in range(len(sources)):


		#yes = False
		#if abs(sources['colc'][i][pointer] - qX - 3) < 3 and abs(sources['rowc'][i][pointer] - qY - 46) < 3:
		#    print(True)
		#    yes = True

		if distance(sources['colc'][i][pointer], sources['rowc'][i][pointer], qX, qY) > 3 and distance(sources['colc'][i][pointer], sources['rowc'][i][pointer], qX, qY) < 350:
			
			#visited = np.zeros((len(data), len(data[0])), dtype=bool)


			# Radial masking with the isophotal radii given is inefficient

			#if sources['iso_a'][i][pointer] > 0:
			# 	data = radiusmask(data, int(sources['colc'][i][pointer] - qX + 50), int(sources['rowc'][i][pointer] - qY + 50), sources['iso_a'][i][pointer], mean, visited)

			
			#print("%f, %f" % (sources['colc'][i][pointer] - qX + 50, sources['rowc'][i][pointer] - qY + 50))
			#else:
			data = floodfill(data, int(sources['colc'][i][pointer] - qX + qX1), int(sources['rowc'][i][pointer] - qY + qY1), mean, mean + stddev, visited)

			#print(True)
			#print("%f, %f" % (sources['colc'][i][pointer] - qX + 50, sources['rowc'][i][pointer] - qY + 50))


			# If no isophotal radii are given, then use radius = 6 pixels as an initial estimate

			"""
			for j in range(-6, 6):
				for k in range(-6, 6):
					if distance(sources['colc'][i][pointer], sources['rowc'][i][pointer], sources['colc'][i][pointer] + j, sources['rowc'][i][pointer] + k) < 6 and int(sources['rowc'][i][pointer] - qY + 50 + k) >= 0 and int(sources['rowc'][i][pointer] - qY + 50 + k) < 101 and int(sources['colc'][i][pointer] - qX + 50 + j) >= 0 and int(sources['colc'][i][pointer] - qX + 50 + j) < 101 and data[int(sources['rowc'][i][pointer] - qY + 50 + k)][int(sources['colc'][i][pointer] - qX + 50 + j)] > mean + 3 * stddev:
						data[int(sources['rowc'][i][pointer] - qY + 50 + k)][int(sources['colc'][i][pointer] - qX + 50 + j)] = mean
						#print(data[int(sources['colc'][i][pointer] - qX + 50 + j)][int(sources['rowc'][i][pointer] - qY + 50 + k)])
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
	if x >= 0 and x < len(data[0]) and y >= 0 and y < len(data) and data[y, x] >= threshold and visited[y, x] == False:
		data[y, x] = mean
		visited[y, x] = True

	else:
		return data

	floodfill(data, x - 1, y, mean, threshold, visited)
	floodfill(data, x + 1, y, mean, threshold, visited)
	floodfill(data, x, y - 1, mean, threshold, visited)
	floodfill(data, x, y + 1, mean, threshold, visited)

	return data




# Specifically checks for cut-off bright sources on the boundaries of the image i.e. centroid not on cutout

def perimeter(data, mean, stddev):
	visited = np.zeros((len(data), len(data[0])), dtype=bool)
	for i in range(len(data)):
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
	if x >= 0 and x < len(data[0]) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False:
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
	if x >= 0 and x < len(data[0]) and y >= 0 and y < len(data) and data[y, x] >= threshold and visited[y, x] == False and distance(x, y, 50, 50) < radius:
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




def begin(index):
	#i = int(index)
	i = int(index.split('-')[0])
	#mgi = int(index.split('-')[1])
	color = index.split('-')[1]
	
	#try:
	print(index)
	#filename = 'Test Data Extract/' + str(i) + '.fit'
	#filename = str(i) + '-g.fit'


	#filename = '/data/marvels/billzhu/2175 Dataset/' + color + '/' + str(index) + '.fit'
	#filename = '/data/marvels/billzhu/2175 Reference Dataset/' + color + '/' + str(index) + '.fit'
	#filename = '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/' + color + '/' + str(index) + '.fit'
	filename = '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '.fit'
	hdulist = fits.open(filename)

	#qlist = fits.open('MG II Test Cut/' + str(i) + '_MG.fit')

	#qlist = fits.open('/data/marvels/billzhu/2175 Quasar Cut/' + color + '/' + str(index) + '_DUST.fit')
	#qlist = fits.open('/data/marvels/billzhu/2175 Reference Quasar Cut/' + color + '/' + str(index) + '_REF.fit')
	#qlist = fits.open('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(index) + '_REF.fit')
	qlist = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')

	qx = qlist[0].header['XCOORD']
	qy = qlist[0].header['YCOORD']
	obj_id = qlist[0].header['ID']
	mean1 = qlist[0].header['MEAN_BKG']
	median1 = qlist[0].header['MED_BKG']
	std1 = qlist[0].header['STD_BKG']
	redshift = qlist[0].header['ZABS']
	
	#print("%f, %f" % (x, y))
	qlist.close()

	#except:
	#    print("No coordinates")
	#    return

	# Save some frickin time

		
	
	scidata = hdulist[0].data.astype(float)
	"""
	if 'SKY' in hdulist[0].header.keys():
		scidata -= float(hdulist[0].header['SOFTBIAS'])
		scidata -= float(hdulist[0].header['SKY'])
	else:
		mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, maxiters=5)
		scidata -= median
		print(str(i) + ' No sky')
		#return
	"""
	

	#print(sigma_clipped_stats(scidata, sigma=3.0, maxiters=5))

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


	#bkg_sigma = mad_std(scidata)

	try:
		#print('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '.fit')
		#obj_table = Table.read('/data/marvels/billzhu/2175 Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1)
		#obj_table = Table.read('/data/marvels/billzhu/2175 Reference Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1)
		obj_table = Table.read('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '_Obj.fit', hdu=1)
		#obj_table = Table.read('/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(i) + '_Obj.fit', hdu=1)
	except:
		print(str(i) + ' No Table')
		return
	
	#line_data = linecache.getline('Full Data.txt', i).split()
	#line_data = linecache.getline('DR12 QSO.txt', i).split()
	#print(len(line_data))
	#obj_id = int(line_data[52])
	quasar = obj_table[obj_id - 1]

	try:
		print("%d, %f" % (obj_id, obj_table['M_rr_cc'][obj_id - 1][pointer]))
	except:
		print("can't print table")
		return





	# If no quasar is found, the field image is deemed corrupt and not used
	  
	if quasar == 0:
		print(str(i) + ' No quasar')
		return


	# Calculate the 18 magnitude threshold

	mag18 = 0
	header = hdulist[0].header


	pstable = Table.read('/data/marvels/billzhu/MG II psField/0.37 - 0.55/' + str(i) + '_psField.fit', hdu=7)
	mag20 = pstable['flux20']
	

	if color == 'g':
		mag18 = mag20[1] * 10 **(8. - 18/2.5)
	if color == 'r':
		mag18 = mag20[2] * 10 **(8. - 18/2.5)
	if color == 'i':
		mag18 = mag20[3] * 10 **(8. - 18/2.5)
	if color == 'z':
		mag18 = mag20[4] * 10 **(8. - 18/2.5)
	

	#mag18 = 10**(4.5/2.5)

	
	#qsocut = fits.open('/data/marvels/billzhu/2175 Quasar Cut/' + color + '/' + str(index) + '_DUST.fit')
	#qsocut = fits.open('/data/marvels/billzhu/2175 Reference Quasar Cut/' + color + '/' + str(index) + '_REF.fit')
	#qsocut = fits.open('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(index) + '_REF.fit')
	qsocut = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')
	qsodata = qsocut[0].data.astype(float)

	
	"""
	if 'SKY' in hdulist[0].header.keys():
		qsodata -= float(hdulist[0].header['SOFTBIAS'])
		qsodata -= float(hdulist[0].header['SKY'])
	else:
		qsodata -= median
	"""

	#qsodata -= obj_table['sky'][obj_id - 1][pointer]
	#qsodata -= 1000


	if qsodata[50, 50] < 5:
		print('bad QSO')
		return

	print('reached')

	"""
	for r in range(len(qsodata)):
		for c in range(len(qsodata)):
			if distance(r, c, 150, 150) > 150:
				background.append(qsodata[r, c])
	"""
	#print(sigma_clipped_stats(qsodata, sigma=3.0, maxiters=5))
	

	# Recalculate the sigma stats after background subtraction
	#mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, maxiters=5)
	#print("%f, %f, %f" % (mean, median, std))

	largearr = []
	stars = []
	chunk_size = 50
	diff_fwhm = 1000000

	#psf_fwhm = 100
	#qsovisited = connected(qsodata, 50, 50, mean + 1 * std, np.zeros((101, 101), dtype=bool))
	#qmax = np.max(qsodata)
	counter = 0
	

	#qsocount = photonCount(50, 50, 15, qsodata)
	#mean1, median1, std1 = sigma_clipped_stats(scidata, sigma=3.0, maxiters=5)



	scale = cosmo.kpc_proper_per_arcmin(redshift) * u.arcmin / u.kiloparsec * 0.396 / 60


	for j in range(len(obj_table)):
		sx = obj_table['colc'][j][pointer]
		sy = obj_table['rowc'][j][pointer]
		flags1 = detflags(obj_table['objc_flags'][j])
		flags2 = detflags(obj_table['objc_flags2'][j])

		#try:
		if obj_table['objc_type'][j] == 6 and flags1[12] == False and flags1[17] == False and flags1[18] == False and flags2[27] == False and distance(sx, sy, qx, qy) > 5 and inbounds(sx + chunk_size + 6, sy + chunk_size + 6) and inbounds(sx - chunk_size - 5, sy - chunk_size - 5) and obj_table['psfCounts'][j][pointer] > mag18 and obj_table['M_rr_cc'][j][pointer] > 0 and abs(obj_table['M_rr_cc'][j][pointer] - obj_table['M_rr_cc'][obj_id - 1][pointer]) < 0.1 * obj_table['M_rr_cc'][obj_id - 1][pointer]:
			
			#try:
			"""
			preshift = scidata[int(sy - 10) : int(sy + 11), int(sx - 10) : int(sx + 11)]
			xc, yc = centroid_2dg(preshift, mask = None)
			xc += quasar['colc'][pointer] - 10
			yc += quasar['rowc'][pointer] - 10
			"""

			xc = obj_table['colc'][j][pointer]
			yc = obj_table['rowc'][j][pointer]

			#preshift = scidata[int(yc - chunk_size - 5) : int(yc + chunk_size + 6), int(xc - chunk_size - 5) : int(xc + chunk_size + 6)]

			xu = xc + 350.0
			xl = xc - 350.0
			yu = yc + 350.0
			yl = yc - 350.0
			xc1 = xc
			yc1 = yc

			if xu >= 2048:
				xu = 2047

			if xl < 0:
				xl = 0
			else:
				xc1 = xc - int(xc) + 350.0

			if yu >= 1489:
				yu = 1488

			if yl < 0:
				yl = 0
			else:
				yc1 = yc - int(yc) + 350.0


			scidata2 = np.array(scidata[int(yl) : int(yu), int(xl) : int(xu)])
			visited = np.zeros((len(scidata2), len(scidata2[0])), dtype=bool)
			scidata2 = checkInner(scidata2, obj_table, xc, yc, xc1, yc1, mean1, std1, visited, pointer)
			preshift = scidata2[int(yc1 - chunk_size - 5) : int(yc1 + chunk_size + 6), int(xc1 - chunk_size - 5) : int(xc1 + chunk_size + 6)]
			preshift -= calc_background(scidata2, xc1, yc1, int(400 / scale), int(500 / scale))
			spline = interpolate.interp2d(np.arange(int(xc1 - chunk_size - 5), int(xc1 + chunk_size + 6)),
										  np.arange(int(yc1 - chunk_size - 5), int(yc1 + chunk_size + 6)),
										  preshift)

			xrang = np.arange(xc1 - chunk_size, xc1 + chunk_size + 1)
			yrang = np.arange(yc1 - chunk_size, yc1 + chunk_size + 1)

			if len(xrang) > 2 * chunk_size + 1:
				xrang = xrang[:-1]
			if len(yrang) > 2 * chunk_size + 1:
				yrang = yrang[:-1]

			shifted1 = spline(xrang, yrang)


			if shifted1[chunk_size][chunk_size] > 20:

				#shifted1 = normalize(shifted1, 50, 50, np.max(shifted1), np.max(qsodata), mean1 + 5 * std1, np.zeros((101, 101), dtype=bool))
				mean, median, std = sigma_clipped_stats(shifted1, sigma=3.0, maxiters=5)

				visited = connected(shifted1, 50, 50, mean + 1 * std, np.zeros((101, 101), dtype=bool), 2.4 * obj_table['iso_a'][j][pointer] / 0.396 + 3)
				#print(visited)
				#shifted1 = gaussian_filter(shifted1, sigma=std)
				#smax = np.max(shifted1)
				#objcount = photonCount(chunk_size, chunk_size, 15, shifted1)

				#fits.writeto('Test' + str(j) + '.fit', shifted1, overwrite=True)
				
				
				for r in range(len(visited)):
					for c in range(len(visited)):
						#if distance(r, c, 50, 50) < 1.2 * obj_table['iso_a'][obj_id - 1][pointer]/2 + 3:
						shifted1[c, r] /= obj_table['psfCounts'][j][pointer]
						shifted1[c, r] *= obj_table['psfCounts'][obj_id - 1][pointer]
				

				#shifted1 /= objcount#obj_table['psfCounts'][j][pointer]
				#shifted1 *= qsocount#obj_table['psfCounts'][obj_id - 1][pointer]
				#fits.writeto('Test0' + str(counter) + '.fit', shifted1, overwrite=True)
				#print("%d, %f, %f" % (counter, obj_table['M_rr_cc'][j][pointer], obj_table['iso_a'][j][pointer]))
				#counter += 1

				#shifted1 -= np.median(shifted1)
				
				
				largearr.append(np.reshape(shifted1, (2 * chunk_size + 1)**2))
				#print(True)

				#stars.append(shifted1)

			else:
				print("ERASED")
		#except:
			#continue

	

	
	largearr = np.array(largearr)
	print(np.shape(largearr))


	# Set number of components in PCA, use incremental PCA (IPCA) due to high efficiency and speed
	
	numcomp = len(largearr)


	# Need a healthy number of sources to make the PSF fitting in order to decrease noise, setting at 5% threshold

	#if len(largearr) < 10:
	#    print('No Sources')
	#    return

	print(numcomp)
	mean_vector = []
 
	#print(np.shape(largearr))

	try:
		for j in range(0, (2 * chunk_size + 1)**2):
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
		new_comp.append(np.reshape(temp[chunk_size - 6 : chunk_size + 7, chunk_size - 6 : chunk_size + 7], 169))

	new_comp = np.array(new_comp)
	new_comp = new_comp.T
	print(np.shape(new_comp))

	ipca_comp = ipca_comp.T

	
	#print(ipca_comp)

	#print(np.shape(largearr[0, :]))
	#print(np.shape(ipca_comp))



	take_final = numcomp


	# Final fitting of the first n components, as determined by take_final, into the quasar to build a PSF fit
	print(np.shape(ipca_comp))
	qsodata1 = np.reshape(qsodata, (2 * chunk_size + 1)**2)
	qsodata1 -= mean_vector
	qsodata1 = np.reshape(qsodata1, (2 * chunk_size + 1, 2 * chunk_size + 1))
	coeff = np.dot(np.reshape(qsodata1[chunk_size - 6 : chunk_size + 7, chunk_size - 6 : chunk_size + 7], 169), new_comp)
	
	#coeff = np.dot(qsodata, ipca_comp)

	final_fit = np.dot(ipca_comp[:, 0:take_final], coeff[0:take_final])
	final_fit += mean_vector
	final_fit = np.reshape(final_fit, (2 * chunk_size + 1, 2 * chunk_size + 1))
	#final_fit /= len(largearr)

	qsodata1 = np.reshape(qsodata1, (2 * chunk_size + 1)**2)
	qsodata1 += mean_vector
	qsodata1 = np.reshape(qsodata1, (2 * chunk_size + 1, 2 * chunk_size + 1))


	"""
	background = []

	for r in range(len(final_fit)):
		for c in range(len(final_fit)):
			if distance(r, c, 50, 50) > 50:
				background.append(final_fit[r, c])
	"""

	#mean, median, stddev = sigma_clipped_stats(final_fit, sigma=3.0, maxiters=10)
	#final_fit -= median
		

	#final_fit = gaussian_filter(final_fit, sigma=stddev)

	"""
	fmax = np.max(final_fit)
	qmax = np.max(qsodata)
	for r in range(len(qsovisited)):
		for c in range(len(qsovisited)):
			if qsovisited[r][c] == True:
				final_fit[r][c] /= fmax
				final_fit[r][c] *= qmax

	print(np.min(final_fit))
	print(np.min(qsodata))
	"""


	
	# Section to normalize the PSF fitting by photon count, but this is unneccesary since CORRECT PCA fits better

	"""
	gauss_fit = photutils.fit_2dgaussian(final_fit[40 : 61, 40 : 61], mask = None)
	fit_fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(gauss_fit.x_stddev)
	#print(fit_fwhm)

	
	print("%f, %f" % (fwhm_x, fit_fwhm))
	print("%f, %f" % (np.max(qsodata), np.max(final_fit)))
	ffwhm = max(fwhm_x, fit_fwhm)
	ffphoton_1sig = photonCount(50, 50, 4 * fit_fwhm, final_fit) 
	qsophoton_1sig = photonCount(50, 50, 4 * fit_fwhm, qsodata)

	print("%f, %f" % (qsophoton_1sig, ffphoton_1sig))
	
	
	for j in range(len(final_fit)):
		for k in range(len(final_fit)):
			if distance(50, 50, j, k) < 4 * fit_fwhm:
				final_fit[j][k] /= ffphoton_1sig
				final_fit[j][k] *= qsophoton_1sig
	"""
	

	print("%f, %f" % (np.max(qsodata), np.max(final_fit)))

	#final_fit /= ffphoton_1sig
	#final_fit *= qsophoton_1sig

	
	"""
	line_data = linecache.getline('Full Data.txt', i).split()

	if color == 'g':
		mag = float(line_data[6])
	if color == 'r':
		mag = float(line_data[8])
	if color == 'i':
		mag = float(line_data[10])
	if color == 'z':
		mag = float(line_data[12])
	if color == 'u':
		mag = float(line_data[4])
		

	#try:
	multiplier = 10**(8 - mag / 2.5) * header['FLUX20']
	multiplier1 = quasar['psfCounts'][pointer]
	print(multiplier)
	print(str(multiplier1 - header['SKY']))
	final_fit /= ffphoton_1sig
	final_fit *= multiplier
	#except:
		#final_fit *= qsodata[50, 50]
		#return
	
	
	print("%f, %f" % (np.max(qsodata), np.max(final_fit)))
	
	"""


	# Final residue from subtraction of PSF from QSO


	residue = qsodata - final_fit
	#mean, median, stddev = sigma_clipped_stats(residue, sigma=3.0, maxiters=10)
	#residue -= median
	

	try:
		#fits.writeto('/data/marvels/billzhu/2175 PSF Cut/' + color + '/'+ str(index) + '_PSF.fit', final_fit, qsocut[0].header, overwrite = True)
		#fits.writeto('/data/marvels/billzhu/2175 PSF Subtract/' + color + '/' + str(index) + '_SUB.fit', residue, qsocut[0].header, overwrite = True)
		#fits.writeto('/data/marvels/billzhu/2175 Reference PSF Cut/' + color + '/'+ str(index) + '_PSF.fit', final_fit, hdulist[0].header, overwrite = True)
		#fits.writeto('/data/marvels/billzhu/2175 Reference PSF Subtract/' + color + '/' + str(index) + '_SUB.fit', residue, hdulist[0].header, overwrite = True)

		#fits.writeto('/data/marvels/billzhu/Reference PSF Cut/0.37 - 0.55/' + color + '/' + str(index) + '_PSF.fit', final_fit, hdulist[0].header, overwrite = True)
		#fits.writeto('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + color + '/' + str(index) + '_SUB.fit', residue, hdulist[0].header, overwrite = True)
		fits.writeto('/data/marvels/billzhu/MG II PSF Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_PSF.fit', final_fit, hdulist[0].header, overwrite = True)
		fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB.fit', residue, hdulist[0].header, overwrite = True)

		
		#fits.writeto('Reference Subtract/' + str(i) + '_SUB.fit', residue, hdulist[0].header, overwrite = True)
		#fits.writeto('Reference PSF Cut/' + str(i) + '_PSF.fit', final_fit, hdulist[0].header, overwrite = True)
		print('\n')

		print("DONE TO BOTTOM")
	except:
		print('HEADER IS CORRUPT')


	
	
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

	#check_dirs = os.listdir('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/r/')
	rangelist = []
	#rangelist.append('99795-r')
	#rangelist.append('89-230830-g')

	#rangelist.append('100728-g')

	#rangelist.append('96034-g')

	#rangelist.append('95823-g')
	#rangelist.append(gdirs[0].split('_')[0])
	

	
	for d in gdirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	
	
	
	for d in rdirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	
	
	
	for d in idirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	

	
	for d in zdirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	
	
	
	print(len(rangelist))
	#try:
	pool = Pool(multiprocessing.cpu_count())
	pool.map(begin, rangelist)
	#except:
	#print("Error: Unable to process file")





