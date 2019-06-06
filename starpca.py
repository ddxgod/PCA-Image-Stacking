#Last updated: 9/23/2017
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
from astropy.cosmology import WMAP9 as cosmo
import astropy.units as u
import astropy.coordinates as coord
import sys
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




# Calculate the background by taking the sigma clipped median value between the radii specified

def calc_background(data, x, y, radius1, radius2):
	bkg_array = []
	for i in range(len(data[0])):
		for j in range(len(data)):
			if abs(x - i) < radius2 and abs(y - j) < radius2 and inbounds(i, j) and (i - x)**2 + (j - y)**2 >= radius1**2 and (i - x)**2 + (j - y)**2 <= radius2**2:
				bkg_array.append(data[j, i])

	true_bkg = sigma_clipped_stats(bkg_array, sigma=3.0, iters=5)[0]
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
			

			#else:
			data = floodfill(data, int(sources['colc'][i][pointer] - qX + qX1), int(sources['rowc'][i][pointer] - qY + qY1), mean, mean + stddev, visited)


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




# Main method for executing PSF Subtraction
# Only searches in a 600x600 box around the quasar for a PSF source
# If none are found then the file is passed on and not used

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
	
	scidata = hdulist[0].data.astype(float)


	#print(sigma_clipped_stats(scidata, sigma=3.0, iters=5))

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


	try:
		#print('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '.fit')
		#obj_table = Table.read('/data/marvels/billzhu/2175 Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1)
		#obj_table = Table.read('/data/marvels/billzhu/2175 Reference Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1)
		obj_table = Table.read('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '_Obj.fit', hdu=1)
		#obj_table = Table.read('/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(i) + '.fit', hdu=1)
	except:
		print(str(i) + ' No Table')
		return
	
	#line_data = linecache.getline('Full Data.txt', i).split()
	#line_data = linecache.getline('DR12 QSO.txt', i).split()
	#print(len(line_data))
	#obj_id = int(line_data[52])

	qlist = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')

	qx = qlist[0].header['XCOORD']
	qy = qlist[0].header['YCOORD']
	obj_id = qlist[0].header['ID']
	mean1 = qlist[0].header['MEAN_BKG']
	median1 = qlist[0].header['MED_BKG']
	std1 = qlist[0].header['STD_BKG']
	
	#print("%f, %f" % (x, y))
	qlist.close()


	print("%d, %f" % (obj_id, obj_table['M_rr_cc'][obj_id - 1][pointer]))



	# Calculate the 18 magnitude threshold

	mag18 = 0
	header = hdulist[0].header


	
	try:
		mag18 = header['FLUX20'] * 10 **(8. - 18/2.5)
	except:
		if color == 'g':
			mag18 = 15000
		if color == 'r':
			mag18 = 10500
		if color == 'i':
			mag18 = 8800
		if color == 'z':
			mag18 = 1900
		print(str(i) + ' MAG20 APPROX = 15000')
	

	print('reached')


	largearr = []
	stars = []
	chunk_size = 50
	diff_fwhm = 1000000


	counter = 0
	

	#qsocount = photonCount(50, 50, 15, star_image)
	#mean1, median1, std1 = sigma_clipped_stats(scidata, sigma=3.0, iters=5)

	chosen = -1

	RA = float(linecache.getline('Full Data.txt', i).split()[1])
	DEC = float(linecache.getline('Full Data.txt', i).split()[2])
	redshift = 0

	for j in range(len(mgtable)):
		if abs(mgtable['RA'][j] - RA) < 0.0001 and mgtable['DEC'][j] - DEC < 0.0001 and mgtable['ZABS'][j][0] >= 0.37 and mgtable['ZABS'][j][0] < 0.55:
			redshift = mgtable['ZABS'][j][0]
			break

	if redshift == 0:
		return

	scale = cosmo.kpc_proper_per_arcmin(redshift) * u.arcmin / u.kiloparsec * 0.396 / 60


	for j in range(len(obj_table)):
		sx = obj_table['colc'][j][pointer]
		sy = obj_table['rowc'][j][pointer]
		flags1 = detflags(obj_table['objc_flags'][j])
		flags2 = detflags(obj_table['objc_flags2'][j])

		#try:
		if j != obj_id - 1 and obj_table['objc_type'][j] == 6 and flags1[12] == False and flags1[17] == False and flags1[18] == False and flags2[27] == False and distance(sx, sy, qx, qy) > 5 and inbounds(sx + chunk_size + 6, sy + chunk_size + 6) and inbounds(sx - chunk_size - 5, sy - chunk_size - 5) and obj_table['psfCounts'][j][pointer] > mag18 and obj_table['M_rr_cc'][j][pointer] > 0 and abs(obj_table['M_rr_cc'][j][pointer] - obj_table['M_rr_cc'][obj_id - 1][pointer]) < 0.1 * obj_table['M_rr_cc'][obj_id - 1][pointer]:
			
			print(j)
			#try:
			qx = obj_table['colc'][j][pointer]
			qy = obj_table['rowc'][j][pointer]
			#mean1, median1, std1 = sigma_clipped_stats(scidata, sigma=3.0, iters=5)

			#preshift = scidata[int(yc - chunk_size - 5) : int(yc + chunk_size + 6), int(xc - chunk_size - 5) : int(xc + chunk_size + 6)]
			xu = qx + 350.0
			xl = qx - 350.0
			yu = qy + 350.0
			yl = qy - 350.0
			xc1 = qx
			yc1 = qy

			if xu + 350 >= 2048:
				xu = 2047

			if xl - 350 < 0:
				xl = 0
			else:
				xc1 = qx - int(qx) + 350.0

			if yu + 350 >= 1489:
				yu = 1488

			if yl - 350 < 0:
				yl = 0
			else:
				yc1 = qy - int(qy) + 350.0


			scidata2 = np.array(scidata[int(yl) : int(yu), int(xl) : int(xu)])
			visited = np.zeros((len(scidata2), len(scidata2[0])), dtype=bool)
			scidata2 = checkInner(scidata2, obj_table, qx, qy, xc1, yc1, mean1, std1, visited, pointer)
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

			star_image = spline(xrang, yrang)

			if star_image[chunk_size, chunk_size] < 20:
				continue


			chosen = j

			print('reached point 1')
			break

			#except:
				#continue

	if chosen == -1:
		print("NO SOURCES FOUND")
		return



	visited = connected(star_image, 50, 50, mean1 + 1 * std1, np.zeros((101, 101), dtype=bool), 3 * obj_table['iso_a'][chosen][pointer] / 0.396 + 3)
	#print(np.shape(star_image))


	for j in range(len(obj_table)):
		sx = obj_table['colc'][j][pointer]
		sy = obj_table['rowc'][j][pointer]
		flags2 = detflags(obj_table['objc_flags2'][j])
		flags1 = detflags(obj_table['objc_flags'][j])

		#try:
		if j != chosen and j != obj_id - 1 and obj_table['objc_type'][j] == 6 and flags1[12] == False and flags1[17] == False and flags1[18] == False and flags2[27] == False and distance(sx, sy, qx, qy) > 5 and inbounds(sx + chunk_size + 6, sy + chunk_size + 6) and inbounds(sx - chunk_size - 5, sy - chunk_size - 5) and obj_table['psfCounts'][j][pointer] > mag18 and obj_table['M_rr_cc'][j][pointer] > 0 and abs(obj_table['M_rr_cc'][j][pointer] - obj_table['M_rr_cc'][chosen][pointer]) < 0.1 * obj_table['M_rr_cc'][chosen][pointer]:


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
			visited1 = np.zeros((len(scidata2), len(scidata2[0])), dtype=bool)
			scidata2 = checkInner(scidata2, obj_table, xc, yc, xc1, yc1, mean1, std1, visited1, pointer)
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


				#shifted1 = normalize(shifted1, 50, 50, np.max(shifted1), np.max(star_image), mean1 + 5 * std1, np.zeros((101, 101), dtype=bool))
				
				#print(visited)
				#shifted1 = gaussian_filter(shifted1, sigma=std)
				#smax = np.max(shifted1)
				#objcount = photonCount(chunk_size, chunk_size, 15, shifted1)

				#fits.writeto('Test' + str(j) + '.fit', shifted1, clobber=True)
				
				
				for r in range(len(visited)):
					for c in range(len(visited)):
						if visited[c, r] == True:
							shifted1[c, r] /= obj_table['psfCounts'][j][pointer]
							shifted1[c, r] *= obj_table['psfCounts'][chosen][pointer]
				

				#shifted1 /= objcount#obj_table['psfCounts'][j][pointer]
				#shifted1 *= qsocount#obj_table['psfCounts'][obj_id - 1][pointer]
				#fits.writeto('Test00' + str(counter) + '.fit', scidata2, clobber=True)
				#print("%d, %f, %f" % (counter, obj_table['M_rr_cc'][j][pointer], obj_table['iso_a'][j][pointer]))
				#scounter += 1

				#shifted1 -= np.median(shifted1)
				
				
				largearr.append(np.reshape(shifted1, (2 * chunk_size + 1)**2))
				#print(True)

				#stars.append(shifted1)


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
	star_image1 = np.reshape(star_image, (2 * chunk_size + 1)**2)
	star_image1 -= mean_vector
	star_image1 = np.reshape(star_image1, (2 * chunk_size + 1, 2 * chunk_size + 1))
	coeff = np.dot(np.reshape(star_image1[chunk_size - 6 : chunk_size + 7, chunk_size - 6 : chunk_size + 7], 169), new_comp)
	
	#coeff = np.dot(star_image, ipca_comp)

	final_fit = np.dot(ipca_comp[:, 0:take_final], coeff[0:take_final])
	final_fit += mean_vector
	final_fit = np.reshape(final_fit, (2 * chunk_size + 1, 2 * chunk_size + 1))
	#final_fit /= len(largearr)

	star_image1 = np.reshape(star_image1, (2 * chunk_size + 1)**2)
	star_image1 += mean_vector
	star_image1 = np.reshape(star_image1, (2 * chunk_size + 1, 2 * chunk_size + 1))
	

	print("%f, %f" % (np.max(star_image), np.max(final_fit)))



	# Final residue from subtraction of PSF from QSO


	residue = star_image - final_fit
	#mean, median, stddev = sigma_clipped_stats(residue, sigma=3.0, iters=10)
	#residue -= median
	

	try:
		#fits.writeto('/data/marvels/billzhu/2175 PSF Cut/' + color + '/'+ str(index) + '_PSF.fit', final_fit, qsocut[0].header, clobber = True)
		#fits.writeto('/data/marvels/billzhu/2175 PSF Subtract/' + color + '/' + str(index) + '_SUB.fit', residue, qsocut[0].header, clobber = True)
		#fits.writeto('/data/marvels/billzhu/2175 Reference PSF Cut/' + color + '/'+ str(index) + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
		#fits.writeto('/data/marvels/billzhu/2175 Reference PSF Subtract/' + color + '/' + str(index) + '_SUB.fit', residue, hdulist[0].header, clobber = True)

		#fits.writeto('/data/marvels/billzhu/Reference PSF Cut/0.37 - 0.55/' + color + '/' + str(index) + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
		#fits.writeto('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + color + '/' + str(index) + '_SUB.fit', residue, hdulist[0].header, clobber = True)
		fits.writeto('/data/marvels/billzhu/Star PCA Cut/' + color + '/' + str(i) + '-' + color + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
		fits.writeto('/data/marvels/billzhu/Star PCA Subtract/' + color + '/' + str(i) + '-' + color + '_SUB.fit', residue, hdulist[0].header, clobber = True)

		
		#fits.writeto('Reference Subtract/' + str(i) + '_SUB.fit', residue, hdulist[0].header, clobber = True)
		#fits.writeto('Reference PSF Cut/' + str(i) + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
		print('\n')

		print("DONE TO BOTTOM")
	except:
		print('HEADER IS CORRUPT')


	
	
# Code that opens up a maximum of 8 processes for concurrent execution
	
if __name__ == '__main__':
	#multiprocessing.set_start_method('spawn')

	mgtable = Table.read('QSObased_Trimmed_SDSS_DR7_107.fits')
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


	#gdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/i/')
	#rdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/r/')
	#idirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/z/')
	
	check_dirs = os.listdir('/data/marvels/billzhu/Star PCA Subtract/r/')
	rangelist = []
	#rangelist.append('99795-r')
	#rangelist.append('89-230830-g')

	#rangelist.append('100728-g')

	#rangelist.append('96034-g')

	#rangelist.append('100009-i')
	#rangelist.append(rdirs[0].split('_')[0])



	
	for d in gdirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	
	"""
	for d in rdirs:
		index = d.split('_')[0]
		if str(index) + '_SUB.fit' not in check_dirs:
			rangelist.append(index)
	
	
	
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





