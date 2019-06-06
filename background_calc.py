from astropy.io import fits
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



# Calculate the background by taking the sigma clipped median value between the radii specified

def calc_background_stats(data, x, y, radius1, radius2):
	bkg_array = []
	for i in range(len(data[0])):
		for j in range(len(data)):
			if abs(x - i) < radius2 and abs(y - j) < radius2 and inbounds(i, j) and (i - x)**2 + (j - y)**2 >= radius1**2 and (i - x)**2 + (j - y)**2 <= radius2**2:
				bkg_array.append(data[j, i])

	bkg_stats = sigma_clipped_stats(bkg_array, sigma=3.0, maxiters=5)

	return bkg_stats




# Calculate the mean, median, and stddev at 1 - 2 arcsec, 2 - 3 sec, etc. in individual frames and write to individual Astropy Table for futhur use

def calc_arcsec_stats(data, x, y, counter, color):
	arcsec_array = [[] for k in range(200)]

	for i in range(len(data[0])):
		for j in range(len(data)):
			#print("%d, %d" % (len(data[0]), len(data)))
			#print("%f, %f, %f, %f, %f" % (i, j, x, y, distance(i, j, x, y)))
			arcsec_array[int(math.floor(distance(i, j, x, y) * 0.396))].append(data[j, i])

	arcsec_stats = []
	mean_array = []
	median_array = []
	stddev_array = []

	for i in range(len(arcsec_array)):
		new_stats = sigma_clipped_stats(arcsec_array[i], sigma=3.0, maxiters=5)
		arcsec_stats.append(new_stats)
		mean_array.append(new_stats[0])
		median_array.append(new_stats[1])
		stddev_array.append(new_stats[2])

	mean_col = fits.Column(name='mean', format='D', unit='counts', array=mean_array)
	median_col = fits.Column(name='median', format='D', unit='counts', array=median_array)
	stddev_col = fits.Column(name='stddev', format='D', unit='counts', array=stddev_array)
	coldefs = fits.ColDefs([mean_col, median_col, stddev_col])
	#hdu = fits.BinTableHDU.from_columns([mean_col, median_col, stddev_col])

	return coldefs




# Checks if a sources > 2 sigma is inside of a certain radius, and if so, masks the source by calling floodfill()

def checkInner(data, sources, qX, qY, qX1, qY1, mean, stddev, visited, pointer):
	count = 0
	for i in range(len(sources)):


		if distance(sources['colc'][i][pointer], sources['rowc'][i][pointer], qX, qY) > 3 and distance(sources['colc'][i][pointer], sources['rowc'][i][pointer], qX, qY) < 350:

			data = floodfill(data, int(sources['colc'][i][pointer] - qX + qX1), int(sources['rowc'][i][pointer] - qY + qY1), mean, mean + stddev, visited)



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





# Main method that begins the processes

def begin(index):
	#i = int(index)
	i = int(index.split('-')[0])
	mgi = int(index.split('-')[1])
	color = index.split('-')[2]
	
	#try:
	print(index)
	#filename = 'Test Data Extract/' + str(i) + '.fit'
	#filename = str(i) + '-g.fit'


	#filename = '/data/marvels/billzhu/2175 Dataset/' + color + '/' + str(index) + '.fit'
	#filename = '/data/marvels/billzhu/2175 Reference Dataset/' + color + '/' + str(index) + '.fit'
	filename = '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/' + color + '/' + str(index) + '.fit'
	#filename = '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '.fit'
	hdulist = fits.open(filename)

	#qlist = fits.open('MG II Test Cut/' + str(i) + '_MG.fit')

	#qlist = fits.open('/data/marvels/billzhu/2175 Quasar Cut/' + color + '/' + str(index) + '_DUST.fit')
	#qlist = fits.open('/data/marvels/billzhu/2175 Reference Quasar Cut/' + color + '/' + str(index) + '_REF.fit')
	qlist = fits.open('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(index) + '_REF.fit')
	#qlist = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')

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
		#obj_table = Table.read('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '_Obj.fit', hdu=1)
		obj_table = Table.read('/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(i) + '_Obj.fit', hdu=1)
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


	
	pstable = Table.read('/data/marvels/billzhu/Reference psField/0.37 - 0.55/' + str(i) + '_psField.fit', hdu=7)
	mag20 = pstable['flux20']
	

	if color == 'g':
		mag18 = mag20[1] * 10 **(8. - 18/2.5)
	if color == 'r':
		mag18 = mag20[2] * 10 **(8. - 18/2.5)
	if color == 'i':
		mag18 = mag20[3] * 10 **(8. - 18/2.5)
	if color == 'z':
		mag18 = mag20[4] * 10 **(8. - 18/2.5)
	

	
	#qsocut = fits.open('/data/marvels/billzhu/2175 Quasar Cut/' + color + '/' + str(index) + '_DUST.fit')
	#qsocut = fits.open('/data/marvels/billzhu/2175 Reference Quasar Cut/' + color + '/' + str(index) + '_REF.fit')
	qsocut = fits.open('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(index) + '_REF.fit')
	#qsocut = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')
	qsodata = qsocut[0].data.astype(float)



	if qsodata[50, 50] < 5:
		return

	print('reached')


	largearr = []
	stars = []
	chunk_size = 50
	diff_fwhm = 1000000

	counter = 0


	scale = cosmo.kpc_proper_per_arcmin(redshift) * u.arcmin / u.kiloparsec * 0.396 / 60

	#filedir = '/data/marvels/billzhu/background/' + color + '/'
	filedir = '/data/marvels/billzhu/Reference Background/' + color + '/'


	for j in range(len(obj_table)):
		sx = obj_table['colc'][j][pointer]
		sy = obj_table['rowc'][j][pointer]
		flags1 = detflags(obj_table['objc_flags'][j])
		flags2 = detflags(obj_table['objc_flags2'][j])

		try:
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

				#print("%f, %f" % (xc, yc))
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

				#print("%f, %f, %f, %f, %f, %f" % (xc1, yc1, xl, yl, xu, yu))
				scidata2 = np.array(scidata[int(yl) : int(yu), int(xl) : int(xu)])
				visited = np.zeros((len(scidata2), len(scidata2[0])), dtype=bool)
				scidata2 = checkInner(scidata2, obj_table, xc, yc, xc1, yc1, mean1, std1, visited, pointer)
				bkg_stats = calc_background_stats(scidata2, xc1, yc1, int(400 / scale), int(500 / scale))
				print(bkg_stats)

				if scidata2[int(yc1), int(xc1)] - bkg_stats[0] > 20:
					bkg_stats_arr.append([bkg_stats])
					#print("%f, %f" % (xc1, yc1))
					coldefs = calc_arcsec_stats(scidata2, xc1, yc1, counter, color)
					bkg_col = fits.Column(name=str(round(float(400/scale * 0.396), 3)) + '-' + str(round(float(500/scale * 0.396), 3)) + ' arcsec bkg stats', format='D', unit='counts', array=bkg_stats)
					coldefs.add_col(bkg_col)

					hdu = fits.BinTableHDU.from_columns(coldefs)
					hdu.writeto(filedir + index + '_' + str(counter) + '.fits', overwrite=True)
					#print(counter)
					counter += 1

				else:
					print("ERASED")
		except:
			print('EXCEPTION')
			continue

	




if __name__ == '__main__':

	#gdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/g/')
	#rdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/r/')
	#idirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/2175 Reference Quasar Cut/z/')
	#dirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/')


	#gdirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/g/')
	#rdirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/r/')
	#idirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/z/')

	"""
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
	

	rangelist = []
	
	bkg_stats_arr = []

	#rangelist.append(gdirs[0].split('_')[0])
	#begin(rangelist[0])	

	"""
	for d in gdirs:
		index = d.split('_')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		rangelist.append(index)
	
	
	
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
	
	
	
	print(len(rangelist))
	#try:
	pool = Pool(multiprocessing.cpu_count())
	pool.map(begin, rangelist)

	"""
	writer = open('Mg II Background Statistics.txt', 'w')

	for i in range(len(bkg_stats_arr)):
		writer.write(bkg_stats_arr[i] + '\n')

	writer.close()
	"""
	#except:
	#print("Error: Unable to process file")