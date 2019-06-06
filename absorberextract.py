#Last updated: 9/23/2017
#Description: Second version of qso cutting that uses the fpObj file and quasar catalogs for optimal centroid measurements and locating in the field images.

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
import ezgal
sys.setrecursionlimit(15000000)

import sdsspy
from sdsspy import atlas
from sdsspy import yanny



# Calculates the distance btween two given points

def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)



def inbounds(x, y):
	return x < 2048 and x >= 0 and y < 1489 and y >= 0



# Method that checks if a source > 1 sigma is outside of a certain radius, and if so, masks it by putting mean value

def checkOutter(data, mean, std):
	count = 0
	for i in range(len(data)):
		for j in range(len(data[0])):
			if data[i][j] > mean + 3 * std and distance(i, j, 50, 50) > 50:
				data[i][j] = mean
	return data



# Checks if a sources > 2 sigma is inside of a certain radius, and if so, masks the source by calling floodfill()

def checkInner(data, sources, qX, qY, mean, stddev, pointer, gthresh):
	count = 0

	for i in range(len(sources)):
		
		"""
		gal = 0
		if sources['type'][i][1] == 3:
			gal += 1
		if sources['type'][i][2] == 3:
			gal += 1
		if sources['type'][i][3] == 3:
			gal += 1


		
		if sources['objc_type'][i] == 3 and sources['psfCounts'][i][pointer] < gthresh:
			#print(gthresh)
			#print(sources['psfCounts'][i][pointer])
			#print(pointer)
			continue
		"""


		if distance(sources['COLC'][i][pointer] - qX + 50, sources['ROWC'][i][pointer] - qY + 50, 50, 50) > 5 and distance(sources['COLC'][i][pointer] - qX + 50, sources['ROWC'][i][pointer] - qY + 50, 50, 50) < 100:
			#print(True)
			sx = int(sources['COLC'][i][pointer])
			sy = int(sources['ROWC'][i][pointer])

			#if sources['psfCounts'][i][pointer] > 10000:
				#print("%f, %f, %f" % (sources['objc_type'][i], sources['COLC'][i][pointer], sources['ROWC'][i][pointer]))


			visited = np.zeros((len(data), len(data[0])), dtype=bool)
			

			# Radial masking with isophotal radii is ineffective, so discontinued
			"""
			if sources['iso_a'][i][pointer] > 0:
				data = radiusmask(data, int(sources['COLC'][i][pointer] - qX + 50), int(sources['ROWC'][i][pointer] - qY + 50), sources['iso_a'][i][pointer], mean, visited)
			
			#print("%f, %f" % (sources['COLC'][i][pointer] - qX + 50, sources['ROWC'][i][pointer] - qY + 50))
			else:
			"""

			data = floodfill(data, sx, sy, mean, mean + stddev, visited)

			#print(True)
			#print("%f, %f" % (sources['COLC'][i][pointer] - qX + 50, sources['ROWC'][i][pointer] - qY + 50))

	return data




# If the isophotal major axis is given, then use it for accurate masking

def radiusmask(data, xc, yc, isophotal_radius, mean, visited):
	if xc < 0 or xc >= len(data) or yc < 0 or yc >= len(data[0]):
		return data

	for r in range(len(data)):
		for c in range(len(data)):
			if distance(r, c, xc, yc) < 1.2 * isophotal_radius / 2 + 3:
				data[c][r] = mean

	return data




# Specifically checks for cut-off bright sources on the boundaries of the image i.e. centroid not on cutout

def perimeter(data, mean, stddev):
	visited = np.zeros((len(data), len(data)), dtype=bool)
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



			
# Main method for masking by redcursively supressing all connecting bright pixels to the mean value

def floodfill(data, x, y, mean, threshold, visited):
	if x >= 0 and x < len(data[0]) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False:
		data[y][x] = mean
		visited[y][x] = True
	else:
		return data

	floodfill(data, x - 1, y, mean, threshold, visited)
	floodfill(data, x + 1, y, mean, threshold, visited)
	floodfill(data, x, y - 1, mean, threshold, visited)
	floodfill(data, x, y + 1, mean, threshold, visited)

	return data





# Stand in for Atlas files, determines which pixels are associated with each source by checking if greater than threshold (2 sigma)
# Returns a boolean array

def connected(data, x, y, threshold, visited, radius):
	if x >= 0 and x < len(data) and y >= 0 and y < len(data) and data[y][x] >= threshold and visited[y][x] == False and distance(x, y, 50, 50) <= radius:
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




def begin(index):
	print(index)
	i = int(index.split('-')[0])
	mgi = int(index.split('-')[1])
	color = index.split('-')[2]

	try:
		hdulist = fits.open('/data/marvels/billzhu/2175 Reference Dataset/' + color + '/' + str(index) + '.fit')
		#hdulist = fits.open('/data/marvels/billzhu/2175 Dataset/' + color + '/' + str(index) + '.fit')
		#hdulist = fits.open('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '.fit')
		#hdulist = fits.open('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + str(mgi) + '-' + color + '.fit')
		scidata = hdulist[0].data.astype(float)

		obj_table = Table.read('/data/marvels/billzhu/2175 Reference Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1) 
		#obj_table = Table.read('/data/marvels/billzhu/2175 Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1) 
		#obj_table = Table.read('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '_Obj.fit', hdu=1)
		
		#obj_table = Table.read('/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(i) + '_Obj.fit', hdu=1)
		test = obj_table['COLC'][0]
	except:
		print(str(i) + ' Can not get table')
		return
	
	#line_data = linecache.getline('Full Data.txt', i).split()
	#obj_id = int(line_data[52])


	
	print(dr12_table[mgi]['OBJ_ID'])
	obj_id = objid_extract(int(dr12_table[mgi]['OBJ_ID']), False)['id']
	print("%f, %f" % (obj_id, len(obj_table)))
	if obj_id >= len(obj_table):
		print('obj_id over table length')
		return
	

	
	#quasar = obj_table[obj_id - 1]

	
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

	chunk_size = 50
	


	#try:
		#imdict = atlas.atlas.read_atlas(filename='/data/marvels/billzhu/MG II Atlas/0.37 - 0.55/' + str(i) + '_Atlas.fit', id=obj_id, trim=False)
		#imdict = atlas.atlas.read_atlas(filename='/data/marvels/billzhu/Reference Atlas/0.37 - 0.55/' + str(i) + '_Atlas.fit', id=obj_id, trim=False)
	#except:
	#	return


	#fits.writeto('AtlasTest.fit', imdict['images'][pointer], clobber=True)
	#print(imdict.keys())


	"""
	poststamp = imdict['images'][pointer].astype(float)

	for r in range(len(poststamp)):
		for c in range(len(poststamp[0])):
			if poststamp[r][c] == 1000:
				poststamp[r][c] = 0


	#fits.writeto('TestAtlas1.fit', poststamp, clobber=True)
	#print(np.shape(poststamp))
	
	#qx, qy = np.unravel_index(poststamp.argmax(), np.shape(poststamp))
	#scidata[]
	row0 = imdict['row0'][pointer]
	col0 = imdict['col0'][pointer]
	#print(imdict)
	row, col = np.shape(poststamp)
	xc = quasar['COLC'][pointer]
	yc = quasar['ROWC'][pointer]
	image = scidata[int(yc - chunk_size - 5) : int(yc + chunk_size + 6), int(xc - chunk_size - 5) : int(xc + chunk_size + 6)]
	mean1, median1, std1 = sigma_clipped_stats(image, sigma=3.0, iters=5)
	scidata[row0 : row0 + row, col0 : col0 + col] -= poststamp
	"""

	#true_bkg = calc_background(scidata, xc, yc, 200, 250)
	#scidata -= true_bkg

	#print("%f, %f" % (quasar['COLC'][pointer], quasar['ROWC'][pointer]))

	#qsofile = fits.open('/data/marvels/billzhu/2175 Reference Quasar Cut/' + color + '/' + str(index) + '_REF.fit')
	#qsodata = qsofile[0].data.astype(float)


	for j in range(len(obj_table)):
		flags1 = detflags(obj_table['OBJC_FLAGS'][j])
		flags2 = detflags(obj_table['OBJC_FLAGS2'][j])

		if obj_table['OBJC_TYPE'][j] == 3 and obj_table['PARENT'][j] == -1 and obj_table['NCHILD'][j] == 0 and obj_table['PSFMAG'][j][pointer] > 19 and obj_table['PSFMAG'][j][pointer] < 22 and obj_table['PETROTH90'][j][pointer] > 0 and obj_table['PETROTH90'][j][pointer] < obj_table['PETROTH90'][obj_id - 1][pointer] and inbounds(obj_table['COLC'][j][pointer] + chunk_size + 6, obj_table['ROWC'][j][pointer] + chunk_size + 6) and inbounds(obj_table['COLC'][j][pointer] - chunk_size - 5, obj_table['ROWC'][j][pointer] - chunk_size - 5):
			
			#image = scidata[int(quasar['ROWC'][pointer] - 10) : int(quasar['ROWC'][pointer] + 11), int(quasar['COLC'][pointer] - 10) : int(quasar['COLC'][pointer] + 11)]

			#xc, yc = centroid_2dg(image, mask = None)
			#print("%f, %f, %f, %f" % (quasar['COLC'][pointer], quasar['ROWC'][pointer], xc, yc))
			#xc += int(quasar['COLC'][pointer] - 10)
			#yc += int(quasar['ROWC'][pointer] - 10)

			try:
				imdict = atlas.atlas.read_atlas(filename='/data/marvels/billzhu/2175 Reference Atlas/' + str(i) + '-' + str(mgi) + '_Atlas.fit', id=j + 1, trim=False)

				poststamp = imdict['images'][pointer].astype(float)

				for r in range(len(poststamp)):
					for c in range(len(poststamp[0])):
						if poststamp[r, c] < 1005:
							poststamp[r, c] = 1000
				poststamp -= 1000


				#fits.writeto('Test2175Atlas.fit', poststamp, clobber=True)
				row0 = imdict['row0'][pointer]
				col0 = imdict['col0'][pointer]
				row, col = np.shape(poststamp)
				#newimage = scidata[row0 : row0 + row, col0 : col0 + col]
				scidata[row0 : row0 + row, col0 : col0 + col] -= poststamp
				print(np.shape(scidata))
				region = np.zeros(np.shape(scidata), dtype=bool)
				print(np.shape(region))
				region[row0 : row0 + row, col0 : col0 + col] = True
			except:
				continue

			"""
			for r in range(len(region):
				for c in range(len(region[0])):
					if region[r, c] == False:
						#print('Reached')
						scidata[r, c] = 0
			"""

			#fits.writeto('TestAtlas1.fit', scidata, clobber=True)
				#print(np.shape(poststamp))
			#except:
				#print('Problem with Atlas')
				#continue

			xc = obj_table['COLC'][j][pointer]
			yc = obj_table['ROWC'][j][pointer]
			image = scidata[int(yc - chunk_size - 5) : int(yc + chunk_size + 6), int(xc - chunk_size - 5) : int(xc + chunk_size + 6)]
			mean1, median1, std1 = sigma_clipped_stats(image, sigma=3.0, iters=5)


			print("%f, %f, %f" % (obj_table['PSFMAG'][j][pointer], obj_table['COLC'][j][pointer], obj_table['ROWC'][j][pointer]))
			print("%f, %f, %f" % (mean1, median1, std1))

			
			# Add back the QSO's poststamp image
			#scidata[row0 : row0 + row, col0 : col0 + col] += poststamp
			scidata = checkInner(scidata, obj_table, xc, yc, mean1, std1, pointer, 100)
			scidata[row0 : row0 + row, col0 : col0 + col] += poststamp

			# Interpolate the image to the centroid of the QSO
			image = scidata[int(yc - chunk_size - 5) : int(yc + chunk_size + 6), int(xc - chunk_size - 5) : int(xc + chunk_size + 6)]
			print(sigma_clipped_stats(image)[1])
			spline = interpolate.interp2d(np.arange(int(xc - chunk_size - 5), int(xc + chunk_size + 6)),
										  np.arange(int(yc - chunk_size - 5), int(yc + chunk_size + 6)),
										  image)
			xrang = np.arange(xc - chunk_size, xc + chunk_size + 1)
			yrang = np.arange(yc - chunk_size, yc + chunk_size + 1)

			if len(xrang) > 2 * chunk_size + 1:
				xrang = xrang[:-1]
			if len(yrang) > 2 * chunk_size + 1:
				yrang = yrang[:-1]

			shifted = spline(xrang, yrang)
			
			shifted = perimeter(shifted, mean1, std1)
			mean1, median1, std1 = sigma_clipped_stats(shifted, sigma=3.0, iters=5)
			shifted -= median1

			for r in range(len(shifted)):
				for c in range(len(shifted[0])):
					if distance(r, c, 50, 50) > 2 * obj_table['PETROTH90'][j][pointer] / 0.396:
						shifted[r][c] = 0

			"""
			for r in range(len(shifted)):
				for c in range(len(shifted)):
					if distance(r, c, 50, 50) > obj_table['PETROTH90'][j][pointer]/0.396 * 1.2 and shifted[r, c] > mean1 + 2 * std1:
						shifted[r, c] = 0
			"""
			# Take out the Quasar's image, and any absorbers associated with it if very close
			# CHANGED: USING ATLAS FILES NOW

			"""
			qsovisited = connected(shifted, chunk_size, chunk_size, mean1 + std1, np.zeros((2 * chunk_size + 1, 2 * chunk_size + 1), dtype=bool), obj_table['iso_a'][obj_id - 1][pointer] * 0.6 + 3)
			poststamp = np.zeros((2 * chunk_size + 1, 2 * chunk_size + 1))
			for r in range(len(shifted)):
				for c in range(len(shifted)):
					if qsovisited[c][r] == True:
						poststamp[c][r] = shifted[c][r]


			#fits.writeto('Test57104pre.fit', shifted, clobber = True)
			#fits.writeto('Test57104.fit', poststamp, clobber = True)            
			shifted -= poststamp
			"""
			try:
				hdulist[0].header.append(('XCOORD', xc, 'x coordinate of quasar in image'), end = True)
				hdulist[0].header.append(('YCOORD', yc, 'y coordinate of quasar in image'), end = True)
				hdulist[0].header.append(('PSFMAG', obj_table['PSFMAG'][j][pointer], 'id of the quasar in the fpObj file'), end = True) 
			except:
				print(str(i) + ' Unable to get coords')
				#return


			try:
				#shifted /= hdulist[0].header['NMGY']
				fits.writeto('/data/marvels/billzhu/Experiment Absorbers/' + color + '/' + str(index) + '_EA.fit', shifted, hdulist[0].header, clobber = True)
				#fits.writeto('/data/marvels/billzhu/2175 Quasar Cut/' + color + '/' + str(index) + '_DUST.fit', shifted, hdulist[0].header, clobber = True)
				#fits.writeto('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit', shifted, hdulist[0].header, clobber = True)
				#fits.writeto('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + str(mgi) + '-' + color + '_REF.fit', shifted, hdulist[0].header, clobber = True)

				return

			except:
				print('Header is corrupt')
				return
		else:
			continue




# Main method that can be wired for multiprocessing purposes using Pool

if __name__ == '__main__':
	#multiprocessing.set_start_method('spawn')
	dust_table = Table.read('final_catalog_full.fit')
	dr12_table = Table.read('DR12Q.fits')


	#mgtable = Table.read('Trimmed_SDSS_DR7_107.fits')
	#reader = open('Full Data.txt', 'r')

	
	# Load ezgal model file with 100 Myr burst, low metallicity
	model = ezgal.model('www.baryons.org/ezgal/models/cb07_burst_0.1_z_0.05_chab.model')
	#model = ezgal.model('www.baryons.org/ezgal/models/bc03_ssp_z_0.02_chab.model')

	# Desired formation redshift
	zf = 4.0
	# Fetch an array of redshifts out to given formation redshift
	zs = np.arange(0.005, zf, 0.005)
	#print(np.shape(zs))
	model.set_normalization('sloan_g', 0.00000000000000001, -24.4)
	thresh_mags = model.get_apparent_mags(zf, filters=['sloan_u', 'sloan_g', 'sloan_r', 'sloan_i', 'sloan_z'], zs=zs)
	

	# Normalize to the SED of a galaxy with rest-frame absolute magnitude = -22.4
	# According to Ned Wright's UCLA CosmoCalc, at redshift 0.0072, the light travel time is 100 Myr, the burst time the model uses
	
	#print(thresh_mags[:, 2])


	
	#dirs = os.listdir('/data/marvels/billzhu/2175 Dataset/')
	#dirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/')

	#gdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/g/')
	#rdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/r/')
	#idirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/z/')
	#udirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/u/')


	#gdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/g/')
	#rdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/r/')
	#idirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/z/')

	
	gdirs = os.listdir('/data/marvels/billzhu/2175 Reference Dataset/g/')
	#rdirs = os.listdir('/data/marvels/billzhu/2175 Dataset/r/')
	#idirs = os.listdir('/data/marvels/billzhu/2175 Dataset/i/')
	#zdirs = os.listdir('/data/marvels/billzhu/2175 Dataset/z/')
	
	
	#gcheck_dirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/g/')
	#icheck_dirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/i/')
	#zcheck_dirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/z/')
	#print(len(check_dirs))
	


	rangelist = []
	#rangelist.append('42725-g')

	#rangelist.append('123-18494-g')
	#begin(rangelist[0])

	#rangelist.append('29066-g')

	
	for d in gdirs:
		index = d.split('.')[0]
		#if d.split('.')[0] + '_MG.fit' not in gcheck_dirs:
		rangelist.append(index)
		#begin(index)
	
	#print(len(rangelist))
	#begin(rangelist[1])
	
		#if d.split('-')[0] + '_MG.fit' not in check_dirs:
		#rangelist.append(index + '-r')
		#rangelist.append(index + '-i')
		#rangelist.append(index + '-z')
		#rangelist.append(index + '-u')


	"""
	for d in rdirs:
		index = d.split('.')[0]
		#if d.split('-')[0] + '_MG.fit' not in rcheck_dirs:
		rangelist.append(index)
		#begin(index)
	
	for d in idirs:
		index = d.split('.')[0]
		#if d.split('-')[0] + '_MG.fit' not in icheck_dirs:
		rangelist.append(index)
		#begin(index)

	for d in zdirs:
		index = d.split('.')[0]
		#if d.split('-')[0] + '_MG.fit' not in zcheck_dirs:
		rangelist.append(index)
		#begin(index)
	"""
	
	
	#for d in udirs:
	#    index = d.split('.')[0]
	#    #if d.split('-')[0] + '_MG.fit' not in check_dirs:
	#    rangelist.append(index)
	

	#print(len(rangelist))
	#try:
	#rlist = []
	#rlist.append('9807-g')
	pool = Pool(multiprocessing.cpu_count())
	pool.map(begin, rangelist)
	#begin(rangelist[0])
	#except:
	#print("Error: Unable to process file")

	print(True)
	
