from astropy.io import fits
import numpy as np
import scipy as sp
from astropy.stats import jackknife_resampling
from astropy.stats import jackknife_stats
import os
import math
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
import multiprocessing
from multiprocessing import Pool
import pickle




# Calculates the distance btween two given points

def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)



# Finds the mean of all photon counts where the value is 3 sigma above the mean

def photoncount(scidata, radius1, radius2):
	flux = 0
	length = 0
	#print(np.shape(scidata))
				
	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if distance(i, j, 50, 50) <= radius1 and distance(i, j, 50, 50) >= radius2:
				flux += scidata[i][j]# / 2000 / 10**8
				length += 1

	return flux / length



# Find the number of points within 10 kpc

def getlength(scidata, radius1, radius2):
	length = 0
	mean, median, stddev = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if distance(i, j, 50, 50) <= radius1 and distance(i, j, 50, 50) >= radius2:
				length += 1

	return length



def photometric(scidata):
	SB_array = []
	outter = []
	scale = cosmo.kpc_proper_per_arcmin(0.48) * u.arcmin / u.kiloparsec * 0.396 / 60
	print(scale)

	boundaries = [3, 10, 13, 16, 19, 27, 37, 51, 67, 100, 140]

	for j in range(len(boundaries) - 1):
		f = photoncount(scidata, boundaries[j + 1] / scale, boundaries[j] / scale)
		#print("%f, %f" % (boundaries[j + 1], boundaries[j]))
		#print(f)
			
		outter.append(boundaries[j + 1] / scale * 0.396)

		#mag = 22.5 - 2.5 * np.log10(f)
		mag = -2.5 * np.log10(f)
		if np.isnan(mag):
			mag = 38
		

		# Calculate the surface brightness by adding the amount of kpc, similar to traditional mg/arcsec^2 formula
		#surface_brightness = mag + 2.5 * math.log10(0.396**2)

		#print(surface_brightness)
		#SB_array.append(surface_brightness)
		SB_array.append(f)
		#print(SB_array)


	return SB_array





def SED(scidata, color):
	nmgy_count = 0

	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if distance(i, j, 50, 50) > 4 and distance(i, j, 50, 50) <= 42:
				nmgy_count += scidata[i, j] / 2000 / 10**8


	factor = 0

	if color == 'g':
		factor = 4770 * 10**-8
	if color == 'r':
		factor = 6231 * 10**-8
	if color == 'i':
		factor = 7625 * 10**-8
	if color == 'z':
		factor = 9134 * 10**-8

	return nmgy_count * 3631 * 10**-23 * 29979245800/(factor)**2 * 10**-8






def begin(i, stacked_data, all_data, length, color):
	resampled_data = (stacked_data - all_data[i, :, :]) / (length - 1)
	SB_array = photometric(resampled_data)
	print(SB_array)
	#SB_matrix.append(SB_array)
	sed = SED(resampled_data, color)
	print(sed)
	#SED_matrix.append(sed)


	return SB_array, sed





if __name__ == '__main__':

	gdirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/g/')
	rdirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/r/')
	idirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/i/')
	zdirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/z/')

	all_data_g = []
	all_data_r = []
	all_data_i = []
	all_data_z = []

	stacked_data_g = np.zeros((101, 101), dtype=float)
	stacked_data_r = np.zeros((101, 101), dtype=float)
	stacked_data_i = np.zeros((101, 101), dtype=float)
	stacked_data_z = np.zeros((101, 101), dtype=float)
	#print(np.shape(stacked_data))

	#print(rdirs)

	counter = 0



	"""
	For each of the four bands, find the QSOs that passed all four bands, then use jackknife resampling to build the SB profiles of each different stack in each band
	At the end of al SB profile calculations, combine the profile data for all four bands

	"""

	for i in range(len(gdirs)):
		index = gdirs[i].split('-')[0]
		#index = str(gdirs[i].split('-')[0] + '-' + gdirs[i].split('-')[1])
		#print(index)

		if index + '-r_SUB.fit' not in rdirs or index + '-i_SUB.fit' not in idirs or index + '-z_SUB.fit' not in zdirs:
			continue

		counter += 1
		hdulist_g = fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/g/' + gdirs[i])
		scidata_g = hdulist_g[0].data
		all_data_g.append(np.array(scidata_g))
		stacked_data_g += np.array(scidata_g)
		del scidata_g
		del hdulist_g[0].data
		hdulist_g.close()



	print(counter)

	stacked_data_g = np.array(stacked_data_g)
	all_data_g = np.array(all_data_g)
	print(np.shape(all_data_g))
	SB_matrix_g = []
	SED_matrix_g = []

	for i in range(counter):
		print("%s, %d" % ('g', i))
		SB_array, sed = begin(i, stacked_data_g, all_data_g, counter, 'g')
		SB_matrix_g.append(SB_array)
		SED_matrix_g.append(sed)

	del stacked_data_g
	del all_data_g
	





	for i in range(len(rdirs)):
		index = rdirs[i].split('-')[0]
		#index = str(rdirs[i].split('-')[0] + '-' + rdirs[i].split('-')[1])

		if index + '-g_SUB.fit' not in gdirs or index + '-i_SUB.fit' not in idirs or index + '-z_SUB.fit' not in zdirs:
			continue

		hdulist_r = fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/r/' + rdirs[i])
		scidata_r = hdulist_r[0].data
		all_data_r.append(np.array(scidata_r))
		stacked_data_r += np.array(scidata_r)
		del scidata_r
		del hdulist_r[0].data
		hdulist_r.close()



	stacked_data_r = np.array(stacked_data_r)
	all_data_r = np.array(all_data_r)
	print(np.shape(all_data_r))
	SB_matrix_r = []
	SED_matrix_r = []

	for i in range(counter):
		print("%s, %d" % ('r', i))
		SB_array, sed = begin(i, stacked_data_r, all_data_r, counter, 'r')
		SB_matrix_r.append(SB_array)
		SED_matrix_r.append(sed)

	del stacked_data_r
	del all_data_r






	for i in range(len(idirs)):
		index = idirs[i].split('-')[0]
		#index = str(idirs[i].split('-')[0] + '-' + idirs[i].split('-')[1])

		if index + '-g_SUB.fit' not in gdirs or index + '-r_SUB.fit' not in rdirs or index + '-z_SUB.fit' not in zdirs:
			continue

		hdulist_i = fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/i/' + idirs[i])
		scidata_i = hdulist_i[0].data
		all_data_i.append(np.array(scidata_i))
		stacked_data_i += np.array(scidata_i)
		del scidata_i
		del hdulist_i[0].data
		hdulist_i.close()



	stacked_data_i = np.array(stacked_data_i)
	all_data_i = np.array(all_data_i)
	print(np.shape(all_data_i))
	SB_matrix_i = []
	SED_matrix_i = []

	for i in range(counter):
		print("%s, %d" % ('i', i))
		SB_array, sed = begin(i, stacked_data_i, all_data_i, counter, 'i')
		SB_matrix_i.append(SB_array)
		SED_matrix_i.append(sed)

	del stacked_data_i
	del all_data_i





	for i in range(len(zdirs)):
		index = zdirs[i].split('-')[0]
		#index = str(zdirs[i].split('-')[0] + '-' + zdirs[i].split('-')[1])

		if index + '-g_SUB.fit' not in gdirs or index + '-r_SUB.fit' not in rdirs or index + '-i_SUB.fit' not in idirs:
			continue

		hdulist_z = fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/z/' + zdirs[i])
		scidata_z = hdulist_z[0].data
		all_data_z.append(np.array(scidata_z))
		stacked_data_z += np.array(scidata_z)
		del scidata_z
		del hdulist_z[0].data
		hdulist_z.close()



	stacked_data_z = np.array(stacked_data_z)
	all_data_z = np.array(all_data_z)
	print(np.shape(all_data_z))
	SB_matrix_z = []
	SED_matrix_z = []

	for i in range(counter):
		print("%s, %d" % ('z', i))
		SB_array, sed = begin(i, stacked_data_z, all_data_z, counter, 'z')
		SB_matrix_z.append(SB_array)
		SED_matrix_z.append(sed)

	del stacked_data_z
	del all_data_z


	


	SB_matrix_g = np.array(SB_matrix_g)
	SB_matrix_r = np.array(SB_matrix_r)
	SB_matrix_i = np.array(SB_matrix_i)
	SB_matrix_z = np.array(SB_matrix_z)

	SED_matrix_g = np.array(SED_matrix_g)
	SED_matrix_r = np.array(SED_matrix_r)
	SED_matrix_i = np.array(SED_matrix_i)
	SED_matrix_z = np.array(SED_matrix_z)

	matrix_full = []

	for i in range(len(SB_matrix_g)):
		matrix_full.append(list(SB_matrix_g[i]) + list(SB_matrix_r[i]) + list(SB_matrix_i[i]) + list(SB_matrix_z[i]))# + list(np.array([SED_matrix_g[i], SED_matrix_r[i], SED_matrix_i[i], SED_matrix_z[i]])))


	matrix_full = np.array(matrix_full)
	print(np.shape(matrix_full))
	matrix_full = matrix_full.T
	cov_matrix = np.cov(matrix_full)

	print(cov_matrix)

	file = open('MG II SB and SED Sample Covariance Matrix 2.txt', 'w')


	for i in range(len(cov_matrix)):
		file.write(str(cov_matrix[i]) + '\n')

	file.close()

	with open('MGJack.txt', 'wb') as fp:
		pickle.dump(cov_matrix, fp);

"""
for c in range(len(data[0])):
	for r in range(len(data[0][0])):
		resampled_data = jackknife_resampling(data[:, c, r])
		#test_statistic = np.mean
		#estimate, bias, stderr, conf_interval = jackknife_stats(data[:, c, r], test_statistic, 0.95)
		#print(stderr / 2000 / 10**8)
"""

