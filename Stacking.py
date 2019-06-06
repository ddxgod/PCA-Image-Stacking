from astropy.io import fits
from astropy.table import Column
from astropy.table import Table
import scipy as sp
import scipy.optimize as opt
import numpy as np
from numpy import *
import string
import decimal
import matplotlib.pyplot as plt
import fileinput
import math
import linecache


from astropy.stats import mad_std
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import os
import os.path




# Calculates the distance btween two given points

def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)



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


 #dirs = os.listdir('MG II Subtract/')


# Binn the data by Signal of 2796 A line

"""
table_list = Table.read('QSObased_Trimmed_SDSS_DR7_107.fits')
dirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/')
refcutdirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/')


for j in range(1, 7):
	bin1 = []
	bin2 = []
	bin3 = []
	refbin1 = []
	refbin2 = []
	refbin3 = []
	lower_bound = round(0.37 + 0.03 * (float(j) - 1.), 2)
	upper_bound = round(0.37 + 0.03 * float(j), 2)

	print(len(table_list))
	for i in range(len(table_list)):
		if table_list['REW_MGII_2796'][i][0] >= 0.8 and table_list['REW_MGII_2796'][i][0] < 1.12 and table_list['ZABS'][i][0] >= lower_bound and table_list['ZABS'][i][0] < upper_bound and (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit') in dirs:

			if 'FLUX20' not in fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'))[0].header.keys():
				continue
				#print(str(table_list['INDEX_QSO'][i] + 1))
			
			bin1.append(i)
			
			for k in range(len(refcutdirs)):
				if str(table_list['INDEX_QSO'][i] + 1) in refcutdirs[k] and 'FLUX20' in fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + refcutdirs[k])[0].header.keys():
					refbin1.append(refcutdirs[k].split('_')[0])


		if table_list['REW_MGII_2796'][i][0] >= 1.12 and table_list['REW_MGII_2796'][i][0] < 1.58 and table_list['ZABS'][i][0] >= lower_bound and table_list['ZABS'][i][0] < upper_bound and (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit') in dirs:

			if 'FLUX20' not in fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'))[0].header.keys():
				continue         

			bin2.append(i)

			for k in range(len(refcutdirs)):
				if str(table_list['INDEX_QSO'][i] + 1) in refcutdirs[k] and 'FLUX20' in fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + refcutdirs[k])[0].header.keys():
					refbin2.append(refcutdirs[k].split('_')[0])


		if table_list['REW_MGII_2796'][i][0] >= 1.58 and table_list['ZABS'][i][0] >= lower_bound and table_list['ZABS'][i][0] < upper_bound and (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit') in dirs:

			if 'FLUX20' not in fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + (str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'))[0].header.keys():
				continue        

			bin3.append(i)

			for k in range(len(refcutdirs)):
				if str(table_list['INDEX_QSO'][i] + 1) in refcutdirs[k] and 'FLUX20' in fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + refcutdirs[k])[0].header.keys():
					refbin3.append(refcutdirs[k].split('_')[0])

					

	scitot080_112 = np.zeros((42, 42, len(bin1)))
	scitot112_158 = np.zeros((42, 42, len(bin2)))
	scitot158 = np.zeros((42, 42, len(bin3)))


	counter = 0
	for i in bin1:
		print(str(table_list['INDEX_QSO'][i] + 1))
		filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'
		hdulist = fits.open(filename)
		scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
		
		scitot080_112[:, :, counter] = scidata
		counter += 1

	counter = 0
	for i in bin2:
		filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'
		hdulist = fits.open(filename)
		scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
		scitot112_158[:, :, counter] = scidata
		counter += 1
			
	counter = 0
	for i in bin3:
		filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + str(table_list['INDEX_QSO'][i] + 1) + '_SUB.fit'
		hdulist = fits.open(filename)
		scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
		scitot158[:, :, counter] = scidata
		counter += 1

	
	reftot080_112 = np.zeros((42, 42, len(refbin1)))
	reftot112_158 = np.zeros((42, 42, len(refbin2)))
	reftot158 = np.zeros((42, 42, len(refbin3)))

	for i in range(len(refbin1)):
		filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + refbin1[i] + '_SUB.fit'
		hdulist = fits.open(filename)
		scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
		reftot080_112[:, :, i] = scidata[0 : 42, 0 : 42]

	for i in range(len(refbin2)):
		filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + refbin2[i] + '_SUB.fit'
		hdulist = fits.open(filename)
		scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
		reftot112_158[:, :, i] = scidata[0 : 42, 0 : 42]

	for i in range(len(refbin3)):
		filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + refbin3[i] + '_SUB.fit'
		hdulist = fits.open(filename)
		scidata = hdulist[0].data.astype(float) / hdulist[0].header['FLUX20'] * 2000.
		reftot158[:, :, i] = scidata[0 : 42, 0 : 42]
	 

	scitot080_112_final = np.zeros((42, 42))
	scitot112_158_final = np.zeros((42, 42))
	scitot158_final = np.zeros((42, 42))
	reftot080_112_final = np.zeros((42, 42))
	reftot112_158_final = np.zeros((42, 42))
	reftot158_final = np.zeros((42, 42))
	
	for i in range(42):
		for k in range(42):
			scitot080_112_final[i][k] = np.mean(scitot080_112[i][k][:])
			scitot112_158_final[i][k] = np.mean(scitot112_158[i][k][:])
			scitot158_final[i][k] = np.mean(scitot158[i][k][:])
			reftot080_112_final[i][k] = np.mean(reftot080_112[i][k][:])
			reftot112_158_final[i][k] = np.mean(reftot112_158[i][k][:])
			reftot158_final[i][k] = np.mean(reftot158[i][k][:])

			
	fits.writeto('080_112_' + str(lower_bound) + '_' + str(upper_bound) + '_MGcomb.fit', scitot080_112_final, clobber = True)
	fits.writeto('112_158_' + str(lower_bound) + '_' + str(upper_bound) + '_MGcomb.fit', scitot112_158_final, clobber = True)
	fits.writeto('158_' + str(lower_bound) + '_' + str(upper_bound) + '_MGcomb.fit', scitot158_final, clobber = True)

	
	fits.writeto('080_112_' + str(lower_bound) + '_' + str(upper_bound) + '_REFcomb.fit', reftot080_112_final, clobber = True)
	fits.writeto('112_158_' + str(lower_bound) + '_' + str(upper_bound) + '_REFcomb.fit', reftot112_158_final, clobber = True)
	fits.writeto('158_' + str(lower_bound) + '_' + str(upper_bound) + '_REFcomb.fit', reftot158_final, clobber = True)

	list1 = list()
	header1 = fits.Header(list1)
	header1.append(('MGIIQSO', len(bin1), 'Number of MG II QSO stacked, not number of absorbers'))
	header1.append(('REFQSO', len(refbin1), 'Number of REF QSO stacked'))
	fits.writeto('/data/marvels/billzhu/stacked/080_112_' + str(lower_bound) + '_' + str(upper_bound) + '_SUBcomb.fit', scitot080_112_final - reftot080_112_final, header1, clobber = True)

	list1 = list()
	header1 = fits.Header(list1)
	header1.append(('MGIIQSO', len(bin2), 'Number of MG II QSO stacked, not number of absorbers'))
	header1.append(('REFQSO', len(refbin2), 'Number of REF QSO stacked'))
	fits.writeto('/data/marvels/billzhu/stacked/112_158_' + str(lower_bound) + '_' + str(upper_bound) + '_SUBcomb.fit', scitot112_158_final - reftot112_158_final, header1, clobber = True)

	list1 = list()
	header1 = fits.Header(list1)
	header1.append(('MGIIQSO', len(bin3), 'Number of MG II QSO stacked, not number of absorbers'))
	header1.append(('REFQS)', len(refbin3), 'Number of REF QSO stacked'))
	fits.writeto('/data/marvels/billzhu/stacked/158_' + str(lower_bound) + '_' + str(upper_bound) + '_SUBcomb.fit', scitot158_final - reftot158_final, header1, clobber = True)
	


		



"""    

dirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/i/')
#dirs = os.listdir('/data/marvels/billzhu/Star PCA Subtract/g/')
#dirs = os.listdir('/data/marvels/billzhu/Simulated PSF Subtract/r14')
#print(dirs)
#check_dirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/g/')
#print(check_dirs)
count = len(dirs)
print(str(count))
#scitot1 = np.zeros((42, 42, len(dirs)))

#scitot = np.zeros((8, 101, 101))
scitot = np.zeros((101, 101))

#scitot2 = np.zeros((42, 42, len(dirs)//2))
j = 0

#tableDR12 = Table.read('DR12Q.fits')

#counter = np.zeros(8)
counter = 0
mean = []
num = 0
for i in range(len(dirs)):
	#print(str(i))

	#print(dirs[i])
	#if dirs[i].split('_')[0] + '_SUB.fit' not in check_dirs:
	#    continue


	#try:
	filename = '/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/i/' + dirs[i]
	#filename = '/data/marvels/billzhu/Star PCA Subtract/g/' + dirs[i]
	#filename = '/data/marvels/billzhu/Simulated PSF Subtract/r14/' + dirs[i]
	hdulist = fits.open(filename)
	#obj_id = hdulist[0].header['ID']
	scidata = hdulist[0].data.astype(float)

	"""
	if np.max(scidata[43 : 58, 43 : 58]) > 25:
		print(dirs[i])
		continue
	
	
	
	bkg_sigma = mad_std(scidata)
	daofind = DAOStarFinder(fwhm = 1.4 / 0.4, threshold=5.*bkg_sigma)
	sources = daofind.find_stars(scidata)

	mean1, median1, stddev1 = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
	visited = np.zeros((len(scidata), len(scidata[0])), dtype=bool)
	print("%s, %d, %f" % (dirs[i], len(sources), bkg_sigma))
	for k in range(len(sources)):
		#if sources['peak'][k] > 3 * bkg_sigma:
		scidata = floodfill(scidata, int(sources['xcentroid'][k]), int(sources['ycentroid'][k]), mean1, mean1 + stddev1, visited)
	


	#if scidata[69, 62] > 50:
	#	print(dirs[i])
	#	continue
	
	
	
	mean1, median1, stddev1 = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
	
	for r in range(len(scidata)):
		for c in range(len(scidata[0])):
			if scidata[r,c] < mean1 - 3 * stddev1 or (scidata[r,c] > mean1 + 3 * stddev1):# and distance(r, c, 50, 50) > 8):
				scidata[r,c] = mean1
	


	
	#if np.mean(scidata) < 0 or scidata[0][23] > 20:
	#	print("%s, %f" % (dirs[i], np.mean(scidata)))
	#	continue
	
	
	
	#index = dirs[i].split('_')[0]
	#print(index)
	#j = int(index.split('-')[0])
	#mgi = int(index.split('-')[1])
	#color = index.split('-')[2]
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

	print(str(j) + '-' + str(mgi) + '.fit')
	obj_table = Table.read('/data/marvels/billzhu/2175 MG II Obj/' + str(j) + '-' + str(mgi) + '.fit', hdu=1)
	print('Reached')
	scidata *= 10**(0.4 * obj_table['EXTINCTION'][obj_id - 1][pointer])
	scitot += scidata
	counter += 1
	"""

	color = 'i'
	table = Table.read('/data/marvels/billzhu/MG II psField/0.37 - 0.55/' + dirs[i].split('-')[0] + '_psField.fit', hdu=7)
	print(counter)
	scidata /= table['flux20'][3]
	mean.append(table['flux20'][3])
	num += 1
	scidata *= 2000
	"""
	if 'FLUX20' not in hdulist[0].header.keys():
		scidata /= 2500
		scidata *= 2000
	else:
		scidata /= hdulist[0].header['FLUX20']
		#print(hdulist[0].header['FLUX20'])
		mean.append(hdulist[0].header['FLUX20'])
		num += 1
		scidata *= 2000
	"""
	#if sigma_clipped_stats(scidata)[0] < 0:
	#    print(dirs[i])
	#    continue

	#scidata *= hdulist[0].header['NMGY']
	
	
	linedata = linecache.getline('Full Data.txt', int(dirs[i].split('-')[0])).split()

	
	

	multiplier = 0
	if color == 'g':
		multiplier = 10 ** (0.4 * 0.736 * float(linedata[14]))
	if color == 'r':
		multiplier = 10 ** (0.4 * 0.534 * float(linedata[14]))
	if color == 'i':
		multiplier = 10 ** (0.4 * 0.405 * float(linedata[14]))
	if color == 'z':
		multiplier = 10 ** (0.4 * 0.287 * float(linedata[14]))

	scidata *= multiplier
	


	"""	
	if '0_SUB' in dirs[i]:
		scitot[0, :, :] += scidata
		counter[0] += 1
	if '1_SUB' in dirs[i]:
		scitot[1, :, :] += scidata
		counter[1] += 1
	if '2_SUB' in dirs[i]:
		scitot[2, :, :] += scidata
		counter[2] += 1
	if '3_SUB' in dirs[i]:
		scitot[3, :, :] += scidata
		counter[3] += 1
	if '4_SUB' in dirs[i]:
		scitot[4, :, :] += scidata
		counter[4] += 1
	if '5_SUB' in dirs[i]:
		scitot[5, :, :] += scidata
		counter[5] += 1
	if '6_SUB' in dirs[i]:
		scitot[6, :, :] += scidata
		counter[6] += 1
	if '7_SUB' in dirs[i]:
		scitot[7, :, :] += scidata
		counter[7] += 1
	"""
	
	#else:
	#    print("%s, %f" % (dirs[i], np.max(scidata)))


	

	
	if np.shape(scitot) != np.shape(scidata):
		continue
	
	scitot += scidata
	counter += 1
	
	
	#except:
		#print("Error: File not found")
		#j += 1

#print(np.median(mean))
#print(counter)


#print(scitot)

#medcomb = np.zeros((8, 101, 101))
medcomb = np.zeros((101, 101))
#medcomb2 = np.zeros((42, 42))


print(counter)

#for k in range(8):
for i in range(len(medcomb[0])):
	for j in range(len(medcomb[0])):
		#mean, median, stddev = sigma_clipped_stats(scitot[i][j][:], sigma_lower=-3.0, iters=5)
		#print("%f, %f, %f" % (mean, median, stddev))
		medcomb[i, j] += scitot[i, j] / counter


print(np.mean(mean))

#fits.writeto('RefHalf1.fit', medcomb, clobber = True)
#fits.writeto('RefHalf2.fit', medcomb2, clobber = True)

#for k in range(8):
fits.writeto('MGSUBComb62-i.fit', medcomb, overwrite = True)


####### NOTEEE
# g-band: 2543.2761, 2478.4668
# r-band: 1969.1451, 1952.2144
# i-band: 1458.7325, 1439.3639
# z-band: 310.18903, 300.73694


