from shutil import copyfile
from astropy.io import fits
from astropy.table import Table
import os
import os.path
import numpy as np
import scipy as sp
import math
import string
import linecache
import multiprocessing
from multiprocessing import Pool

"""
for i in range(8373, 50001):
	filename = 'MG II Dataset/' + str(i) + '.fit'
	try:
		hdulist = fits.open(filename)
		counter = 0
		while counter < 4:
			
	except:
		print(str(i))
		try:
			copyfile('Final Data Extract/' + str(i) + '.fit', 'Reference Dataset/' + str(i) + '_REF.fit')
		except:
			print(str(i))
"""




def begin(index):
	print(index)
	color = index.split('-')[1]
	index = int(index.split('-')[0])
	line = linecache.getline('Full Data.txt', index)
	mgdata = line.split()
	redshift = float(mgdata[3])
	mgumag = float(mgdata[4])
	mggmag = float(mgdata[6])
	mgrmag = float(mgdata[8])
	mgimag = float(mgdata[10])
	mgzmag = float(mgdata[12])
	count = 0

	for i in range(1, 105784):
		if i == index or ismg[i] == True:
			continue

		
		line = linecache.getline('Full Data.txt', i)
		refdata = line.split()
		z = float(refdata[3])
		refumag = float(refdata[4])
		refgmag = float(refdata[6])
		refrmag = float(refdata[8])
		refimag = float(refdata[10])
		refzmag = float(refdata[12])

		if abs(redshift - z) < 0.1 and abs(mgumag - refumag) < 0.5 and abs(mggmag - refgmag) < 0.5 and abs(mgrmag - refrmag) < 0.5 and abs(mgimag - refimag) < 0.5 and abs(mgzmag - refzmag) < 0.5 and (i not in usedref):
			
			try:
				#if 'FLUX20' not in fits.open('/data/marvels/jerryxu/dr7/catalog/' + str(i) + '-g.fit')[0].header.keys():
				#    continue

				writer.write("%d %d\n" % (i, index))
				#copyfile('/data/marvels/jerryxu/dr7/raw_catalog/' + str(i) + '.fit', '/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(i) + '.fit')

				#copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(i) + '-g.fit', '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/g/' + str(i) + '-' + str(index) + '-g.fit')
				

				# Download from SDSS DAS the griz bands of data


				"""
				link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-g%d-%04d.fit.gz' % (int(line[44]), int(line[49]), int(line[50]), int(line[44]), int(line[50]), int(line[51]))
				urlretrieve(link, '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/g/' + str(k) + '-g.fit')

				link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-r%d-%04d.fit.gz' % (int(line[44]), int(line[49]), int(line[50]), int(line[44]), int(line[50]), int(line[51]))
				urlretrieve(link, '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/r/' + str(k) + '-r.fit')
				
				link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-i%d-%04d.fit.gz' % (int(line[44]), int(line[49]), int(line[50]), int(line[44]), int(line[50]), int(line[51]))
				urlretrieve(link, '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/i/' + str(k) + '-i.fit')
				
				link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-z%d-%04d.fit.gz' % (int(line[44]), int(line[49]), int(line[50]), int(line[44]), int(line[50]), int(line[51]))
				urlretrieve(link, '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/z/' + str(k) + '-z.fit')
				


				# Download the fpObj binary tables and fpAtlas files

				link = 'http://das.sdss.org/imaging/%d/%d/objcs/%d/fpAtlas-%06d-%d-%04d.fit' % (int(line[44]), int(line[49]), int(line[50]), int(line[44]), int(line[50]), int(line[51]))
				urlretrieve(link, '/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(k) + '_Obj.fit')

				link = 'http://das.sdss.org/imaging/%d/%d/objcs/%d/fpAtlas-%06d-%d-%04d.fit' % (int(line[44]), int(line[49]), int(line[50]), int(line[44]), int(line[50]), int(line[51]))
				urlretrieve(link, '/data/marvels/billzhu/MG II Atlas/0.37 - 0.55/' + str(k) + '_Atlas.fit')
				"""

				"""
				if (str(i) + '_SUB.fit') in refcutdirs:
					print(i)
					os.rename('/data/marvels/billzhu/Reference PSF Subtract 2/0.37 - 0.55/' + str(i) + '_SUB.fit', '/data/marvels/billzhu/Reference PSF Subtract 2/0.37 - 0.55/' + str(i) + '_' + str(index) + '_SUB.fit')

				"""
				# Rename the reference PSF final, rescaled images to include their MG II QSO index for easy access later


				count += 1
				usedref.append(i)

			except:
				continue

		if count == 4:
			return
		





def getmg(i):
	line = linecache.getline('Full Data.txt', i)
	refdata = line.split()
	z = float(refdata[3])
	refgmag = float(refdata[6])

	for j in range(len(mgtable)):
		if abs(mgtable['RA'][j] - float(refdata[1])) < 0.0001 and abs(mgtable['DEC'][j] - float(refdata[2])) < 0.0001:
			ismg[i] = True
			print(i)
			break

		
if __name__ == '__main__':
	#multiprocessing.set_start_method('spawn')

	#try:

	dirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/g/')
	mgtable = Table.read('QSObased_Trimmed_SDSS_DR7_107.fits', hdu=1)

	usedmg = []
	#rangelist = []
	usedref = []
	
	writer = open('Reference QSO NEW.txt', 'w')
	#mgdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/g/')
	#refdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/')
	#refcutdirs = os.listdir('/data/marvels/billzhu/Reference PSF Subtract 2/0.37 - 0.55')
	
	#print(refcutdirs)
	ismg = np.zeros(105784, dtype=bool)


	
	for i in range(1, 105784):
		line = linecache.getline('Full Data.txt', i)
		refdata = line.split()
		z = float(refdata[3])
		#refgmag = float(refdata[6])

		for j in range(len(mgtable)):
			if abs(mgtable['RA'][j] - float(refdata[1])) < 0.0001 and abs(mgtable['DEC'][j] - float(refdata[2])) < 0.0001:
				ismg[i] = True
				print(i)
				#print(ismg[i])
				break
	


	#pool = Pool(multiprocessing.cpu_count())
	#pool.map(getmg, np.arange(1, 105784))


	writer.write("%d\n" % (sum(ismg)))

	for d in dirs:
		index = d.split('_')[0]
		if index not in usedmg: 
			#rangelist.append(index)
			begin(index)
			usedmg.append(index)
	
	writer.close()


	# Not using multiprocessing due to shared data, will resolve later

	"""
	pool = Pool(os.cpu_count())
	pool.map(begin, rangelist)
	#except:
	#print("Error: Unable to process file")     
	"""
	
	
