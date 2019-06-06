from astropy.table import Table
from astropy.io import fits
from shutil import copyfile
import numpy as np
np.set_printoptions(threshold=np.inf)
import scipy as sp
import math
import os
import urllib
from urllib.request import urlretrieve
import multiprocessing
from multiprocessing import Pool
import linecache



def distance(ra, dec, xra, xdec):
	return math.sqrt((ra - xra)**2 + (dec - xdec)**2)



# Run through all the QSOs that have Mg II Absorbers
# USE THE FRICKIN QSO CATALOG
# Since there may be multiple absorbers, use the lowest redshift one as marker

# Instead of interating over all absorbers, since it's already sorted by redshift, choose the lowest one i.e. [0]
# This ensures that there will be no low redshift absorbers in each bin, only possibly high redshift ones

def begin(j):
	
	#print(j)
	try:
		#print(mgtable['ZABS'][j])
		if mgtable['NABS'][j] == 1 and mgtable['ZABS'][j][0] >= 0.37 and mgtable['ZABS'][j][0] < 0.55 and mgtable['REW_MGII_2796'][j][0] > 0.8:
			#print(mgtable['ZABS'][j])
			#copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(mgtable['INDEX_QSO'][j] + 1) + '-g.fit', '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + str(mgtable['INDEX_QSO'][j] + 1) + '-g.fit')


			for k in range(1, 105784):
				data = linecache.getline('Full Data.txt', k).split()
				#print(k)
				if abs(float(data[1]) - mgtable['RA'][j]) < 0.0001 and abs(float(data[2]) - mgtable['DEC'][j]) < 0.0001:
					print(mgtable['ZABS'][j][0])
					#print(data)

					
					# Download from SDSS DAS the griz bands of data

					"""
					link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-g%d-%04d.fit.gz' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
					urlretrieve(link, '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/g/' + str(k) + '-g.fit')

					link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-r%d-%04d.fit.gz' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
					urlretrieve(link, '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/r/' + str(k) + '-r.fit')
					
					link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-i%d-%04d.fit.gz' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
					urlretrieve(link, '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/i/' + str(k) + '-i.fit')
					
					link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-z%d-%04d.fit.gz' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
					urlretrieve(link, '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/z/' + str(k) + '-z.fit')
					"""


					# Download the fpObj binary tables and fpAtlas files
					"""
					link = 'http://das.sdss.org/imaging/%d/%d/objcs/%d/fpObjc-%06d-%d-%04d.fit' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
					urlretrieve(link, '/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(k) + '_Obj.fit')

					
					link = 'http://das.sdss.org/imaging/%d/%d/objcs/%d/fpAtlas-%06d-%d-%04d.fit' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
					urlretrieve(link, '/data/marvels/billzhu/MG II Atlas/0.37 - 0.55/' + str(k) + '_Atlas.fit')
					break
					
					link = 'http://das.sdss.org/imaging/%d/%d/calibChunks/%d/tsField-%06d-%d-%d-%04d.fit' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[49]), int(data[51]))
					urlretrieve(link, '/data/marvels/billzhu/MG II tsField/0.37 - 0.55/' + str(k) + '_tsField.fit')
					"""

					link = 'http://das.sdss.org/imaging/%d/%d/objcs/%d/psField-%06d-%d-%04d.fit' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
					urlretrieve(link, '/data/marvels/billzhu/MG II psField/0.37 - 0.55/' + str(k) + '_psField.fit')
					break
					


		
	except:
		return




if __name__ == '__main__':
	#multiprocessing.set_start_method('spawn')
	#f = fits.open('QSObased_Trimmed_SDSS_DR7_107.fits')
	mgtable = Table.read('QSObased_Trimmed_SDSS_DR7_107.fits', hdu=1)
	reader = open('Full Data.txt', 'r')
	count = 0
	
	
	#try:
	pool = Pool(os.cpu_count())
	print(len(mgtable))
	pool.map(begin, np.arange(len(mgtable)))

	print(len(mgtable))
	for j in range(len(mgtable)):
		if mgtable['ZABS'][j][0] >= 0.37 and mgtable['ZABS'][j][0] < 0.55:
			count += 1

	print(count)
	#except:
	#print("Error: Unable to process file")
