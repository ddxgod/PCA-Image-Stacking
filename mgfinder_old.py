# cross-checks Dr. Ben Guangtun Zhu's MG II catalog with the DR7 catalog to retrieve the correct files




from astropy.table import Table
from astropy.io import fits
from shutil import copyfile
import numpy as np
np.set_printoptions(threshold=np.inf)
import scipy as sp
import math
import os
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
        if mgtable['ZABS'][j] >= 0.37 and mgtable['ZABS'][j] < 0.55:
            print(mgtable['ZABS'][j])
            #copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(mgtable['INDEX_QSO'][j] + 1) + '-g.fit', '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + str(mgtable['INDEX_QSO'][j] + 1) + '-g.fit')


            for k in range(1, 105784):
                temp = linecache.getline('Full Data.txt', k).split()
                #print(k)
                if abs(float(temp[1]) - mgtable['RA'][j]) < 0.0001 and abs(float(temp[2]) - mgtable['DEC'][j]) < 0.0001:
                    print(k)
                    copyfile('/data/marvels/jerryxu/dr7/raw_catalog/' + str(k) + '.fit', '/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(k) + '.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-g.fit', '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/g/' + str(k) + '-g.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-r.fit', '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/r/' + str(k) + '-r.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-i.fit', '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/i/' + str(k) + '-i.fit')
                    #copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(i) + '-u.fit', '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + str(i) + '-u.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-z.fit', '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/z/' + str(k) + '-z.fit')
                    break

                
        """
        if mgtable['ZABS'][j][0] >= 0.55 and mgtable['ZABS'][j][0] < 0.76:

            for k in range(1, 105784):
                temp = linecache.getline('Full Data.txt', k).split()
                #print(k)
                if abs(float(temp[1]) - mgtable['RA'][j]) < 0.001:
                    print(k)
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-g.fit', '/data/marvels/billzhu/MG II Dataset/0.55 - 0.76/g/' + str(k) + '-g.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-r.fit', '/data/marvels/billzhu/MG II Dataset/0.55 - 0.76/r/' + str(k) + '-r.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-i.fit', '/data/marvels/billzhu/MG II Dataset/0.55 - 0.76/i/' + str(k) + '-i.fit')
                    #copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(i) + '-u.fit', '/data/marvels/billzhu/MG II Dataset/0.55 - 0.76/' + str(i) + '-u.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-z.fit', '/data/marvels/billzhu/MG II Dataset/0.55 - 0.76/z/' + str(k) + '-z.fit')
                    copyfile('/data/marvels/jerryxu/dr7/raw_catalog/' + str(k) + '.fit', '/data/marvels/billzhu/MG II Obj/0.55 - 0.76/' + str(k) + '.fit')
                    break

        if mgtable['ZABS'][j][0] >= 0.76 and mgtable['ZABS'][j][0] < 1.00:
            for k in range(1, 105784):
                temp = linecache.getline('Full Data.txt', k).split()
                #print(k)
                #print(temp)
                #print(mgtable['RA'][j])
                if abs(float(temp[1]) - mgtable['RA'][j]) < 0.001:
                    print(k)
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-g.fit', '/data/marvels/billzhu/MG II Dataset/0.76 - 1.00/g/' + str(k) + '-g.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-r.fit', '/data/marvels/billzhu/MG II Dataset/0.76 - 1.00/r/' + str(k) + '-r.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-i.fit', '/data/marvels/billzhu/MG II Dataset/0.76 - 1.00/i/' + str(k) + '-i.fit')
                    #copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(i) + '-u.fit', '/data/marvels/billzhu/MG II Dataset/0.76 - 1.00/' + str(i) + '-u.fit')
                    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(k) + '-z.fit', '/data/marvels/billzhu/MG II Dataset/0.76 - 1.00/z/' + str(k) + '-z.fit')
                    copyfile('/data/marvels/jerryxu/dr7/raw_catalog/' + str(k) + '.fit', '/data/marvels/billzhu/MG II Obj/0.76 - 1.00/' + str(k) + '.fit')
                    break
        """
        
    except:
        return




if __name__ == '__main__':
    #multiprocessing.set_start_method('spawn')
    #f = fits.open('QSObased_Trimmed_SDSS_DR7_107.fits')
    mgtable = Table.read('Trimmed_SDSS_DR7_107.fits', hdu=1)
    reader = open('Full Data.txt', 'r')
    count = 0
    
    
    #try:
    pool = Pool(os.cpu_count())
    print(len(mgtable))
    pool.map(begin, np.arange(len(mgtable)))

    print(len(mgtable))
    for j in range(len(mgtable)):
        if mgtable['ZABS'][j] >= 0.37 and mgtable['ZABS'][j] < 0.55:
            count += 1

    #print(count)
    #except:
    #print("Error: Unable to process file")

