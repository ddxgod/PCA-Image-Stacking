from astropy.io import fits
from astropy.table import Table
import numpy as np
np.set_printoptions(threshold=np.inf)
import string
import urllib
from urllib.request import urlretrieve
import os
import bz2
import multiprocessing
from multiprocessing import Pool
import linecache


#f = fits.open('final_catalog_full.fit')
#table2175 = f[1].data
#tableDR12 = Table.read('DR12Q.fits')
#table2175 = Table.read('final_catalog_full.fit')


def begin(i):
    for j in range(len(tableDR12)):
        if table2175['RA'][0][i] == tableDR12['RA'][j] and table2175['DEC'][0][i] == tableDR12['DEC'][j]:
            print("%f, %f" % (i, j))
            
            """
            link = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%d/%d/frame-g-%06d-%d-%04d.fits.bz2' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
            urlretrieve(link, '/data/marvels/billzhu/2175 Dataset/g/' + str(i) + '-' + str(j) + '-g.fit')
            zipfile = bz2.BZ2File('/data/marvels/billzhu/2175 Dataset/g/' + str(i) + '-' + str(j) + '-g.fit')
            data = zipfile.read()
            open('/data/marvels/billzhu/2175 Dataset/g/' + str(i) + '-' + str(j) + '-g.fit', 'wb').write(data)

            link = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%d/%d/frame-r-%06d-%d-%04d.fits.bz2' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
            urlretrieve(link, '/data/marvels/billzhu/2175 Dataset/r/' + str(i) + '-' + str(j) + '-r.fit')
            zipfile = bz2.BZ2File('/data/marvels/billzhu/2175 Dataset/r/' + str(i) + '-' + str(j) + '-r.fit')
            data = zipfile.read()
            open('/data/marvels/billzhu/2175 Dataset/r/' + str(i) + '-' + str(j) + '-r.fit', 'wb').write(data)

            link = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%d/%d/frame-i-%06d-%d-%04d.fits.bz2' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
            urlretrieve(link, '/data/marvels/billzhu/2175 Dataset/i/' + str(i) + '-' + str(j) + '-i.fit')
            zipfile = bz2.BZ2File('/data/marvels/billzhu/2175 Dataset/i/' + str(i) + '-' + str(j) + '-i.fit')
            data = zipfile.read()
            open('/data/marvels/billzhu/2175 Dataset/i/' + str(i) + '-' + str(j) + '-i.fit', 'wb').write(data)

            link = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%d/%d/frame-z-%06d-%d-%04d.fits.bz2' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
            urlretrieve(link, '/data/marvels/billzhu/2175 Dataset/z/' + str(i) + '-' + str(j) + '-z.fit')
            zipfile = bz2.BZ2File('/data/marvels/billzhu/2175 Dataset/z/' + str(i) + '-' + str(j) + '-z.fit')
            data = zipfile.read()
            open('/data/marvels/billzhu/2175 Dataset/z/' + str(i) + '-' + str(j) + '-z.fit', 'wb').write(data)
            """

            link = 'https://data.sdss.org/sas/dr12/boss/photoObj/301/%d/%d/photoObj-%06d-%d-%04d.fits' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
            urlretrieve(link, '/data/marvels/billzhu/2175 Obj/' + str(i) + '-' + str(j) + '.fit')
            
            return
            #except:
            #print("Unable to download")

    print(str(i) + 'No file found')





# Code that opens up a maximum of 8 processes for concurrent execution
    
if __name__ == '__main__':
    tableDR12 = Table.read('DR12Q.fits')
    table2175 = Table.read('final_catalog_full.fit')

    rangelist = np.arange(len(table2175['PLATE'][0]))
    #rangelist = []
    #rangelist.append(0)
    #print(len(table2175['PLATE'][0]))

    #try:
    print(os.cpu_count())
    pool = Pool(os.cpu_count())
    #print(rangelist)
    pool.map(begin, rangelist)
    #except:
    #print("Error: Unable to process file")

