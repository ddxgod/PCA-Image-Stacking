from astropy.table import Table
import numpy as np
import os
import linecache
import urllib
from urllib.request import urlretrieve
import os
import bz2
import multiprocessing
from multiprocessing import Pool


def begin(i):
    line = linecache.getline('2175 Reference QSO.txt', i)
    line_data = line.split()
    dust_index = int(line_data[0])
    catalog_index = int(line_data[1])
    j = int(line_data[2])
    print("%d, %d, %d" % (dust_index, catalog_index, j))

    if str(dust_index) + '-' + str(j) + '_Atlas.fit' in check_dirs:
        return

    
    link = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%d/%d/frame-g-%06d-%d-%04d.fits.bz2' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
    urlretrieve(link, '/data/marvels/billzhu/2175 Reference Dataset/g/' + str(dust_index) + '-' + str(j) + '-g.fit')
    zipfile = bz2.BZ2File('/data/marvels/billzhu/2175 Reference Dataset/g/' + str(dust_index) + '-' + str(j) + '-g.fit')
    data = zipfile.read()
    open('/data/marvels/billzhu/2175 Reference Dataset/g/' + str(dust_index) + '-' + str(j) + '-g.fit', 'wb').write(data)

    link = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%d/%d/frame-r-%06d-%d-%04d.fits.bz2' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
    urlretrieve(link, '/data/marvels/billzhu/2175 Reference Dataset/r/' + str(dust_index) + '-' + str(j) + '-r.fit')
    zipfile = bz2.BZ2File('/data/marvels/billzhu/2175 Reference Dataset/r/' + str(dust_index) + '-' + str(j) + '-r.fit')
    data = zipfile.read()
    open('/data/marvels/billzhu/2175 Reference Dataset/r/' + str(dust_index) + '-' + str(j) + '-r.fit', 'wb').write(data)

    link = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%d/%d/frame-i-%06d-%d-%04d.fits.bz2' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
    urlretrieve(link, '/data/marvels/billzhu/2175 Reference Dataset/i/' + str(dust_index) + '-' + str(j) + '-i.fit')
    zipfile = bz2.BZ2File('/data/marvels/billzhu/2175 Reference Dataset/i/' + str(dust_index) + '-' + str(j) + '-i.fit')
    data = zipfile.read()
    open('/data/marvels/billzhu/2175 Reference Dataset/i/' + str(dust_index) + '-' + str(j) + '-i.fit', 'wb').write(data)

    link = 'https://dr12.sdss.org/sas/dr12/boss/photoObj/frames/301/%d/%d/frame-z-%06d-%d-%04d.fits.bz2' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
    urlretrieve(link, '/data/marvels/billzhu/2175 Reference Dataset/z/' + str(dust_index) + '-' + str(j) + '-z.fit')
    zipfile = bz2.BZ2File('/data/marvels/billzhu/2175 Reference Dataset/z/' + str(dust_index) + '-' + str(j) + '-z.fit')
    data = zipfile.read()
    open('/data/marvels/billzhu/2175 Reference Dataset/z/' + str(dust_index) + '-' + str(j) + '-z.fit', 'wb').write(data)
    

    link = 'https://data.sdss.org/sas/dr12/boss/photoObj/301/%d/%d/photoObj-%06d-%d-%04d.fits' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
    urlretrieve(link, '/data/marvels/billzhu/2175 Reference Obj/' + str(dust_index) + '-' + str(j) + '.fit')


    link = 'https://data.sdss.org/sas/dr12/boss/photo/redux/301/%d/objcs/%d/fpAtlas-%06d-%d-%04d.fit' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
    urlretrieve(link, '/data/marvels/billzhu/2175 Reference Atlas/' + str(dust_index) + '-' + str(j) + '_Atlas.fit')




if __name__ == '__main__':
    tableDUST = Table.read('/data/marvels/billzhu/final_catalog_full.fit', hdu=1)
    tableDR12 = Table.read('/data/marvels/billzhu/DR12Q.fits', hdu=1)
    #multiprocessing.Array([tableDUST, tableDR12], 2)

    rangelist = np.arange(1, sum(1 for line in open('2175 Reference QSO.txt')) + 1)
    check_dirs = os.listdir('/data/marvels/billzhu/2175 Reference Atlas/')


    print(rangelist)
    pool = Pool(os.cpu_count())
    print('reached')
    pool.map(begin, rangelist)
