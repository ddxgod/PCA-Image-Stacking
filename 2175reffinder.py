from astropy.io import fits
from astropy.table import Table
import numpy as np
import urllib
from urllib.request import urlretrieve
import os
import bz2
from shutil import copyfile
import string

tableDUST = Table.read('/data/marvels/billzhu/final_catalog_full.fit', hdu=1)
tableDR12 = Table.read('/data/marvels/billzhu/DR12Q.fits', hdu=1)
dirs = os.listdir('/data/marvels/billzhu/2175 Dataset/g/')
print(np.shape(tableDUST[0]))
print(np.shape(tableDR12))
used = []

writer = open('2175 Reference QSO.txt', 'w')

for i in range(len(dirs)):
    dust_index = int(dirs[i].split('-')[0])
    catalog_index = int(dirs[i].split('-')[1])
    count = 0
    #print(dust_index)
    for j in range(len(tableDR12)):
        if tableDR12['Z_MGII'][j] == -1 and abs(tableDUST[0]['Z'][dust_index] - tableDR12['Z_VI'][j]) < 0.1 and abs(tableDR12['PSFMAG'][catalog_index][0] - tableDR12['PSFMAG'][j][0]) < 0.5 and abs(tableDR12['PSFMAG'][catalog_index][1] - tableDR12['PSFMAG'][j][1]) < 0.5 and abs(tableDR12['PSFMAG'][catalog_index][2] - tableDR12['PSFMAG'][j][2]) < 0.5 and abs(tableDR12['PSFMAG'][catalog_index][3] - tableDR12['PSFMAG'][j][3]) < 0.5 and abs(tableDR12['PSFMAG'][catalog_index][4] - tableDR12['PSFMAG'][j][4]) < 0.5 and (j not in used):
            writer.write("%d %d %d\n" % (dust_index, catalog_index, j))
            used.append(j)
            count += 1

            
            print("%f, %f" % (dust_index, j))
            """
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


            link = 'https://data.sdss.org/sas/dr12/boss/photo/redux/301/%d/objcs/%d/fpObjc-%06d-%d-%04d.fit' % (tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['RUN_NUMBER'][j], tableDR12['COL_NUMBER'][j], tableDR12['FIELD_NUMBER'][j])
            urlretrieve(link, '/data/marvels/billzhu/2175 Reference Obj/' + str(dust_index) + '-' + str(j) + '.fit')
            """

        if count == 4:
            break


writer.close()
