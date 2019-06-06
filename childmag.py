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
sys.setrecursionlimit(15000000)


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



obj_id = objid_extract(1237658425174393019, False)['id']
print(obj_id)



"""
dustdirs = os.listdir('2175 Dataset/g/')
dust_table = Table.read('final_catalog_full.fit')
dr12_table = Table.read('DR12Q.fits')
numchild = 0

for j in range(len(dustdirs)):
    print(j)
    i = int(dustdirs[j].split('-')[0])
    mgi = int(dustdirs[j].split('-')[1])

    obj_table = Table.read('/data/marvels/billzhu/2175 Obj/' + str(i) + '-' + str(mgi) + '.fit', hdu=1) 
    obj_id = objid_extract(int(dr12_table[mgi]['OBJ_ID']), False)['id']
    
    if obj_table['PARENT'][obj_id - 1] != -1 or obj_table['NCHILD'][obj_id - 1] > 0:
        print(obj_table['PSFMAG'][obj_id - 1][0])
        numchild += 1

    if obj_table['PARENT'][obj_id - 1] != -1:
        print(obj_table['PSFMAG'][obj_table['PARENT'][obj_id - 1]][0])

    if obj_table['NCHILD'][obj_id - 1] > 0:
        counter = obj_table['NCHILD'][obj_id - 1]

        while counter > 0:
            print(obj_table['PSFMAG'][obj_id - 1 + counter][0])
            counter -= 1

print(numchild)
"""