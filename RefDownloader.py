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




# Download the griz band field image files, fpObj binary table files, and fpAtlas image poststamps

def begin(j):
    try:
        
        line_data = linecache.getline('Reference QSO NEW.txt', j)
        index = int(line_data.split()[0])
        mgi = int(line_data.split()[1])
        data = linecache.getline('Full Data.txt', index).split()
        print(index)


        # Download from SDSS DAS the griz bands of data

        
        link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-g%d-%04d.fit.gz' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
        urlretrieve(link, '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/g/' + str(index) + '-' + str(mgi) + '-g.fit')

        link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-r%d-%04d.fit.gz' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
        urlretrieve(link, '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/r/' + str(index) + '-' + str(mgi) + '-r.fit')
        
        link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-i%d-%04d.fit.gz' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
        urlretrieve(link, '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/i/' + str(index) + '-' + str(mgi) + '-i.fit')
        
        link = 'http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%06d-z%d-%04d.fit.gz' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
        urlretrieve(link, '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/z/' + str(index) + '-' + str(mgi) + '-z.fit')
        


        # Download the fpObj binary tables and fpAtlas files
        
        link = 'http://das.sdss.org/imaging/%d/%d/objcs/%d/fpObjc-%06d-%d-%04d.fit' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
        urlretrieve(link, '/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(index) + '_Obj.fit')

        
        link = 'http://das.sdss.org/imaging/%d/%d/objcs/%d/psField-%06d-%d-%04d.fit' % (int(data[44]), int(data[49]), int(data[50]), int(data[44]), int(data[50]), int(data[51]))
        urlretrieve(link, '/data/marvels/billzhu/Reference psField/0.37 - 0.55/' + str(index) + '_psField.fit') 
                   
    except:
        return




if __name__ == '__main__':
    #multiprocessing.set_start_method('spawn')
    #f = fits.open('QSObased_Trimmed_SDSS_DR7_107.fits')
    #mgtable = Table.read('Trimmed_SDSS_DR7_107.fits', hdu=1)
    reader = open('Full Data.txt', 'r')
    count = 0
    
    num_lines = sum(1 for line in open('Reference QSO NEW.txt', 'r'))
    rangelist = np.arange(1, num_lines + 1)
    
    #try:
    pool = Pool(os.cpu_count())
    pool.map(begin, rangelist)


    print(num_lines)
    #except:
    #print("Error: Unable to process file")
