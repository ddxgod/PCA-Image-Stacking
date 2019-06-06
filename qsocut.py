#!/c:/Anaconda3

from astropy.io import fits
from astropy.table import Table
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from astropy import wcs
import photutils
from photutils import centroid_2dg
from photutils import centroid_com
from photutils import centroid_1dg
from photutils import DAOStarFinder
from photutils.psf import FittableImageModel
import scipy as sp
from scipy import interpolate
import numpy as np
import string
import decimal
import matplotlib.pyplot as plt
import fileinput
import math
import os
import time
import linecache
import sewpy

import multiprocessing
from multiprocessing import Pool
import fileinput

#writer = open('Pixel Coordinates 50000')
#writer = open('MG II Pixel Coordinates.txt', 'a')
reader = open('Full Data.txt', 'r')



# Class the defines the point object, which specifies ra, dec, and image pixel coordinates

class point:
    ra = 0
    dec = 0
    x = 0
    y = 0

    def __init__(self, ra, dec, x, y):
        self.ra = ra
        self.dec = dec
        self.x = x
        self.y = y

    def calculate(other):
        return math.sqrt((self.x - other.x) ** 2 + (self.y - other.y) ** 2)

    

    
# Class that defines the source object

class source:
    def __init__(self, i):
        self.i = i

    
    def inbounds(x, y):
        if x > 0 and x < 2048 and y > 0 and y < 1489:
            return True
        return False



    
# Calculates the distance btween two given points

def distance(x, y, x1, y1):
    return math.sqrt((x - x1)**2 + (y - y1)**2)




# Method that checks if a source > 1 sigma is outside of a certain radius, and if so, masks it by putting mean value

def checkOutter(data, mean, std):
    count = 0
    for i in range(len(data)):
        for j in range(len(data[0])):
            if data[i][j] > mean + std and distance(i, j, 50, 50) > 25:
                data[i][j] = mean
    return data



# Method that checks if a sources > 3 sigma is inside of a certain radius, and if so, rejects the image

def checkInner(data, sources):
    count = 0
    for i in range(len(sources)):
        if distance(sources['xcentroid'][i], sources['ycentroid'][i], 50, 50) < 30:
            count += 1
        if count >= 2:
            return True

    return False




# Method that carries out the the splicing of the quasar cutout, reading data from the complete DR7 quasar catalog and from the header of the corresponding fits file that contains the quasar
# Calculates the centroid of the quasar and uses that as center of an 2x upsampled 2d cutout of the quasar

def begin(index):
    print(index)
    #mgi = int(index.split('-')[1])
    i = int(index.split('-')[0])
    #mgi = int(index.split('-')[1])
    color = index.split('-')[1]
    #print("%f, %f" % (i, mgi))
    #i = int(index)

    line = linecache.getline('Full Data.txt', i)
    #filename = 'Final Data Extract/' + str(i) + '.fit'
    #filename = 'MG II Dataset/' + str(i) + '.fit'


    #f1 = '/data/marvels/billzhu/2175 Dataset/' + str(i) + '-' + str(mgi) + '-g.fit'
    f1 = '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '.fit'
    #f1 = '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + str(mgi) + '-' + color + '.fit'


    hdulist = 0
    try:
        hdulist = fits.open(f1)
    except:
        print("FILE DOES NOT EXIST")
        return
    
    prihdr = hdulist[0].header

    #raDegColPix = 0
    #raDegRowPix = 0
    #decDegColPix = 0
    #decDegRowPix = 0


    # If statements for finding out ra and dec pixel specifications

    """
    if 'RA' in prihdr.comments['CD1_1']:
        raDegColPix = prihdr['CD1_1']
        raDegRowPix = prihdr['CD1_2']
        decDegColPix = prihdr['CD2_1']
        decDegRowPix = prihdr['CD2_2']
    else:
        decDegColPix = prihdr['CD1_1']
        decDegRowPix = prihdr['CD1_2']
        raDegColPix = prihdr['CD2_1']
        raDegRowPix = prihdr['CD2_2']

    ra1 = 0
    dec1 = 0
    if 'RA' in prihdr.comments['CRVAL1']:
        ra1 = prihdr['CRVAL1']
        dec1 = prihdr['CRVAL2']
    else:
        ra1 = prihdr['CRVAL2']
        dec1 = prihdr['CRVAL1']

    x1 = 0
    y1 = 0
    if 'X' in prihdr.comments['CRPIX1']:
        x1 = prihdr['CRPIX1']
        y1 = prihdr['CRPIX2']
    else:
        x1 = prihdr['CRPIX2']
        y1 = prihdr['CRPIX1']

    refpix = point(ra1, dec1, x1, y1)
    #print("%f, %f, %f, %f" % (x1, y1, ra1, dec1))
    """

    # Calculates the coordinates given RA and DEC and the data of a reference pixel

    data = line.split()
    qRA = float(data[1])
    qDEC = float(data[2])

    #bigtable = Table.read('final_catalog_full.fit')
    #qRA = bigtable['RA'][0][i]
    #qDEC = bigtable['DEC'][0][i]

    #print("%f, %f, %f" % (i, qRA, qDEC))
    #print("%f, %f, %f, %f" % (raDegColPix, raDegRowPix, decDegColPix, decDegRowPix))


    
    # Package function yields much higher precision than hard-coded function

    print("%f, %f" % (qRA, qDEC))
    wcstransform = wcs.WCS(prihdr)
    x_test, y_test = wcstransform.wcs_world2pix(qRA, qDEC, 0)
    print("%f, %f" % (x_test, y_test))

    


    

    """
    p = qRA - refpix.ra + raDegColPix * refpix.x + raDegRowPix * refpix.y
    q = qDEC - refpix.dec + decDegColPix * refpix.x + decDegRowPix * refpix.y
    qX = int((p * decDegRowPix - q * raDegRowPix) / (raDegColPix * decDegRowPix - raDegRowPix * decDegColPix))
    qY = int((q * raDegColPix - p * decDegColPix) / (raDegColPix * decDegRowPix - raDegRowPix * decDegColPix))


    # Checks for RA and DEC boundary issues
    
    if source.inbounds(qX, qY) == False:
        #print(str(i))
        if qRA < refpix.ra and refpix.ra - qRA > 358:
            qRA += 360
        else:
            if qRA > refpix.ra and qRA - refpix.ra > 358:
                refpix.ra += 360
        if qDEC < refpix.dec and refpix.dec - qDEC > 178:
            qDEC += 90                      
        else:
             if qDEC > refpix.dec and qDEC - refpix.dec > 178:
                refpix.dec += 90
        p = qRA - refpix.ra + raDegColPix * refpix.x + raDegRowPix * refpix.y
        q = qDEC - refpix.dec + decDegColPix * refpix.x + decDegRowPix * refpix.y
        qX = int((p * decDegRowPix - q * raDegRowPix) / (raDegColPix * decDegRowPix - raDegRowPix * decDegColPix))
        qY = int((q * raDegColPix - p * decDegColPix) / (raDegColPix * decDegRowPix - raDegRowPix * decDegColPix))
    
       
    print("%f,  %f" % (qX, qY))
    """

    """
    scidata = hdulist[0].data.astype(float)
    obj_table = Table.open('/data/marvels/billzhu/MG II Obj/0.37 - 0.55/' + str(i) + '.fit')[0]

    pointer = 0
    if color == 'g':
        pointer = 1
    if color == 'r':
        pointer = 2
    if color == 'i':
        pointer = 3
    if color == 'z':
        pointer = 4
    if color = 'u':
        pointer = 0

    quasar = 0
    for j in range(len(obj_table)):
        if abs(obj_table['rowc'][j] - y_test) < 5 and abs(obj_table['colc'][j] - x_test) < 5 and obj_table['nchild'][j] > 0:
            quasar = obj_table[j]
            break

    children = []
    for j in range(len(obj_table)):
        if obj_table['parent'][j] == quasar['id']:
            children.append(obj_table[j])
    """
    
    try:
        qX = int(x_test)
        qY = int(y_test)
    except:
        return

    half = 15
    yl = qY - half
    yu = qY + half
    xl = qX - half
    xu = qX + half

    if yl < 0:
        yl = 0
    if yu > len(scidata):
        yu = len(scidata)
    if xl < 0:
        xl = 0
    if xu > len(scidata[0]):
        xu = len(scidata[0])
    image = scidata[yl : yu, xl : xu]

    bkg_sigma = mad_std(scidata)
    daofind = DAOStarFinder(fwhm = 2., threshold = 5.*bkg_sigma)
    sources = daofind(image)
    sources['xcentroid'] += xl
    sources['ycentroid'] += yl


    print(sources)

    """
    deg_diff = 1000000
    quasar = 0
    ra11 = 0
    dec11 = 0
    for j in range(len(sources)):
        dra_maybe = abs(sources['xcentroid'][j] - refpix.x) * raDegColPix + abs(sources['ycentroid'][j] - refpix.y) * raDegRowPix
        ddec_maybe = abs(sources['xcentroid'][j] - refpix.x) * decDegColPix + abs(sources['ycentroid'][j] - refpix.y) * decDegRowPix
        dra_real = abs(qRA - refpix.ra)
        ddec_real = abs(qDEC - refpix.dec)

        #print(math.sqrt((dra_maybe - dra_real)**2 + (ddec_maybe - ddec_real)**2))
        
        if math.sqrt((dra_maybe - dra_real)**2 + (ddec_maybe - ddec_real)**2) < deg_diff:
            deg_diff = math.sqrt((dra_maybe - dra_real)**2 + (ddec_maybe - ddec_real)**2)
            quasar = sources[j]
    """

    diff = 1000000
    quasar = 0
    for j in range(len(sources)):
        dist1 = distance(sources['xcentroid'][j], sources['ycentroid'][j], qX, qY)
        if dist1 < diff:
            diff = dist1
            quasar = sources[j]

    if quasar == 0:
        #print(deg_diff)
        print("ERROR NO QUASAR")
        return

    print(quasar)
    
    
                
    # Write the right coordinates to fits file header for access later for PSF subtraction
    
    hdulist[0].header.append(('XCOORD', int(quasar['xcentroid']), 'x coordinate of quasar in image'), end = True)
    hdulist[0].header.append(('YCOORD', int(quasar['ycentroid']), 'y coordinate of quasar in image'), end = True)

    
    # Code for shifting the quasar to the actual centroid
    
    chunk_size = 50
    if source.inbounds(quasar['xcentroid'] + chunk_size, quasar['ycentroid'] + chunk_size) and source.inbounds(quasar['xcentroid'] - chunk_size, quasar['ycentroid'] - chunk_size):
        #print("SUCCESS")

        
        """
        #try:
        preshift = scidata[qYmax - int(chunk_size) : qYmax + int(chunk_size) + 1, qXmax - int(chunk_size) : qXmax + int(chunk_size) + 1]

        mean, median, std = sigma_clipped_stats(preshift, sigma=3.0, iters=5)

        qXc, qYc = centroid_1dg(preshift, mask = None)
        #print("%f  %f" % (qXc, qYc))
        qXc += qXmax - chunk_size
        qYc += qYmax - chunk_size
        #print("%f  %f" % (qXc, qYc))
        """        

        qXc = quasar['xcentroid']
        qYc = quasar['ycentroid']

        
        preshift = scidata[int(qYc) - chunk_size - 5: int(qYc) + chunk_size + 6, int(qXc) - chunk_size - 5 : int(qXc) + chunk_size + 6]

        
        """
        plt.imshow(preshift, origin='lower', interpolation='nearest', cmap='viridis')
        marker = '+'
        ms, mew = 30, 2.
        #plt.plot(qXc - qXma + chunk_size, qYc - qYmax + chunk_size, color='#1f77b4', marker=marker, ms=ms, mew=mew)
        plt.show()
        plt.pause(3)
        """
        

        xr = np.arange(int(qXc) - chunk_size - 5, int(qXc) + chunk_size + 6, 1)
        yr = np.arange(int(qYc) - chunk_size - 5, int(qYc) + chunk_size + 6, 1)


        # Shifts the data to center around the centroid using 2d interpolation

        try:
            if(len(xr) == len(preshift[0])):
                #print("SUCCESS 2.0")
                shifted = []

                spline = interpolate.interp2d(xr, yr, preshift)

                xrf = np.arange(qXc - chunk_size, qXc + chunk_size + 1, 1)
                yrf = np.arange(qYc - chunk_size, qYc + chunk_size + 1, 1)

                if len(xrf) > 101:
                    xrf = xrf[:-1].copy()
                if len(yrf) > 101:
                    yrf = yrf[:-1].copy()

                shifted = spline(xrf, yrf)

                """
                mean, median, stddev = sigma_clipped_stats(preshift, sigma=3.0, iters=5)
                preshift -= median
                check_sources = daofind.find_stars(preshift)

                if len(check_sources) > 1:
                    print(len(check_sources))
                    return
                """

                # If the source has a weird shape i.e. due to gravitational lensing, then check if the maximum pixel is within 3 pixels of the center to ensure consistency
                daofind = DAOStarFinder(fwhm = 2, threshold=3.0*bkg_sigma)
                sources = daofind.find_stars(shifted)
                #sources['xcentroid'] += xl
                #sources['ycentroid'] += yl
                cont = checkInner(shifted, sources)

                if cont == True:
                    return
                mean1, median1, std1 = sigma_clipped_stats(shifted, sigma=3.0, iters=5)
                #print("%f, %f" % (mean1, std1))
                shifted = checkOutter(shifted, mean1, std1)
                


                """
                max_coords = np.unravel_index(shifted.argmax(), shifted.shape)
                max_coords = list(max_coords)
                print(max_coords)
                for k in range(len(max_coords)//2):
                    yt = max_coords[2 * k]
                    xt = max_coords[2 * k + 1]
                    print("%f, %f" % (xt, yt))
                    if distance(xt, yt, 40, 40) > 3:
                        print('MAX TOO FAR')
                        return
                """

                print('NO FAIL YET')
                #fits.writeto('/data/marvels/billzhu/2175 Quasar Cut/' + str(i) + '-' + str(mgi) + '_DUST.fit', shifted, hdulist[0].header, clobber = True)
                #fits.writeto('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + str(i) + '-' + str(mgi) + '_REF.fit', shifted, hdulist[0].header, clobber = True)
                #fits.writeto('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit', shifted, hdulist[0].header, clobber = True)
                fits.writeto('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + str(mgi) + '-' + color + '_REF.fit', shifted, hdulist[0].header, clobber = True)

        except:
            print(False)
            return

        return




# Main method that can be wired for multiprocessing purposes using Pool

if __name__ == '__main__':
    multiprocessing.set_start_method('spawn')

    """
    for j in range(1, 11):
        #print(j)


        try:
            begin(j)
        except:

            # Write new line to keep consistency in retrieving data based on line number

            writer.write('\n')
            print("Error: Unable to process")

    """

    exist = []
    #dirs = os.listdir('/data/marvels/billzhu/2175 Dataset/')
    #dirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/')

    #gdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/g/')
    #rdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/r/')
    #idirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/i/')
    #zdirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/z/')
    #udirs = os.listdir('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/u/')

    gdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/g/')

    #check_dirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/')
    rangelist = []
    #rangelist.append()

    #print(rangelist)
    for d in gdirs:
        index = d.split('.')[0]
        rangelist.append(index)
        
        #if d.split('-')[0] + '_MG.fit' not in check_dirs:
        #rangelist.append(index + '-r')
        #rangelist.append(index + '-i')
        #rangelist.append(index + '-z')
        #rangelist.append(index + '-u')

    """
    for d in rdirs:
        index = d.split('.')[0]
        #if d.split('-')[0] + '_MG.fit' not in check_dirs:
        rangelist.append(index)

    for d in idirs:
        index = d.split('.')[0]
        #if d.split('-')[0] + '_MG.fit' not in check_dirs:
        rangelist.append(index)

    for d in zdirs:
        index = d.split('.')[0]
        #if d.split('-')[0] + '_MG.fit' not in check_dirs:
        rangelist.append(index)

    for d in udirs:
        index = d.split('.')[0]
        #if d.split('-')[0] + '_MG.fit' not in check_dirs:
        rangelist.append(index)
    """

    #try:
    #rlist = []
    #rlist.append('9807-g')
    pool = Pool(os.cpu_count())
    pool.map(begin, rangelist)
    #except:
    #print("Error: Unable to process file")

    print(True)
        
