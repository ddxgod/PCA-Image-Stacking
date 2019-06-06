from astropy.io import fits
from astropy.table import Column
import scipy as sp
import scipy.optimize as opt
from scipy import interpolate
import numpy as np
np.set_printoptions(threshold=np.inf)
from numpy import *
import functools
import string
import decimal
import matplotlib.pyplot as plt
import fileinput
import math
import os
import photutils
from photutils import DAOStarFinder
from photutils import IRAFStarFinder
from photutils import subtract_psf
from photutils import centroid_2dg
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from functools import cmp_to_key

import multiprocessing
from multiprocessing import Pool
import linecache

# Checks if the given coordinates are within the bounds of the image

def inbounds(x, y):
    if x > 0 and x < 2048 and y > 0 and y < 1489:
        return True
    return False



# Calculates the distance btween two given points

def distance(x, y, x1, y1):
    return math.sqrt((x - x1)**2 + (y - y1)**2)



# Method that checks if a particular residue is noise or not

def checkNoise(x, y, qX, qY, data):
    halo = []
    for i in range(42):
        for j in range(42):
            if abs(distance(i, j, qX, qY) - distance(x, y, qX, qY)) <= 2 and distance(i, j, x, y) >= 2:
                halo.append(data[j, i])
    mean, median, std = sigma_clipped_stats(halo, sigma=3.0, iters=5)

    if data[y, x] > mean + 3 * std:
        return True
                

    
    
# Method that calculates the total photon count within 3 sigma of the quasar centroid

def photonCount(xc, yc, sigma, data):
    count = 0
    for i in range(len(data)):
        for j in range(len(data)):
            if distance(i, j, xc, yc) <= 5 * sigma and data[i][j] > 0:
                count += data[i][j]
    return count



# Method that checks if a bright source is outside of a certain radius

def checkOutter(data, mean, stddev):
    for i in range(len(data)):
        for j in range(len(data[0])):
            if data[i][j] > mean + stddev and distance(i, j, 50, 50) > 25:
                data[i][j] = mean

    return data

    

# Class that assists in choosing the optimal PSF sources by storing index and Z value

class tuplet:
    def __init__(self, i, z):
        self.i = i
        self.z = z
        
    def getSelf(z):
        return self

    def getIndex():
        return i

    def getZ():
        return z


#reader = open('Pixel Coordinates 50000.txt', 'r')
#reader = open('MG II Pixel Coordinates.txt', 'r')



# Main method for executing PSF Subtraction
# Only searches in a 600x600 box around the quasar for a PSF source
# If none are found then the file is passed on and not used

def begin(index):
    chunk_size = 50
    i = int(index.split('-')[0])
    color = index.split('-')[1]
    
    #line = linecache.getline('Full Data.txt', i)
    #num, x, y = linecache.getline('Pixel Coordinates 50000.txt', i).split()
    #x = int(x)
    #y = int(y)
    #print("%f   %f" % (x, y))
    

    try:
        #filename = 'Test Data Extract/' + str(i) + '.fit'
        filename = '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '.fit'

        #filename = 'Reference Dataset/' + str(i) + '_REF.fit'
        hdulist = fits.open(filename)

        qlist = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')

        #qlist = fits.open('Reference Quasar Cut/' + str(i) + '_REF.fit')
        x = qlist[0].header['XCOORD']
        y = qlist[0].header['YCOORD']
    except:
        print("No coordinates")
        return


    
    #half = 500
    scidata = hdulist[0].data.astype(float)
    mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
    change_top = False
    
    if x + chunk_size > 2048:
        change_top = True
        filler = np.array([float(median)] * len(scidata))
        for j in range(10):
            scidata = np.insert(scidata, len(scidata[0]), filler, 1)

    if x - chunk_size < 0:
        x += chunk_size
        filler = np.array([float(median)] * len(scidata))
        for j in range(10):
            scidata = np.insert(scidata, 0, filler, 1)

    if y + chunk_size > 1489:
        change_top = True
        filler = np.array([float(median)] * len(scidata[0]))
        for j in range(10):
            scidata = np.insert(scidata, len(scidata), filler, 0)

    if y - chunk_size < 0:
        y += chunk_size
        filler = np.array([float(median)] * len(scidata[0]))
        for j in range(10):
            scidata = np.insert(scidata, 0, filler, 0)

            
    scidata -= median    
    psfindex = -1
    quasar = 0
    bkg_sigma = mad_std(scidata)
    daofind = DAOStarFinder(fwhm = 2., threshold = 5.*bkg_sigma)
    sources = daofind(scidata[int(y - 10) : int(y + 10), int(x - 10) : int(x + 10)])

    
    # Update coordinates of the sources
    
    sources['xcentroid'] += x - 10
    sources['ycentroid'] += y - 10

    
    # Create new column that contains the FWHM of each source for comparison later
    
    FWHM = np.empty([len(sources)])
    column = Column(FWHM, name = 'FWHM')
    sources.add_column(column)


    
    # Find the quasar and calculate its FWHM

    qsigma = 0
    
    for j in range(len(sources)):
        if abs(round(sources['xcentroid'][j]) - x) < 2 and abs(round(sources['ycentroid'][j]) - y) < 2:
            quasar = sources[j]
            width = int(np.sqrt(sources['npix'][j]))
            #print("%d   %d   %f   %f" % (j, width, sources['xcentroid'][j], sources['ycentroid'][j]))
            data = scidata[int(sources['ycentroid'][j] - width/2) : int(sources['ycentroid'][j] + width/2) + 1, int(sources['xcentroid'][j] - width/2) : int(sources['xcentroid'][j] + width/2) + 1]


            """
            plt.imshow(data, origin='lower', interpolation='nearest', cmap='viridis')
            plt.show()
            plt.pause(3)
            """
            
            gauss = 0
            
            if(np.ma.count(data) >= 7):
                gauss = photutils.fit_2dgaussian(data, mask = None)
            
            fwhm = 0
            if gauss != 0:
                fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(gauss.x_stddev)
                quasar['FWHM'] = fwhm
                #qsigma = np.sqrt(gauss.x_stddev**2 + gauss.y_stddev**2)
                #print(quasar['FWHM'])
                break

            
    ztot = 10000
    print(quasar)

    
    # If no quasar is found, the field image is deemed corrupt and not used
    
    if quasar == 0:
        return


    # Define cutout image limits, adjusted to field image boundaries as necessary i.e. < 0 or > max x/y values
    """
    yl = y - half
    yu = y + half
    xl = x - half
    xu = x + half

    if yl < 0:
        yl = 0
    if yu > len(scidata):
        yu = len(scidata)
    if xl < 0:
        xl = 0
    if xu > len(scidata[0]):
        xu = len(scidata[0])
    image = scidata[yl : yu, xl : xu]
    """
    bkg_sigma = mad_std(scidata)
    daofind = DAOStarFinder(fwhm = 0.8 * quasar['FWHM'], threshold=5.*bkg_sigma, roundlo = -0.15, roundhi = 0.15)
    sources = daofind.find_stars(scidata)


    # Shift the source coordinates to the actual image coordinates

    #sources['xcentroid'] += xl
    #sources['ycentroid'] += yl


    # If no sources found, skip iteration
    
    if len(sources) <= 0:
        return

    print(len(sources))

    
    # Calculate the FWHM of each identified source, and append them into a column that is added to the source table
    # Splices the data array for the quasar, with the alleged centroid at the center
    # Fits a 2D Gaussian curve onto the array, and uses the relation between sigma and fwhm

    FWHM = []
    stddev_list = []
    
    for j in range(len(sources)):
        width = int(np.sqrt(sources['npix'][j]))
        #print("%d   %d   %f   %f" % (j, width, sources['xcentroid'][j], sources['ycentroid'][j]))
        data = scidata[int(sources['ycentroid'][j] - width/2) : int(sources['ycentroid'][j] + width/2) + 1, int(sources['xcentroid'][j] - width/2) : int(sources['xcentroid'][j] + width/2) + 1]


        """
        plt.imshow(data, origin='lower', interpolation='nearest', cmap='viridis')
        plt.show()
        plt.pause(3)
        """

        gauss = 0

        if(np.ma.count(data) >= 7):
            gauss = photutils.fit_2dgaussian(data, mask = None)

        fwhm = 0
        if gauss != 0:
            fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(gauss.x_stddev)
        FWHM.append(fwhm)

        if gauss == 0:
            stddev_list.append(0)
        else:
            stddev_list.append(np.sqrt(gauss.x_stddev))

    column = Column(FWHM, name = 'FWHM')
    sources.add_column(column)
    column = Column(stddev_list, name = 'stddev')
    sources.add_column(column)
    

    #print(sources)


    # Helper methods for determining differences between PSF source and QSO
    
    def distance1(x, y):
        return math.sqrt((x - quasar['xcentroid']) ** 2 + (y - quasar['ycentroid']) ** 2)

    def fwhmdiff(fwhm):
        return (fwhm - quasar['FWHM'])

    def lumdiff(flux):
        return (flux - quasar['flux'])

    def xdiff(x):
        return quasar['xcentroid'] - x

    def ydiff(y):
        return quasar['ycentroid'] - y



    distset = []
    fwhmset = []
    lumset = []
    

    # Standardize the sources and calculate the best source by combining distance, fwhm difference, and peak flux difference

    for j in range(len(sources)):
        d = distance1(sources['xcentroid'][j], sources['ycentroid'][j])
        distset.append(d)
        f = fwhmdiff(sources['FWHM'][j])
        fwhmset.append(f)
        l = lumdiff(sources['flux'][j])
        lumset.append(l)


    dstd = np.std(distset)
    fstd = np.std(fwhmset)
    lstd = np.std(lumset)


    
    # Weight of the three variables places FWHM difference as most important, flux difference as next important, and distance as least important

    psflist = []
    indexlist = []
    zlist = []

    for j in range(len(sources)):
        z = 1/3 * abs(distance1(sources['xcentroid'][j], sources['ycentroid'][j])/(dstd)) + 2 * abs(fwhmdiff(sources['FWHM'][j])/(fstd)) + 2/3 * abs(lumdiff(sources['flux'][j])/(lstd))

        """
        tempcut = scidata[int(sources['ycentroid'][j]) - chunk_size : int(sources['ycentroid'][j]) + chunk_size, int(sources['xcentroid'][j]) - chunk_size : int(sources['xcentroid'][j]) + chunk_size]


        s1 = daofind.find_stars(tempcut)

        if len(s1) > 1:
            print(len(s1))
            continue
        """
        
        #print(str(z))
        if z > 0 and math.sqrt((sources['xcentroid'][j] - quasar['xcentroid'])**2 + (sources['ycentroid'][j] - quasar['ycentroid'])**2) > 1 and sources['xcentroid'][j] - chunk_size >= 0 and sources['ycentroid'][j] - chunk_size >= 0 and sources['xcentroid'][j] + chunk_size < 2048 and sources['ycentroid'][j] + chunk_size < 1489:
            #print(str(z))
            ztot = z
            psfindex = j


            # If the list contains less than 5 suitable sources that satisfy all the conditions, then directly add to list
            # If the list already contains 5 sources, then check if the current source has a lower Z value than the source with the largest Z value
            
            if len(psflist) < 9 and z < 3:
                #print(True)
                psflist.append(tuplet(j, z))
                indexlist.append(j)
                zlist.append(z)
            else:
                if len(psflist) > 9 and z < max(zlist) and z < 3:
                    faker = psflist.remove(psf.getZ(max(zlist)))
                    indexlist.remove(faker.getIndex())
                    zlist.remove(faker.getZ())
                    psflist.append(tuplet(j, z))
                    indexlist.append(j)
                    zlist.append(z)

                

    # If no suitable PSF sources are found that satisfy the boundaries, then the file is not used
    
    if len(psflist) == 0:
        print("FAILURE")
        return

    
    PSFlist = []
    FWHMlist = []
    psf = 0

    fwhmdiff = 10000000
    for j in indexlist:
        if abs(sources['FWHM'][j] - quasar['FWHM']) < fwhmdiff:
            psf = sources[j]
            fwhmdiff = abs(sources['FWHM'][j] - quasar['FWHM'])

    """
        PSFlist.append(sources[j])
        FWHMlist.append(sources['FWHM'][j])


    medFWHM = 0
    if len(FWHMlist) % 2 == 0:
        medFWHM = FWHMlist[len(FWHMlist) // 2 - 1]
    else:
        medFWHM = np.median(FWHMlist)


    for k in PSFlist:
        if k['FWHM'] == medFWHM:
            print(True)
            psf = k
            break

    """

    # Interpolates the PSF source to its actual centroid, upsampling it 2x

    #try:
    chunk_size = 50

    print(psf)

    pXc = psf['xcentroid']
    pYc = psf['ycentroid']
    print("%f,  %f" % (pXc, pYc))
    xr = np.arange(int(pXc) - chunk_size - 5, int(pXc) + chunk_size + 6)
    yr = np.arange(int(pYc) - chunk_size - 5, int(pYc) + chunk_size + 6)


    preshift = scidata[int(pYc) - chunk_size - 5 : int(pYc) + chunk_size + 6, int(pXc) - chunk_size - 5: int(pXc) + chunk_size + 6]
    print(np.shape(preshift))

    shifted = []
    spline = interpolate.interp2d(xr, yr, preshift)
    xrf = np.arange(pXc - chunk_size, pXc + chunk_size + 1, 1)
    yrf = np.arange(pYc - chunk_size, pYc + chunk_size + 1, 1)

    if len(xrf) > 101:
        xrf = xrf[:-1].copy()
    if len(yrf) > 101:
        yrf = yrf[:-1].copy()

    shifted = spline(xrf, yrf)

    mean1, median1, std1 = sigma_clipped_stats(shifted, sigma=3.0, iters=5)
    shifted = checkOutter(shifted, mean1, std1)


    #qsocut = fits.open('c:/Research Project/Final Quasar Cut/' + str(i) + '_QSO.fit')
    qsocut = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color +  '_MG.fit')
    #qsocut = fits.open('c:/Research Project/Reference Quasar Cut/' + str(i) + '_REF.fit')
    qsodata = qsocut[0].data.astype(float)
    qsodata -= median

    #qsocount = photoncount(21, 21, qsigma, qsodata)
    #qsodata /= qsocount  #quasar['peak']

    #psfcount = photoncount(21, 21, psf['stddev'], shifted)
    #shifted /= np.max(shifted)  #psf['peak']
    #shifted *= np.max(qsodata)



    gauss_fit = photutils.fit_2dgaussian(shifted[40 : 61, 40 : 61], mask = None)
    fit_fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(abs(gauss_fit.x_stddev))
    #print(fit_fwhm)

    print("%f, %f" % (quasar['FWHM'], fit_fwhm))
    ffwhm = max(quasar['FWHM'], fit_fwhm)
    ffphoton_5sig = photonCount(50, 50, ffwhm, shifted)

    """
    qsophoton_4sig = photonCount(50, 50, ffwhm, qsodata)
    for j in range(len(shifted)):
        for k in range(len(shifted)):
            if distance(50, 50, j, k) < 4 * ffwhm:
                shifted[j][k] /= ffphoton_4sig
                shifted[j][k] *= qsophoton_4sig

    """


    line_data = linecache.getline('Full Data.txt', i).split()
    gmag = float(line_data[6])

    
    try:
        multiplier = 10**(gmag / (-2.5)) * 10**8 * hdulist[0].header['FLUX20']
        shifted /= ffphoton_5sig
        shifted *= multiplier
    except:
        #final_fit *= qsodata[50, 50]
        return
    
    residue = qsodata - shifted

    """
    mean, med, std = sigma_clipped_stats(residue, sigma=3.0, iters=5)

    check = False
    print("%f, %f" % (med, std))
    threshold = mean + 3 * std
    for j in range(42):
        for k in range(42):
            if shifted[k, j] > threshold:
                #print("Over")
                check = checkNoise(j, k, 21, 21, residue)
                if check == True:
                    #print("Over True")
                    return

    """ 

    #fits.writeto('Reference PSF Cut/' + str(i) + '_PSF.fit', shifted, hdulist[0].header, clobber = True)
    #fits.writeto('Test PSF Subtract/' + str(i) + '_1.fit', residue, hdulist[0].header, clobber = True)
    fits.writeto('/data/marvels/billzhu/MG II PSF Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_PSF.fit', shifted, hdulist[0].header, clobber = True)
    fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB.fit', residue, hdulist[0].header, clobber = True)
    #fits.writeto('Reference Subtract/' + str(i) + '_SUB.fit', residue, hdulist[0].header, clobber = True)
    print('\n')



    #except:
    #    return

    return



    
# Code that opens up a maximum of 8 processes for concurrent execution
    
if __name__ == '__main__':
    #multiprocessing.set_start_method('spawn')
    jobs = []
    """
    for j in range(1, 501):
        
        print(j)
        try:
            begin(j)
        except:
            print("Error: Unable to process")
        
        

        print(j)
        try:
            p1 = Process(target = begin, args = (j,))
            jobs.append(p1)
            p1.start()
            if len(job) > 8:
                time.sleep(3)
            #print('parent process:', os.getppid())
            #print('process id:', os.getpid())
        except:
            print("Error: Unable to process")
        
    """


    #dirs = os.listdir('/data/marvels/billzhu/2175 Quasar Cut/')
    #dirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/')

    gdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/g/')
    rdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/r/')
    idirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/i/')
    zdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/z/')
    udirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/u/')
    
    #check_dirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/')
    rangelist = []
    #rangelist.append(gdirs[0].split('_')[0])

    
    for d in gdirs:
        index = d.split('_')[0]
        #if str(index) + '_SUB.fit' not in check_dirs:
        rangelist.append(index)

    """
    for d in rdirs:
        index = d.split('_')[0]
        #if str(index) + '_SUB.fit' not in check_dirs:
        rangelist.append(index)

    for d in idirs:
        index = d.split('_')[0]
        #if str(index) + '_SUB.fit' not in check_dirs:
        rangelist.append(index)

    for d in zdirs:
        index = d.split('_')[0]
        #if str(index) + '_SUB.fit' not in check_dirs:
        rangelist.append(index)

    for d in udirs:
        index = d.split('_')[0]
        #if str(index) + '_SUB.fit' not in check_dirs:
        rangelist.append(index)
    """

    print(len(rangelist))
    #try:
    pool = Pool(os.cpu_count())
    pool.map(begin, rangelist)
    #except:
    #print("Error: Unable to process file")






