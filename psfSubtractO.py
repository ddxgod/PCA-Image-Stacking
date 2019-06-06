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
import multiprocessing

import photutils
from photutils import DAOStarFinder
from photutils import IRAFStarFinder
from photutils import subtract_psf
from photutils import centroid_2dg
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from functools import cmp_to_key


def inbounds(x, y):
    if x > 0 and x < 2048 and y > 0 and y < 1489:
        return True
    return False

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


PSF = []
reader = open('Pixel Coordinates 50000.txt', 'r')

def begin(i):
    chunk_size = 10

    #for j in range(3):
    #    num, x, y = reader.readline().split()

    num, x, y = reader.readline().split()
    x = int(x)
    y = int(y)
    print("%f   %f" % (x, y))
    print(str(i))
    
    filename = 'Test Data Extract/' + str(i) + '.fit'
    hdulist = fits.open(filename)
    half = 300
    scidata = hdulist[0].data.astype(float)
    mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
    
    if x + half + chunk_size > 2048:
        filler = np.array([float(median)] * len(scidata))
        for j in range(10):
            scidata = np.insert(scidata, len(scidata[0]), filler, 1)

    if x - half - chunk_size < 0:
        x += chunk_size
        filler = np.array([float(median)] * len(scidata))
        for j in range(10):
            scidata = np.insert(scidata, 0, filler, 1)

    if y + half + chunk_size > 1489:
        filler = np.array([float(median)] * len(scidata[0]))
        for j in range(10):
            scidata = np.insert(scidata, len(scidata), filler, 0)

    if y - half - chunk_size < 0:
        y += chunk_size
        filler = np.array([float(median)] * len(scidata[0]))
        for j in range(10):
            scidata = np.insert(scidata, 0, filler, 0)

    scidata -= median    
    psfindex = -1
    quasar = 0
    bkg_sigma = mad_std(scidata)
    daofind = DAOStarFinder(fwhm = 2., threshold = 3.*bkg_sigma)
    sources = daofind(scidata[y - 10 : y + 10, x - 10 : x + 10])

    sources['xcentroid'] += x - 10
    sources['ycentroid'] += y - 10
    
    FWHM = np.empty([len(sources)])
    column = Column(FWHM, name = 'FWHM')
    sources.add_column(column)

    # Find the quasar and calculate its FWHM
    
    for j in range(len(sources)):
        if abs(round(sources['xcentroid'][j]) - x) < 2 and abs(round(sources['ycentroid'][j]) - y) < 2:
            quasar = sources[j]
            width = int(np.sqrt(sources['npix'][j]))
            #print("%d   %d   %f   %f" % (j, width, sources['xcentroid'][j], sources['ycentroid'][j]))
            data = scidata[int(sources['ycentroid'][j] - width/2) : int(sources['ycentroid'][j] + width/2), int(sources['xcentroid'][j] - width/2) : int(sources['xcentroid'][j] + width/2)]
            gauss = 0
            
            if(np.ma.count(data) >= 7):
                gauss = photutils.fit_2dgaussian(data, mask = None)
            
            fwhm = 0
            if gauss != 0:
                fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(gauss.x_stddev**2 + gauss.y_stddev**2)
                quasar['FWHM'] = fwhm
                #print(quasar['FWHM'])
                break

    ztot = 10000
    print(quasar)

    if quasar == 0:
        return


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

    bkg_sigma = mad_std(scidata)
    daofind = DAOStarFinder(fwhm = quasar['FWHM'], threshold=3.*bkg_sigma, roundlo = -0.15, roundhi = 0.15)
    sources = daofind.find_stars(image)


    # Shift the source coordinates to the actual image coordinates

    sources['xcentroid'] += xl
    sources['ycentroid'] += yl

    #print(sources)

    # If no sources found, go to next iteration with larger dimensions
    if len(sources) <= 0:
        return

    # Calculate the FWHM of each identified source, and append them into a column that is added to the source table
    # Splices the data array for the quasar, with the alleged centroid at the center
    # Fits a 2D Gaussian curve onto the array, and uses the relation between sigma and fwhm

    FWHM = []
    for j in range(len(sources)):
        width = int(np.sqrt(sources['npix'][j]))
        #print("%d   %d   %f   %f" % (j, width, sources['xcentroid'][j], sources['ycentroid'][j]))
        data = scidata[int(sources['ycentroid'][j] - width/2) - 1 : int(sources['ycentroid'][j] + width/2) + 1, int(sources['xcentroid'][j] - width/2) - 1 : int(sources['xcentroid'][j] + width/2) + 1]
        gauss = 0

        if(np.ma.count(data) >= 7):
            gauss = photutils.fit_2dgaussian(data, mask = None)

        fwhm = 0
        if gauss != 0:
            fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(gauss.x_stddev**2 + gauss.y_stddev**2)
        FWHM.append(fwhm)

    column = Column(FWHM, name = 'FWHM')
    sources.add_column(column)


    def distance(x, y):
        return math.sqrt((x - quasar['xcentroid']) ** 2 + (y - quasar['ycentroid']) ** 2)

    def fwhmdiff(fwhm):
        return (fwhm - quasar['FWHM'])

    def lumdiff(flux):
        return (flux - quasar['flux'])

    def xdiff(x):
        return quasar['xcentroid'] - x

    def ydiff(y):
        return quasar['ycentroid'] - y


    davg = 0
    favg = 0
    lavg = 0
    distset = []
    fwhmset = []
    lumset = []

    # Standardize the sources and calculate the best source by combining distance, fwhm difference, and peak flux difference

    for j in range(len(sources)):
        d = distance(sources['xcentroid'][j], sources['ycentroid'][j])
        distset.append(d)
        davg += d
        f = fwhmdiff(sources['FWHM'][j])
        fwhmset.append(f)
        favg += f
        l = lumdiff(sources['flux'][j])
        lumset.append(l)
        lavg += l

    davg /= len(sources)
    favg /= len(sources)
    lavg /= len(sources)
    dstd = np.std(distset)
    fstd = np.std(fwhmset)
    lstd = np.std(lumset)

    # Weight of the three variables places FWHM difference as most important, flux difference as next important, and distance as least important

    psflist = []
    indexlist = []
    zlist = []

    for j in range(len(sources)):
        z = 1/2 * abs(distance(sources['xcentroid'][j], sources['ycentroid'][j])/(dstd)) + 4/3 * abs(fwhmdiff(sources['FWHM'][j])/(fstd)) + 2/3 * abs(lumdiff(sources['flux'][j])/(lstd))
        #print(str(z) + " " + str(abs(sources['peak'][j] - quasar['peak'])))
        if z > 0 and sources['peak'][j] > 0.7 * quasar['peak'] and inbounds(sources['xcentroid'][j], sources['ycentroid'][j]) and math.sqrt((sources['xcentroid'][j] - quasar['xcentroid'])**2 + (sources['ycentroid'][j] - quasar['ycentroid'])**2) > 1:
            #print(str(z))
            ztot = z
            psfindex = j

            if len(psflist) < 5 and z < 4:
                psflist.append(tuplet(j, z))
                indexlist.append(j)
                zlist.append(z)
            else:
                if len(psflist) > 5 and z < max(zlist) and z < 4:
                    faker = psflist.remove(psf.getZ(max(zlist)))
                    indexlist.remove(faker.getIndex())
                    zlist.remove(faker.getZ())
                    psflist.append(tuple(j, z))
                    indexlist.append(j)
                    zlist.append(z)

                

    if len(psflist) == 0:
        return


    stdev = 10000000
    residue = 0
    cutout = 0

    print(indexlist)
    for j in indexlist:       
        psf = sources[j]
        chunk_size = 10

        # Find the actual centroid of the PSF using 2D Gaussian Fitting, since DAOStarFinder is inaccurate

        print(psf)
        """
        preshift = scidata[int(psf['ycentroid'] - chunk_size) : int(psf['ycentroid'] + chunk_size + 1), int(psf['xcentroid'] - chunk_size) : int(psf['xcentroid'] + chunk_size + 1)]
        mean, med2, std = sigma_clipped_stats(preshift, sigma=3.0, iters=5)
        mask = [[False for x in range(int(chunk_size*2) + 1)] for y in range(int(chunk_size*2) + 1)] 
        for j in range(0, int(chunk_size*2 + 1)):
            for k in range(0, int(chunk_size*2 + 1)):
                if scidata[int(psf['ycentroid'] + k - chunk_size), int(psf['xcentroid'] + j - chunk_size)] < med2:
                    mask[j][k] = True

        #pXc, pYc = centroid_2dg(preshift, mask = mask, error = None)
        #pXc += int(psf['xcentroid']) - chunk_size
        #pYc += int(psf['ycentroid']) - chunk_size
        """
        
        pXc = psf['xcentroid']
        pYc = psf['ycentroid']
        print("%f,  %f" % (pXc, pYc))
        xr = np.arange(int(pXc) - chunk_size, int(pXc) + chunk_size + 1)
        yr = np.arange(int(pYc) - chunk_size, int(pYc) + chunk_size + 1)

        
        preshift = scidata[int(pYc) - chunk_size : int(pYc) + chunk_size + 1, int(pXc) - chunk_size : int(pXc) + chunk_size + 1]

        shifted = []
        spline = interpolate.interp2d(xr, yr, preshift)
        xrf = np.arange(pXc - chunk_size, pXc + chunk_size + 1, 0.5)
        yrf = np.arange(pYc - chunk_size, pYc + chunk_size + 1, 0.5)

        if len(xrf) > 42:
            xrf = xrf[:-1].copy()
        if len(yrf) > 42:
            yrf = yrf[:-1].copy()

        shifted = spline(xrf, yrf)


        qsocut = fits.open('c:/Research Project/Final Quasar Cut/' + str(i) + '_QSO.fit')
        qsodata = qsocut[0].data.astype(float)
        qsodata -= median
        qsodata /= qsodata[21, 21]  #quasar['peak']
        shifted /= shifted[21, 21]  #psf['peak']

        res = qsodata - shifted

        mean, med, std = sigma_clipped_stats(res, sigma=3.0, iters=5)

        if std < stdev:
            residue = res
            stdev = std
            cutout = shifted

        
        print(std)
    
    fits.writeto('Test PSF Cut/' + str(i) + '_PSF.fit', cutout, hdulist[0].header, clobber = True)
    fits.writeto('Test PSF Subtract/' + str(i) + '_1.fit', residue, hdulist[0].header, clobber = True)
    #print(stdev)
    print('\n')

    

    PSF.append(psf)
    #sources = sorted(sources, key = lambda a : 1/2*abs(distance(a['xcentroid'], a['ycentroid'])/(dstd)) + abs(fwhmdiff(a['FWHM'])/(fstd)) + abs(lumdiff(a['flux'])/(lstd)))


    
# Code that opens up a maximum of 8 processes for concurrent execution
    
if __name__ == '__main__':
    multiprocessing.set_start_method('spawn')
    jobs = []
    for j in range(1, 11):
        
        print(j)
        try:
            begin(j)
        except:
            print("Error: Unable to process")
        
        """

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
        
    
writer = open('psf sources 1.txt', 'w')
for psf in PSF:
    writer.write(str(psf))
writer.close()





