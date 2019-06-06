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
from sklearn.decomposition import NMF, PCA, IncrementalPCA
import astropy.units as u
import astropy.coordinates as coord

#from dustmaps.sfd import SFDQuery

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



# Method that checks if a particular residue is noise or not by drawing a ring of similar distance as the "bright" point around the centroid

def checkNoise(x, y, qX, qY, data):
    halo = []
    for i in range(42):
        for j in range(42):
            if abs(distance(i, j, qX, qY) - distance(x, y, qX, qY)) <= 2 and distance(i, j, x, y) >= 2:
                halo.append(data[j, i])
    mean, median, std = sigma_clipped_stats(halo, sigma=3.0, iters=5)

    if data[y, x] > mean + 3 * std:
        return True
    return False
                

    
# Method that calculates the total photon count within 3 sigma of the quasar centroid

def photonCount(xc, yc, sigma, data):
    count = 0
    for i in range(len(data)):
        for j in range(len(data)):
            if distance(i, j, xc, yc) <= sigma:
                count += data[i][j]
    return count



# Method that checks if a bright source is outside of a certain radius

def checkOutter(data, mean, stddev):
    for i in range(len(data)):
        for j in range(len(data[0])):
            if data[i][j] > mean + stddev and distance(i, j, 50, 50) > 25:
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
    #i = int(index)
    i = int(index.split('-')[0])
    mgi = int(index.split('-')[1])
    color = index.split('-')[2]
    #mgi = int(index.split('-')[1])
    
    try:
        print(index)
        #filename = 'Test Data Extract/' + str(i) + '.fit'
        #filename = str(i) + '-g.fit'


        #filename = '/data/marvels/billzhu/2175 Dataset/' + str(index) + '-g.fit'
        filename = '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/' + color + '/' + str(index) + '.fit'
        #print(filename)
        #filename = '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '.fit'
        hdulist = fits.open(filename)

        #qlist = fits.open('MG II Test Cut/' + str(i) + '_MG.fit')

        #qlist = fits.open('/data/marvels/billzhu/2175 Quasar Cut/' + str(index) + '_DUST.fit')
        qlist = fits.open('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(index) + '_REF.fit')
        #qlist = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')
        x = qlist[0].header['XCOORD']
        y = qlist[0].header['YCOORD']
        #print("%f, %f" % (x, y))
        qlist.close()

    except:
        print("No coordinates")
        return

    # Save some frickin time

        
    
    half = 700
    scidata = hdulist[0].data.astype(float)
    mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
    
    if x + chunk_size > 2048:
        filler = np.array([float(median)] * len(scidata))
        for j in range(chunk_size):
            scidata = np.insert(scidata, len(scidata[0]), filler, 1)

    if x - chunk_size < 0:
        x += chunk_size
        filler = np.array([float(median)] * len(scidata))
        for j in range(chunk_size):
            scidata = np.insert(scidata, 0, filler, 1)

    if y + chunk_size > 1489:
        filler = np.array([float(median)] * len(scidata[0]))
        for j in range(chunk_size):
            scidata = np.insert(scidata, len(scidata), filler, 0)

    if y - chunk_size < 0:
        y += chunk_size
        filler = np.array([float(median)] * len(scidata[0]))
        for j in range(chunk_size):
            scidata = np.insert(scidata, 0, filler, 0)


    #if 'SKY' in hdulist[0].header.keys():
    #    scidata -= float(hdulist[0].header['SOFTBIAS'])
    #    scidata -= float(hdulist[0].header['SKY'])
    #else:
    scidata -= median
        
    psfindex = -1
    quasar = 0
    bkg_sigma = mad_std(scidata)

    
    # DAOStarFinder algorithm that finds all sources greater than 3 sigma above the background value, with strict roundness parameters
    
    daofind = DAOStarFinder(fwhm = 2., threshold = 5.*bkg_sigma)
    #print("%f, %f" % (x, y))
    sources = daofind(scidata[y - 10 : y + 10, x - 10 : x + 10])

    
    # Update coordinates of the sources
    
    sources['xcentroid'] += x - 10
    sources['ycentroid'] += y - 10

    
    #print(sources)

    
    # Create new column that contains the FWHM of each source for comparison later
    
    FWHM = np.empty([len(sources)])
    column = Column(FWHM, name = 'FWHM')
    sources.add_column(column)


    
    # Find the quasar and calculate its FWHM
    
    for j in range(len(sources)):
        if abs(round(sources['xcentroid'][j]) - x) < 3.0 and abs(round(sources['ycentroid'][j]) - y) < 3.0:
            quasar = sources[j]
            width = int(np.sqrt(sources['npix'][j]))
            #print("%d   %d   %f   %f" % (j, width, sources['xcentroid'][j], sources['ycentroid'][j]))
            data = scidata[int(sources['ycentroid'][j] - width - 1) : int(sources['ycentroid'][j] + width) + 2, int(sources['xcentroid'][j] - width - 1) : int(sources['xcentroid'][j] + width) + 2]


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
                qsigma = np.sqrt(gauss.x_stddev**2 + gauss.y_stddev**2)
                print(quasar['FWHM'])
                break

            
    ztot = 10000
    print(quasar)

    
    # If no quasar is found, the field image is deemed corrupt and not used
    
    if quasar == 0:
        return


    # Define cutout image limits, adjusted to field image boundaries as necessary i.e. x, y < 0 or > max x/y values
    
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
    daofind = DAOStarFinder(fwhm = quasar['FWHM'], threshold=7.*bkg_sigma, roundlo = -0.20, roundhi = 0.20)
    sources = daofind.find_stars(scidata)
    #print(len(sources))


    #qsocut = fits.open('/data/marvels/billzhu/2175 Quasar Cut/' + str(index) + '_DUST.fit')
    qsocut = fits.open('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/' + color + '/' + str(index) + '_REF.fit')
    #qsocut = fits.open('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_MG.fit')
    qsodata = qsocut[0].data.astype(float)


    

    # Shift the source coordinates to the actual image coordinates

    #sources['xcentroid'] += xl
    #sources['ycentroid'] += yl


    # If no sources found, skip iteration
    
    if len(sources) <= 0:
        return


    # Interpolates the PSF sources to their actual centroids, upsampling it 2x, and adds it to large array for PCA

    temp = 1
    largearr = []
    for j in range(len(sources)):
        if abs(sources['xcentroid'][j] - quasar['xcentroid']) < 2 and abs(sources['ycentroid'][j] - quasar['ycentroid']) < 2:
            continue
        
        chunk_size = 50

        pXc = sources['xcentroid'][j]
        pYc = sources['ycentroid'][j]
        #print("%f,  %f" % (pXc, pYc))

        xr = np.arange(int(pXc) - chunk_size - 5, int(pXc) + chunk_size + 6)
        yr = np.arange(int(pYc) - chunk_size - 5, int(pYc) + chunk_size + 6)


        preshift = scidata[int(pYc) - chunk_size - 5 : int(pYc) + chunk_size + 6, int(pXc) - chunk_size - 5: int(pXc) + chunk_size + 6]

        shifted = []
        try:
            spline = interpolate.interp2d(xr, yr, preshift)
            xrf = np.arange(pXc - chunk_size, pXc + chunk_size + 1, 1)
            yrf = np.arange(pYc - chunk_size, pYc + chunk_size + 1, 1)
        except:
            #print("ERROR")
            continue
        
        if len(xrf) > 101:
            xrf = xrf[:-1].copy()
        if len(yrf) > 101:
            yrf = yrf[:-1].copy()

        shifted = spline(xrf, yrf)
        cont = False


        
        
        # Safety that screens out images with multiple sources by checking incremental means
        # CHECK DISCONTINUED DUE TO INCOMPLETENESS / SOME ERRORS

        """
        meanarr = []
        for k in range(0, 5):
            tempcut = list(shifted[20 - 4 * k : 21 + 4 * k, 20 - 4 * k : 21 + 4 * k])
            mean1 = np.mean(tempcut)
            #print(mean1)
            if len(meanarr) > 0 and mean1 > meanarr[len(meanarr) - 1]:
                cont = True
                #print(temp)
                #fits.writeto(str(temp) + '.fit', shifted, clobber = True)
                #temp += 1
                break
            
            meanarr.append(mean1)
        """

        

        # Originally discontinued, but upon closer inspection, the same source finder parameters is used as in the original source search, thus sources found
        # will be the same, i.e. ideal source finder
        
        #check_source = daofind.find_stars(preshift)
        #if len(check_source) > 1:
        #    continue

        
        # If the source has a weird shape i.e. due to gravitational lensing, then check if the maximum pixel is within 2 pixels of the center to ensure consistency

        mean1, median1, std1 = sigma_clipped_stats(shifted, sigma=3.0, iters=5)
        daofind = DAOStarFinder(fwhm = 2, threshold=3.0*bkg_sigma)
        sources1 = daofind.find_stars(shifted)
        cont = checkInner(shifted, sources1)

        if cont == True:
            continue
        
        shifted = checkOutter(shifted, mean1, std1)

        """
        max_coords = np.unravel_index(shifted.argmax(), shifted.shape)
        max_coords = list(max_coords)
        #print(max_coords)
        for k in range(len(max_coords)//2):
            yt = max_coords[2 * k]
            xt = max_coords[2 * k + 1]
            #print("%f, %f" % (xt, yt))
            if distance(xt, yt, 20, 20) > 4:
                cont = True
                break
        """

        


        #fits.writeto(str(temp) + '.fit', shifted, clobber = True)
        #print(temp)
        #print(meanarr)
        #shifted /= np.max(shifted)
        #shifted *= np.max(qsodata)
        largearr.append(np.reshape(shifted, 10201))

    largearr = np.array(largearr)
    print(np.shape(largearr))


    # Set number of components in PCA, use incremental PCA (IPCA) due to high efficiency and speed
    
    numcomp = 20


    # Need a healthy number of sources to make the PSF fitting in order to decrease noise, setting at 5% threshold
    
    if len(largearr) < 8:
        return

    print(numcomp)
    mean_vector = []

    #print(np.shape(largearr))

    try:
        for j in range(0, 10201):
            mean_vector.append(np.mean(largearr[:, j]))
    except:
        print("NO SOURCE FOUND")
        return

    largearr -= mean_vector
        
    ipca = IncrementalPCA(n_components=numcomp)
    ipca.fit(largearr)
    ipca_comp = ipca.components_
    #print(np.shape(ipca_comp))
    ipca_comp = ipca_comp.T
    #print(ipca_comp)

    #print(np.shape(largearr[0, :]))
    #print(np.shape(ipca_comp))
    total_res = 0
    max_median = 10000000



    """
    # Calculate optimal number of coefficients to be taken

    for p, take in enumerate([12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, 120]):
        if take > len(largearr)//4 * 3 or take < len(largearr) //2:
            continue
        
        totalres = 0
        for j in range(len(largearr)): 
            coeff = np.dot(largearr[j, :], ipca_comp[:, 0:take])
            fit = np.dot(ipca_comp[:, 0:take], coeff[0:take])
            resfit = largearr[j, :] - fit
            total_res += resfit

        m1, m2, s1 = sigma_clipped_stats(total_res, sigma=3.0, iters=5)

        if m2 < max_median:
            max_median = m1
            take_final = take

            # Lowest mean means lowest residual
            
            
            plt.imshow(np.reshape(fit, (42, 42)), origin='lower', interpolation='nearest', cmap='viridis')
            plt.show()
            plt.pause(2)
            plt.close()
    """


    

    #if 'SKY' in hdulist[0].header.keys():
    #    qsodata -= float(hdulist[0].header['SOFTBIAS'])
    #    qsodata -= float(hdulist[0].header['SKY'])
    #else:
    qsodata -= median


    take_final = 4
    

    # Final fitting of the first n components, as determined by take_final, into the quasar to build a PSF fit
    
    qsodata = np.reshape(qsodata, 10201)
    coeff = np.dot(qsodata, ipca_comp[:, 0:take_final])
    final_fit = np.dot(ipca_comp[:, 0:take_final], coeff[0:take_final])
    final_fit += mean_vector
    final_fit = np.reshape(final_fit, (101, 101))
    #final_fit /= len(largearr)
    qsodata = np.reshape(qsodata, (101, 101))

    """
    qx = np.arange(0, len(qsodata))
    qy = np.arange(0, len(qsodata))
    spline = interpolate.interp2d(qx, qy, qsodata)
    qxf = np.arange(0, len(qsodata), 0.1)
    qyf = np.arange(0, len(qsodata), 0.1)
    qsodata = spline(qxf, qyf)
    spline = interpolate.interp2d(qx, qy, final_fit)
    final_fit = spline(qxf, qyf)
    """

    gauss_fit = photutils.fit_2dgaussian(final_fit[40 : 61, 40 : 61], mask = None)
    fit_fwhm = 2*np.sqrt(2*np.log(2))*np.sqrt(gauss_fit.x_stddev)
    #print(fit_fwhm)

    
    #print("%f, %f" % (quasar['FWHM'], fit_fwhm))
    ffwhm = max(quasar['FWHM'], fit_fwhm)
    ffphoton_1sig = photonCount(50, 50, 2 * ffwhm, final_fit) 
    #qsophoton_1sig = photonCount(50, 50, 6, qsodata)

    """
    for j in range(len(final_fit)):
        for k in range(len(final_fit)):
            if distance(50, 50, j, k) < 3:
                final_fit[j][k] /= ffphoton_1sig
                final_fit[j][k] *= qsophoton_1sig
    """
    

    
    #final_fit /= ffphoton_1sig
    #final_fit *= qsophoton_1sig
    
    line_data = linecache.getline('Full Data.txt', i).split()

    if color == 'g':
        mag = float(line_data[6])
    if color == 'r':
        mag = float(line_data[8])
    if color == 'i':
        mag = float(line_data[10])
    if color == 'z':
        mag = float(line_data[12])
    if color == 'u':
        mag = float(line_data[4])
        

    
    try:
        multiplier = 10**(mag / (-2.5)) * 10**8 * hdulist[0].header['FLUX20']
        final_fit /= ffphoton_1sig
        final_fit *= multiplier
    except:
        #final_fit *= qsodata[50, 50]
        return
    

        
    """
    header = hdulist[0].header
    mag20 = header['flux20'] - median
    
    plt.imshow(qsodata, origin='lower', interpolation='nearest', cmap='viridis')
    plt.show()

    plt.imshow(totalfit, origin='lower', interpolation='nearest', cmap='viridis')
    plt.show()
    #plt.pause(3)
    """


    print("%f, %f" % (qsodata[50][50], final_fit[50][50]))
    residue = qsodata - final_fit

    """
    residue /= mag20

    f1 = fits.open('/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/4-g.fit')
    h1 = f1[0].header
    mean1, median1, stddev1 = sigma_clipped_stats(f1[0].data.astype(float), sigma=3.0, iters=5)
    mag20_1 = h1['flux20'] - median1
    residue *= mag20_1

    qdata = linecache.readline('Full Data.txt', i)
    c = coord.SkyCoord(ra = float(qdata[1]), dec = float(qdata[2]))

    sfd = SFDQuery()
    residue *= 10**(0.4 * sfd(c))
    """

    #plt.imshow(residue, origin='lower', interpolation='nearest', cmap='viridis')
    #plt.show()



    # Only used for reference quasars

    """
    for j in range(42):
        for k in range(42):
            if shifted[k, j] > threshold:
                #print("Over")
                check = checkNoise(j, k, 21, 21, residue)
                if check == True:
                    #print("Over True")
                    counter += 1

                if counter > 10:
                    return
    """

    try:
        #fits.writeto('/data/marvels/billzhu/2175 PSF Cut/' + str(index) + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
        #fits.writeto('/data/marvels/billzhu/2175 PSF Subtract/' + str(index) + '_SUB.fit', residue, hdulist[0].header, clobber = True)

        fits.writeto('/data/marvels/billzhu/Reference PSF Cut/0.37 - 0.55/' + color + '/' + str(index) + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
        fits.writeto('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + color + '/' + str(index) + '_SUB.fit', residue, hdulist[0].header, clobber = True)
        #fits.writeto('/data/marvels/billzhu/MG II PSF Cut/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
        #fits.writeto('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + color + '/' + str(i) + '-' + color + '_SUB.fit', residue, hdulist[0].header, clobber = True)

        #fits.writeto('Reference Subtract/' + str(i) + '_SUB.fit', residue, hdulist[0].header, clobber = True)
        #fits.writeto('Reference PSF Cut/' + str(i) + '_PSF.fit', final_fit, hdulist[0].header, clobber = True)
        print('\n')

        print("DONE TO BOTTOM")
    except:
        print('HEADER IS CORRUPT')



    
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

    #gdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/g/')
    #rdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/r/')
    #idirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/i/')
    #zdirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/z/')
    #udirs = os.listdir('/data/marvels/billzhu/MG II Quasar Cut/0.37 - 0.55/u/')


    gdirs = os.listdir('/data/marvels/billzhu/Reference Quasar Cut/0.37 - 0.55/g/')
    
    #check_dirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/')
    rangelist = []
    #rangelist.append('96624-g')

    
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





