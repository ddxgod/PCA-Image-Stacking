from astropy.io import fits
from astropy.table import Table
import scipy as sp
import scipy.optimize as opt
from scipy import interpolate
import numpy as np
np.set_printoptions(threshold=np.inf)
import string
import matplotlib.pyplot as plt
import math
import os
import linecache
import photutils
from photutils import DAOStarFinder
from photutils import subtract_psf
from photutils import centroid_2dg
from astropy.stats import mad_std
from astropy.stats import sigma_clipped_stats
from astropy.cosmology import WMAP9 as cosmo
from astropy.table import Table
from astropy.io import fits
from shutil import copyfile
import numpy as np
import scipy as sp
import math
import os
import multiprocessing
from multiprocessing import Pool



f1 = fits.open('QSObased_Trimmed_SDSS_DR7_107.fits')
mgtable = Table.read('QSObased_Trimmed_SDSS_DR7_107.fits')


def distance(ra, dec, xra, xdec):
    return math.sqrt((ra - xra)**2 + (dec - xdec)**2)



# Find the quasar with the lowest absorber redshift

mgdirs = os.listdir('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/')
#mgdirs = os.listdir('MG II Test Cut/')
lowest_qso = 0
lowest_redshift = 10000
lindex = 0

counter = 0

for f in mgdirs:
    index = int(f.split('_')[0])
    print(index)
    line = linecache.getline('Full Data.txt', index)
    ra = float(line.split()[1])
    dec = float(line.split()[2])
    #print("%f, %f" % (ra, dec))
    z = 100000
    
    for j in range(len(mgtable)):
        if mgtable['INDEX_QSO'][j] + 1 == index:
            z = mgtable['ZABS'][j][0]
            counter += 1
            print(z)
            break
    
    if z < lowest_redshift:
        lowest_redshift = z
        lindex = index
        lowest_qso = fits.open('/data/marvels/billzhu/MG II PSF Subtract/0.37 - 0.55/' + f)
        #lowest_qso = fits.open('MG II Test Cut/' + f)

print(counter)
prihdr = lowest_qso[0].header

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



tot_ra = raDegRowPix * 42
tot_dec = decDegColPix * 42

scale = cosmo.kpc_proper_per_arcmin(0.37) * 60.
x_kpc = scale * tot_ra
y_kpc = scale * tot_dec

print(lowest_redshift)

refdirs = os.listdir('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/')
        
# For each quasar, both absorption and reference, calculate the distance, kpc across the image, and scale

for f in refdirs:
    index = int(f.split('_')[0])

    if index == lindex:
        continue
    
    z = float(linecache.getline('Full Data.txt', index).split()[3])
    hdulist = fits.open('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/' + f)
    #hdulist = fits.open('MG II Test Cut/' + f)
    prihdr = hdulist[0].header

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

    dra = raDegRowPix * 42
    ddec = decDegColPix * 42
    qscale = cosmo.kpc_proper_per_arcmin(z) * 60
    x1_kpc = qscale * dra
    y1_kpc = qscale * ddec
    x_pix = 42 * x_kpc / x1_kpc
    y_pix = 42 * y_kpc / y1_kpc
    
    xr = np.arange(0, 42, 1)
    yr = np.arange(0, 42, 1)

    data = hdulist[0].data.astype(float)
    spline = interpolate.interp2d(xr, yr, data)
    xr1 = np.arange(0, 42, 0.05)
    yr1 = np.arange(0, 42, 0.05)
    dat1 = spline(xr1, yr1)
    print(np.shape(dat1))
    
    #print("%f, %f" % (x_pix, y_pix))
    try:
        dat1 = dat1[int(419 - 10 * y_pix) : int(419 + 10 * y_pix), int(419 - 10 * x_pix) : int(419 + 10 * x_pix)]

        print(np.shape(dat1))
        print("%f, %f" % (419 - 10 * x_pix, 419 - 10 * y_pix))
        xr1 = np.arange(0, 20 * x_pix, 1)
        yr1 = np.arange(0, 20 * y_pix, 1)
        print("%f, %f" % (10 * x_pix, 10 * y_pix))

        if len(xr1) > len(dat1):
            xr1 = xr1[:-1]
        if len(yr1) > len(dat1[0]):
            yr1 = yr1[:-1]

        print(np.shape(xr1))
        print(np.shape(yr1))
        print(np.shape(dat1))
        spline = interpolate.interp2d(xr1, yr1, dat1)
        xrf = np.arange(0, 20 * x_pix, 20 * x_pix / 42)
        yrf = np.arange(0, 20 * y_pix, 20 * y_pix / 42)

        if len(xrf) > 42:
            xrf = xrf[:-1]
        if len(yrf) > 42:
            yrf = yrf[:-1]
        #print(np.shape(xrf))
        #print(np.shape(yrf))
        #print(np.shape(dat1))
        cut = spline(xrf, yrf)

        fits.writeto('/data/marvels/billzhu/Reference PSF Subtract 2/0.37 - 0.55/' + f, cut, prihdr, clobber = True)
        #fits.writeto('MG II Rescale/' + f, cut, prihdr, clobber = True)
    except:
        print("NOPE")
    
    
