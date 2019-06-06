from astropy.io import fits
import scipy as sp
import numpy as np
import glob
import os
import string
import decimal
import matplotlib.pyplot as plt
import fileinput
import math
import photutils

from astropy.table import Table
from photutils.datasets import make_random_gaussians, make_noise_image
from photutils.datasets import make_gaussian_sources
from photutils.detection import IRAFStarFinder
from photutils.psf import IntegratedGaussianPRF, DAOGroup
from photutils.background import MMMBackground, MADStdBackgroundRMS
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.stats import gaussian_sigma_to_fwhm
from astropy.stats import sigma_clipped_stats
from photutils.psf import IterativelySubtractedPSFPhotometry


hdulist = fits.open('96998-g_MG.fit')
image = hdulist[0].data.astype(float)
mean, median, stddev = sigma_clipped_stats(image, sigma=3.0, iters=5)
image -= median
sigma_psf = 1.1
bkgrms = MADStdBackgroundRMS()
std = bkgrms.calc_background_rms(image)
iraffind = IRAFStarFinder(threshold=5*std,
                          fwhm=sigma_psf*gaussian_sigma_to_fwhm,
                          minsep_fwhm=0.01, roundhi=5.0, roundlo=-5.0,
                          sharplo=0.0, sharphi=2.0)

gauss = 0
            
if(np.ma.count(image) >= 7):
    gauss = photutils.fit_2dgaussian(image[40 : 61, 40 : 61], mask = None)

fwhm = 0
if gauss != 0:
    fwhm = abs(gauss.x_stddev) * gaussian_sigma_to_fwhm

print(fwhm)

    
daogroup = DAOGroup(5.0*sigma_psf*gaussian_sigma_to_fwhm)
mmm_bkg = MMMBackground()
fitter = LevMarLSQFitter()
psf_model = IntegratedGaussianPRF(sigma=abs(gauss.x_stddev))

fitshape = (int(3. * fwhm), int(3. * fwhm))
if int(3. * fwhm) % 2 == 0:
    fitshape = (int(3. * fwhm) + 1, int(3. * fwhm) + 1)

photometry = IterativelySubtractedPSFPhotometry(finder=iraffind,
                                                group_maker=daogroup,
                                                bkg_estimator=mmm_bkg,
                                                psf_model=psf_model,
                                                fitter=LevMarLSQFitter(),
                                                niters=1, fitshape=fitshape)
result_tab = photometry(image=image)
residual_image = photometry.get_residual_image()

#print(result_tab)

#plt.imshow(result_tab.astype(float), cmap='viridis', aspect=1, interpolation='nearest',
              # origin='lower')
#plt.show()
hdu = fits.PrimaryHDU(image)
fits.writeto('LastDitch96998.fit', residual_image, hdulist[0].header, clobber = True)
