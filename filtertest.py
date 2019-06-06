from astropy.io import fits
import scipy as sp 
import numpy as np 

from scipy.ndimage.filters import gaussian_filter


f = fits.open('88773-r_MG.fit')
data = f[0].data
data = gaussian_filter(data, sigma=3)

fits.writeto('GaussianTest88773.fit', data, f[0].header, clobber=True)
