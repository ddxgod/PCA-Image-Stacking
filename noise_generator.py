# make Gaussian and Poisson noise images
from photutils.datasets import make_noise_image
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from astropy.io import fits
from astropy.table import Table
import multiprocessing

import astropy.units as u
import astropy.coordinates as coord
import sys


def gaussian_2d(x, y, x0, y0, xsig, ysig):
    return np.exp(-0.5*(((x-x0) / xsig)**2 + ((y-y0) / ysig)**2))


impact_parameter=[0.2,0.8,1.2,1.6,2.0]
FWHM_qso = 1.57
plate_scale = 1
for i in range(0,5):

    x_qso= 50
    y_qso= 50
    x_host= x_qso+impact_parameter[i]
    y_host= 51
    x = np.arange(0, 101, plate_scale)
    y = np.arange(0, 101, plate_scale)
    X, Y = np.meshgrid(x, y)
    QSO=gaussian_2d(X, Y, x_qso, y_qso, FWHM_qso/2.355, FWHM_qso/2.355)
    host=gaussian_2d(X, Y, x_host, y_host, FWHM_qso/2.355, FWHM_qso/2.355)
    plt.imshow(QSO + make_noise_image((101, 101), type='poisson', mean=0.001))
    print(x_qso)
    print(y_qso)
    print(X)
    print(Y)


shape = (101, 101)
image1 = make_noise_image(shape, type='gaussian', mean=0., stddev=5.)
image2 = make_noise_image(shape, type='poisson', mean=5.)

# plot the images
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
ax1.imshow(image1, origin='lower', interpolation='nearest')
ax1.set_title('Gaussian noise ($\mu=0$, $\sigma=5.$)')
ax2.imshow(image2, origin='lower', interpolation='nearest')
ax2.set_title('Poisson noise ($\mu=5$)')


plt.show()