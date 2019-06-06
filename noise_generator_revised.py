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
import itertools


def gaussian_2d(x, y, x0, y0, xsig, ysig):
    return np.exp(-0.5*(((x-x0) / xsig)**2 + ((y-y0) / ysig)**2))


impact_parameter=[0.2,0.8,1.2,1.6,2.0]
FWHM_qso = 1.57 / 0.396
plate_scale = 1
for i in range(0,5):

    x_qso= 50
    y_qso= 50
    x_host= x_qso+impact_parameter[i]
    y_host= 51
    x = np.arange(0, 101, plate_scale)
    y = np.arange(0, 101, plate_scale)
    X, Y = np.meshgrid(x, y)
    QSO=gaussian_2d(X, Y, x_qso, y_qso, FWHM_qso/2.355, FWHM_qso/2.355)*10
    host=gaussian_2d(X, Y, x_host, y_host, FWHM_qso/2.355, FWHM_qso/2.355)*10
    #plt.imshow(QSO + make_noise_image((101, 101), type='poisson', mean=0.001))
    '''
    add sky background and noise here
    '''
    sky = 1.0
    final_image = QSO+host+sky
    plt.imshow(final_image)
    plt.show()
    print(host.shape[0])
    print(host.shape[1])
    for k, j in itertools.product(range(final_image.shape[0]), range(final_image.shape[1])):
        arr=np.random.normal(final_image[k,j], final_image[k,j]**0.5, 1000)
        index_need=np.random.randint(arr.size, size=1)
        final_image[k, j]=arr[index_need]


    print(np.shape(final_image))
    if i == 4:
        fits.writeto('Sim.fit', np.array(final_image), clobber=True)
    #plt.show()
