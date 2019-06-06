import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

def gaussian_2d(x, y, x0, y0, xsig, ysig):
    return np.exp(-0.5*(((x-x0) / xsig)**2 + ((y-y0) / ysig)**2))

#calculate the brightness of the background QSO and the host galaxy
F_qso=91*10**(-17)
F_host=18*10**(-17)
central_wave=7530
photo=1.98*10**(-8)/central_wave
D_tel=10

Ntotal_host=F_host*3.14*(D_tel*100/2)**2/photo
Ntotal_qso=F_qso*3.14*(D_tel*100/2)**2/photo
#CCD information
QE=0.85
trans=0.85
DC=3/3600
RN=5
#size of the object in the image
FWHM_qso=1.2
FWHM_host=1.2
plate_scale=0.25
Npix_host=3.14*(2*FWHM_host/2.355)**2/(plate_scale**2)
Npix_qso=3.14*(2*FWHM_qso/2.355)**2/(plate_scale**2)
#calculate teh brightness of the sky
F_sky=403*10**(-17)
Ntotal_sky=F_sky*3.14*(D_tel*100/2)**2*(plate_scale**2)/photo

#calculate the signal to noise for certain exposure time.
overall_eff=0.4

t=3000
Total_host=Ntotal_host*t*QE*trans*overall_eff
Total_qso=Ntotal_qso*t*QE*trans*overall_eff
SN=Total_host/(Total_host+Npix_host*(Ntotal_sky*t*QE*trans+DC*t+RN**2))**0.5

print(SN)


impact_parameter=[0.2,0.8,1.2,1.6,2.0]

for i in range(0,5):

    x_qso= -1
    y_qso= -1
    x_host= x_qso+impact_parameter[i]
    y_host= -1
    x = np.arange(-5.0, 5.0, plate_scale)
    y = np.arange(-5.0, 5.0, plate_scale)
    X, Y = np.meshgrid(x, y)
    QSO=gaussian_2d(X, Y, x_qso, y_qso, FWHM_qso/2.355, FWHM_qso/2.355)
    host=gaussian_2d(X, Y, x_host, y_host, FWHM_qso/2.355, FWHM_qso/2.355)

    print(x_qso)
    print(y_qso)
    print(X)
    print(Y)


#calculate the scale factor for two objects
    aa=np.where(np.abs(x-x_qso) <= FWHM_qso/2.355*3)
    bb=np.where(np.abs(y-y_qso) <= FWHM_qso/2.355*3)
    cc=np.sum(QSO[10:22,10:22])

    scale_host=Total_host/cc
    scale_qso=Total_qso/cc

#psf for two objects
    PSF_QSO=scale_qso*gaussian_2d(X, Y, x_qso, y_qso, FWHM_qso/2.355, FWHM_qso/2.355)
    PSF_host=scale_host*gaussian_2d(X, Y, x_host, y_host, FWHM_qso/2.355, FWHM_qso/2.355)

    ZZ=PSF_QSO+PSF_host+Ntotal_sky*t*QE*trans+DC*t+RN**2

    plt.imshow(ZZ, vmin=ZZ.min(), vmax=ZZ.max(), origin='lower',
           extent=[x.min(), x.max(), y.min(), y.max()])
    #cbar=plt.colorbar()
    plt.savefig("2d_simulation_EW15_"+str(impact_parameter[i]*10).zfill(2)+".png")
    ims=fits.PrimaryHDU(ZZ)
    opfnam='2d_simulation_EW15_%d.fit'  % i
    ims.writeto(opfnam, clobber=True)
