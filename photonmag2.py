# Create the SED


from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
import numpy as np
import math
import string
import os
import linecache
from scipy import interpolate
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from decimal import Decimal
from scipy.constants import c, pi



# Calculates the distance btween two given points

def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)


def inbounds(x, y):
	return x < 2048 and x >= 0 and y < 1489 and y >= 0


# Finds the mean of all photon counts where the value is 3 sigma above the mean

def photoncount(scidata, radius1, radius2):
	flux = 0
	length = 0
	#print(np.shape(scidata))
				
	mean, median, std = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
	
	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if scidata[i][j] > 0 and distance(i, j, 500, 500) <= radius1 and distance(i, j, 500, 500) >= radius2:
				flux += scidata[i][j]

	return flux



# Calculate the background by taking the sigma clipped median value between the radii specified

def calc_background(data, x, y, radius1, radius2):
	bkg_array = []
	for i in range(len(data[0])):
		for j in range(len(data)):
			if abs(x - i) < radius2 and abs(y - j) < radius2 and inbounds(i, j) and (i - x)**2 + (j - y)**2 >= radius1**2 and (i - x)**2 + (j - y)**2 <= radius2**2:
				bkg_array.append(data[j, i])

	true_bkg = sigma_clipped_stats(bkg_array, sigma=3.0, iters=5)[0]
	print(true_bkg)
	return true_bkg



def sedfitter(flux, wavelength):
	flux = 10**(flux / -2.5)
	return flux * 3631 * 10**(-23) * 29979245800/(wavelength * 10**-8)**2 * 10**-8


def errorcalc(filename):
	band = fits.open(filename)
	scidata = band[0].data.astype(float)
	outside = []
	count = 0
	total_flux = 0
	
	for i in range(len(scidata)):
		for j in range(len(scidata)):
			if distance(i, j, 50, 50) >= 50 and distance(i, j, 50, 50) <= 60:
				outside.append(scidata[i][j])
				count += 1

	mean, median, stddev = sigma_clipped_stats(outside, sigma=3.0, iters=5)
	return stddev * math.sqrt(count)



def lambda_flambda_to_fnu(wavelength, flambda):
    """
    Convert a Fλ vs λ spectrum to Fν vs λ

    Parameters
    ----------
    wavelength: list-like of floats
        The wavelengths in nm.
    flambda: list-like of floats
        Fλ flux density in W/m²/nm (or Lλ luminosity density in W/nm).

    Returns
    -------
    fnu: array of floats
        The Fν flux density in mJy (or the Lν luminosity density in
        1.e-29 W/Hz).

    """
    wavelength = np.array(wavelength, dtype=float)
    flambda = np.array(flambda, dtype=float)

    # Factor 1e+29 is to switch from W/m²/Hz to mJy
    # Factor 1e-9 is to switch from nm to m (only one because the other nm
    # wavelength goes with the Fλ in W/m²/nm).
    fnu = 1e+29 * 1e-9 * flambda * wavelength * wavelength / c

    return fnu




total_flux = []

hdulist = fits.open('TotComb62-g.fit')
scidata = hdulist[0].data.astype(float)
#scidata -= calc_background(scidata, 50, 50, 65, 70)
nmgy_count = 0

for i in range(len(scidata)):
	for j in range(len(scidata[0])):
		if distance(i, j, 50, 50) > 4 and distance(i, j, 50, 50) <= 42:
			nmgy_count += scidata[i][j] / 2000 / 10**8

gmag = -2.5 * math.log10(nmgy_count)# - 1.12
print(gmag)
total_flux.append(10**(gmag/-2.5) * 3631 * 10**-23 * 29979245800/(4770 * 10**-8)**2 * 10**-8)
print(10**(gmag/-2.5) * 3631 * 10**3)


hdulist = fits.open('TotComb62-r.fit')
scidata = hdulist[0].data.astype(float)
#scidata -= calc_background(scidata, 50, 50, 65, 70)
nmgy_count = 0

for i in range(len(scidata)):
	for j in range(len(scidata[0])):
		if distance(i, j, 50, 50) > 4 and distance(i, j, 50, 50) <= 42:
			nmgy_count += scidata[i][j] / 2000 / 10**8

rmag = -2.5 * math.log10(nmgy_count)# - 0.49
print(rmag)
total_flux.append(10**(rmag/-2.5) * 3631 * 10**-23 * 29979245800/(6231 * 10**-8)**2 * 10**-8)
print(10**(rmag/-2.5) * 3631 * 10**3)

hdulist = fits.open('TotComb62-i.fit')
scidata = hdulist[0].data.astype(float)
#scidata -= calc_background(scidata, 50, 50, 65, 70)
nmgy_count = 0

for i in range(len(scidata)):
	for j in range(len(scidata[0])):
		if distance(i, j, 50, 50) > 4 and distance(i, j, 50, 50) <= 42:
			nmgy_count += scidata[i][j] / 2000 / 10**8

imag = -2.5 * math.log10(nmgy_count)# - 0.08
print(imag)
total_flux.append(10**(imag/-2.5) * 3631 * 10**-23 * 29979245800/(7625 * 10**-8)**2 * 10**-8)
print(10**(imag/-2.5) * 3631 * 10**3)

hdulist = fits.open('TotComb62-z.fit')
scidata = hdulist[0].data.astype(float)
#scidata -= calc_background(scidata, 50, 50, 65, 70)
nmgy_count = 0

for i in range(len(scidata)):
	for j in range(len(scidata[0])):
		if distance(i, j, 50, 50) > 4 and distance(i, j, 50, 50) <= 42:
			nmgy_count += scidata[i][j] / 2000 / 10**8

zmag = -2.5 * math.log10(nmgy_count)# - 0.08
print(zmag)
total_flux.append(10**(zmag/-2.5) * 3631 * 10**-23 * 29979245800/(9134 * 10**-8)**2 * 10**-8)
print(10**(zmag/-2.5) * 3631 * 10**3)


gflux = gmag# - 1.12
rflux = rmag# - 0.49
iflux = imag# - 0.08
zflux = zmag# - 0.08

gsed = total_flux[0]
rsed = total_flux[1]
ised = total_flux[2]
zsed = total_flux[3]

"""
gerror_low = -2.5 * math.log10(10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-g.fit')) + 22.5) / -2.5) + 10**(gflux / -2.5))
rerror_low = -2.5 * math.log10(10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-r.fit')) + 22.5) / -2.5) + 10**(rflux / -2.5))
ierror_low = -2.5 * math.log10(10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-i.fit')) + 22.5) / -2.5) + 10**(iflux / -2.5))
zerror_low = -2.5 * math.log10(10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-z.fit')) + 22.5) / -2.5) + 10**(zflux / -2.5))
gerror_high = -2.5 * math.log10(-10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-g.fit')) + 22.5) / -2.5) + 10**(gflux / -2.5))
rerror_high = -2.5 * math.log10(-10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-r.fit')) + 22.5) / -2.5) + 10**(rflux / -2.5))
ierror_high = -2.5 * math.log10(-10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-i.fit')) + 22.5) / -2.5) + 10**(iflux / -2.5))
zerror_high = -2.5 * math.log10(-10**((-2.5 * math.log10(errorcalc('2175REFSUBComb13-z.fit')) + 22.5) / -2.5) + 10**(zflux / -2.5))

print("%s: %f, %f, %f, %s" % ('g-band', gflux, gerror_low, gerror_high, '%.4E' % Decimal(gsed)))
print("%s: %f, %f, %f, %s" % ('r-band', rflux, rerror_low, rerror_high, '%.4E' % Decimal(rsed)))
print("%s: %f, %f, %f, %s" % ('i-band', iflux, ierror_low, ierror_high, '%.4E' % Decimal(ised)))
print("%s: %f, %f, %f, %s" % ('z-band', zflux, zerror_low, zerror_high, '%.4E' % Decimal(zsed)))
"""



"""
4770
6231
7625
9134
"""


fig, ax = plt.subplots(1, 2)

gerror_low = -2.5 * math.log10((10**(gflux/-2.5) * 2478 * 10**8 + 3 * math.sqrt(10**(gflux/-2.5) * 2478 * 10**8)) / 2478 / 10**8)
rerror_low = -2.5 * math.log10((10**(rflux/-2.5) * 1952 * 10**8 + 3 * math.sqrt(10**(rflux/-2.5) * 1952 * 10**8)) / 1952 / 10**8)
ierror_low = -2.5 * math.log10((10**(iflux/-2.5) * 1439 * 10**8 + 3 * math.sqrt(10**(iflux/-2.5) * 1439 * 10**8)) / 1439 / 10**8)
zerror_low = -2.5 * math.log10((10**(zflux/-2.5) * 301 * 10**8 + 3 * math.sqrt(10**(zflux/-2.5) * 301 * 10**8)) / 301 / 10**8)
gerror_high = -2.5 * math.log10((10**(gflux/-2.5) * 2478 * 10**8 - 3 * math.sqrt(10**(gflux/-2.5) * 2478 * 10**8)) / 2478 / 10**8)
rerror_high = -2.5 * math.log10((10**(rflux/-2.5) * 1952 * 10**8 - 3 * math.sqrt(10**(rflux/-2.5) * 1952 * 10**8)) / 1952 / 10**8)
ierror_high = -2.5 * math.log10((10**(iflux/-2.5) * 1439 * 10**8 - 3 * math.sqrt(10**(iflux/-2.5) * 1439 * 10**8)) / 1439 / 10**8)
zerror_high = -2.5 * math.log10((10**(zflux/-2.5) * 301 * 10**8 - 3 * math.sqrt(10**(zflux/-2.5) * 301 * 10**8)) / 301 / 10**8)

print("%s: %f, %f, %f, %s" % ('g-band', gflux, gerror_low, gerror_high, '%.4E' % Decimal(gsed)))
print("%s: %f, %f, %f, %s" % ('r-band', rflux, rerror_low, rerror_high, '%.4E' % Decimal(rsed)))
print("%s: %f, %f, %f, %s" % ('i-band', iflux, ierror_low, ierror_high, '%.4E' % Decimal(ised)))
print("%s: %f, %f, %f, %s" % ('z-band', zflux, zerror_low, zerror_high, '%.4E' % Decimal(zsed)))

print(10**(gmag/-2.5) * 3631 * 10**3)
print(10**(gerror_low/-2.5) * 3631 * 10**3 - 10**(gmag/-2.5) * 3631 * 10**3)
print(10**(rmag/-2.5) * 3631 * 10**3)
print(10**(rerror_low/-2.5) * 3631 * 10**3 - 10**(rmag/-2.5) * 3631 * 10**3)
print(10**(imag/-2.5) * 3631 * 10**3)
print(10**(ierror_low/-2.5) * 3631 * 10**3 - 10**(imag/-2.5) * 3631 * 10**3)
print(10**(zmag/-2.5) * 3631 * 10**3)
print(10**(zerror_low/-2.5) * 3631 * 10**3 - 10**(zmag/-2.5) * 3631 * 10**3)


error = [[sedfitter(gerror_low, 4770) - gsed, sedfitter(rerror_low, 6231) - rsed, sedfitter(ierror_low, 7625) - ised, sedfitter(zerror_low, 9134) - zsed], [gsed - sedfitter(gerror_high, 4770), rsed - sedfitter(rerror_high, 6231), ised - sedfitter(ierror_high, 7625), zsed - sedfitter(zerror_high, 9134)]]
#yerr = np.multiply(error, [[6, 6, 6, 6], [6, 6, 6, 6]])
yerr = error

#ax.errorbar([4770, 6231, 7625, 9134], [gsed, rsed, ised, zsed], xerr = xerr, fmt='bo')

xerr = [[770, 731, 775, 634], [730, 619, 875, 616]]
(_, caps, _) = ax[0].errorbar([4770, 6231, 7625, 9134], [gsed, rsed, ised, zsed], xerr = xerr, yerr = yerr, capsize=5, elinewidth=2, fmt='bo')

for cap in caps:
	cap.set_color('blue')
	cap.set_markeredgewidth(1.5)


error2 = [[sedfitter(22.28, 4770) - sedfitter(22.49, 4770), sedfitter(21.14, 6231) - sedfitter(21.25, 6231), sedfitter(20.65, 7625) - sedfitter(20.77, 7625), sedfitter(20.03, 9134) - sedfitter(20.23, 9134)], [sedfitter(22.49, 4770) - sedfitter(22.76, 4770), sedfitter(21.25, 6231) - sedfitter(21.37, 6231), sedfitter(20.77, 7625) - sedfitter(20.91, 7625), sedfitter(20.23, 9134) - sedfitter(20.47, 9134)]]


(_, caps, _) = ax[0].errorbar([4770, 6231, 7625, 9134], [sedfitter(22.49, 4770), sedfitter(21.25, 6231), sedfitter(20.77, 7625), sedfitter(20.23, 9134)], yerr = error2, capsize=5, elinewidth=2, fmt='ro')

for cap in caps:
	cap.set_color('red')
	cap.set_markeredgewidth(1.5)

print(error)

majorLocator = MultipleLocator(5 * 10**-18)
#majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(1 * 10**-18)

ax[0].yaxis.set_major_locator(majorLocator)
#ax.yaxis.set_major_formatter(majorFormatter)
ax[0].yaxis.set_minor_locator(minorLocator)

#ax.errorbar([4770, 6231, 7625, 9134], [gsed, rsed, ised, zsed], yerr = yerr, fmt='bo')
ax[0].axis([3000, 10000, 0, 1.6 * 10**-17])
ax[0].set_xlabel('Wavelength [$\AA$]', fontsize=12)
ax[0].set_ylabel('Flux Density [erg sec$^{-1}$ cm$^{-2}$ $\AA^{-1}$]', fontsize=12)
#fig.suptitle('SED of Mg II Absorbers, 0.37 $\leq z_{abs}$ < 0.55', fontsize=20)#, fontname='Times New Roman')
MG_patch = mpatches.Patch(color='red', label='Zibetti et al. (2007) SED')
MY_patch = mpatches.Patch(color='blue', label='Our Mg II $z_{abs}$ = 0.37 - 0.55 SED')
ax[0].legend(handles=[MG_patch, MY_patch])






table = Table.read('MGIIAbsorbers_best_model.fits', hdu=1)
#fluxdata = file[0].data

#for i in range(len(fluxdata)):
#	fluxdata[i] *= 


ax[1].plot(10 * table['wavelength'][400 : 1120], table['Fnu'][400 : 1120], 'black')

#nebular_continuum_old = lambda_flambda_to_fnu(table['wavelength'], table['nebular.lines_old'])
#plt.plot(10 * table['wavelength'][237 : 1115], nebular_continuum_old[237 : 1115])
#plt.plot([486, 628, 770, 913], [0.01028078, 0.01721973, 0.0178659554, 0.0328115124], 'ro')

xerr = [[770, 731, 775, 634], [730, 619, 875, 616]]
error2 = [[0.0007596202484747563, 0.0013467644293272667, 0.001973498688289742, 0.005766794212355562], [0.0007596202484747563, 0.0013467644293272667, 0.001973498688289742, 0.005766794212355562]]



# K-corrected: [0.01028078, 0.01721973, 0.0178659554, 0.0328115124], [0.00141492, 0.00166133, 0.00172368, 0.0055615314]

(_, caps, _) = ax[1].errorbar([4770, 6231, 7625, 9134], [0.004414325116227801, 0.01054552285379697, 0.01668525934622801, 0.029511966374186765], xerr = xerr, yerr = error2, capsize=5, elinewidth=2, fmt='bo')

for cap in caps:
    cap.set_color('blue')
    cap.set_markeredgewidth(1.5)




#plt.axes().yaxis.set_tick_params(which='minor', right = 'off')
ax[1].set_xscale('log')
ax[1].set_xlabel('Wavelength [$\AA$]', fontsize = 12)
ax[1].set_ylabel('Flux Density [mJy]', fontsize = 12)
#ax[1].set_title('Mg II Absorber Host Galaxy SED 0.37 $\leq$ $z_{abs}$ < 0.55', fontsize=20, y=1.04)

MG_patch = mpatches.Patch(color='blue', label='SDSS Photometric Fluxes')
NET_patch = mpatches.Patch(color='black', label='CIGALE SED Fitting')
ax[1].legend(handles=[MG_patch, NET_patch])


plt.show()



"""
spline = interpolate.interp2d(np.arange(len(scidata)), np.arange(len(scidata)), scidata)
scidata = spline(np.arange(0, len(scidata), 0.1), np.arange(0, len(scidata), 0.1))
#scidata *= 1.15
SBarray = []
outter = []
scale = cosmo.kpc_proper_per_arcmin(1.00) * u.arcmin / u.kiloparsec * 0.396 / 60
print(scale)

for j in range(1, 9):
	#print(5 * j / scale)
	f = photoncount(scidata, 60 * (25 / 3)**((1.0 / 8)*j) / scale, 60 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale)
		
	print("%f, %f" % (60 * (25 / 3)**((1.0 / 8)*j) / scale, 60 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale))
	outter.append(6 * (25 / 3)**((1.0 / 8)*j))
	f /= 100
	mag = -2.5 * np.log10(f) + 22.5
	#print(mag)

	# Calculate the surface brightness by adding the amount of kpc, similar to traditional mg/arcsec^2 formula
	surface_brightness = mag + 2.5 * np.log10(math.pi * ((60 * (25 / 3)**((1.0 / 8) * j))**2 - (60 * (25 / 3)**((1.0 / 8) * (j - 1)))**2) / 100)
	print(surface_brightness)
	SBarray.append(surface_brightness)

print(SBarray)


plt.plot(np.unique(outter), np.poly1d(np.polyfit(outter, SBarray, 3))(np.unique(outter)))
plt.plot(outter, SBarray, 'ro')

plt.axis([6, 60, 38, 28])
plt.xscale('log')
plt.xlabel('R [kpc]')
plt.ylabel('Total SB [mag/kpc^2]')
plt.show()
"""
