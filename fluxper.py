# Calculate and plot the percentage of flux at different impact parameters




from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import scipy as sp
import os
import linecache
import math
from scipy import interpolate
from scipy.stats import powerlaw
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from decimal import Decimal
import powerlaw


# Calculates the distance btween two given points

def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)




def gaussian_2d(x, y, x0, y0, xsig, ysig):
	return np.exp(-0.5*(((x-x0) / xsig)**2 + ((y-y0) / ysig)**2))


# Method that calculates the total photon count within 3 sigma of the quasar centroid

def photonCount(xc, yc, sigma, data):
	count = 0
	for i in range(len(data)):
		for j in range(len(data)):
			if distance(i, j, xc, yc) <= sigma:
				count += data[i, j]
	return count



# Create the individual PSF frames

def calc_galflux(fwhm, comb_flux):
	FWHM = fwhm / 0.396
	plate_scale = 1
	
	x_qso = 50
	y_qso = 50
	x = np.arange(0, 101, plate_scale)
	y = np.arange(0, 101, plate_scale)
	X, Y = np.meshgrid(x, y)
	QSO = gaussian_2d(X, Y, x_qso, y_qso, FWHM/2.355, FWHM/2.355)
	QSO *= comb_flux / sum(sum(QSO))
	#plt.imshow(QSO + make_noise_image((101, 101), type='poisson', mean=0.001))
	

	return photonCount(50, 50, FWHM, QSO)






majorLocator = MultipleLocator(0.4)
majorFormatter = FormatStrFormatter('%.1f')
minorLocator = MultipleLocator(10)
minorFormatter = FormatStrFormatter('%f')


fig, ax = plt.subplots(1, 3)
hdulist = fits.open('SimulatedPSFSUBComb14-r5.fit')
data = hdulist[0].data.astype(float)
#data -= sigma_clipped_stats(data[0], sigma=3.0, iters=5)[0]

galflux = 10**(-1.2/2.5)
#galflux = calc_galflux(1.4, 10**(-1.2/2.5))
print(calc_galflux(1.4, 10**(-1.2/2.5)) / galflux)
totalflux = []
yerr = []
signal2 = []


#for i in range(len(data)):
#	data[i] -= sigma_clipped_stats(data[i], sigma=3.0, iters=5)[1]

for i in range(len(data)):
	#totalflux = sum(sum(data[i]))
	total_flux = 0
	#data[i] -= sigma_clipped_stats(data[i], sigma=3.0, iters=5)[1]
	signal = []

	for r in range(len(data[i])):
		for c in range(len(data[i][0])):
			if distance(r, c, 50, 50 + i) <= 1.4 / 0.396 * 1:		
				total_flux += data[i, r, c]

	counter2 = 0
	for r in range(len(data[0])):
		for c in range(len(data[0][0])):
			if distance(r, c, 75, 90) <= 1.4 / 0.396 * 1:
				signal.append(data[0, r, c])
				counter2 += 1

	#total_flux
	#print(counter2)
	#mean, median, std = sigma_clipped_stats(data[0])
	#yerr.append(3 * std / total_flux * 100)
	totalflux.append(total_flux)
	#signal2.append(3 * std)
	#signal -= np.mean(signal)
	yerr.append((sum(signal) + 3 * np.sqrt(np.mean(signal)**2)) / total_flux * 100)
	signal2.append((sum(signal) + 3 * np.sqrt(np.mean(signal)**2)))

	#total_flux /= len(data[0])**2

	"""
	fluxarr = []
	flux_radii = 0
	for r in range(15):
		flux_radii = photoncount(data[i], r + 1, r) / galflux
		fluxarr.append(flux_radii)
	
	"""

"""
	radii_range = np.arange(15)

	if i == 0:
		plt.plot(radii_range, fluxarr, 'o-', color='red')
	if i == 1:
		plt.plot(radii_range, fluxarr, 'o-', color='orange')
	if i == 2:
		plt.plot(radii_range, fluxarr, 'o-', color='gold')
	if i == 3:
		plt.plot(radii_range, fluxarr, 'o-', color='lime')
	if i == 4:
		plt.plot(radii_range, fluxarr, 'o-', color='green')
	if i == 5:
		plt.plot(radii_range, fluxarr, 'o-', color='navy')
	if i == 6:
		plt.plot(radii_range, fluxarr, 'o-', color='blue')
	if i == 7:
		plt.plot(radii_range, fluxarr, 'o-', color='magenta')
	



red_patch = mpatches.Patch(color='red', label='0 arcsec')
orange_patch = mpatches.Patch(color='orange', label='0.4 arcsec')
gold_patch = mpatches.Patch(color='gold', label='0.8 arcsec')
lime_patch = mpatches.Patch(color='lime', label='1.2 arcsec')
green_patch = mpatches.Patch(color='green', label='1.6 arcsec')
navy_patch = mpatches.Patch(color='navy', label='2.0 arcsec')
blue_patch = mpatches.Patch(color='blue', label='2.4 arcsec')
magenta_patch = mpatches.Patch(color='magenta', label='2.8 arcsec')

plt.legend(handles=[red_patch, orange_patch, gold_patch, lime_patch, green_patch, navy_patch, blue_patch, magenta_patch])
plt.xlabel('Radius [pixel ~ 0.396 arcsec]', fontsize=12)
plt.ylabel('Average pixel intensity [% of total]', fontsize=12)
plt.title('Average pixel intensity at increasing radii', y=1.05, fontsize=18)
plt.show()
"""

#totalflux /= totalflux[7]
#print(10**(-1.2/2.5))
print(totalflux)
ax[0].plot(np.arange(0.0, 3.2, 0.4), 100 * np.array(totalflux) / galflux, 'bo-')
ax[0].set_xlabel('Impact Parameter [arcsec]', fontsize=16)
ax[1].set_xlabel('Impact Parameter [arcsec]', fontsize=16)
ax[2].set_xlabel('Impact Parameter [arcsec]', fontsize=16)
ax[0].set_ylabel('Total reside flux ratio [% of simulated absorber flux]', fontsize=16)
ax[0].set_title('Seeing FWHM = 1.4 arcsec', y=1.03, fontsize=20)
ax[0].xaxis.set_major_locator(majorLocator)
ax[0].xaxis.set_major_formatter(majorFormatter)
ax[0].yaxis.set_minor_locator(minorLocator)
#ax[0].xaxis.set_minor_locator(minotFormatter)

#yerr = [[20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621], [20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621]]
print(yerr)


ax[0].axhline(y=100 * totalflux[1] / galflux + yerr[1], color='red', linestyle='--')
ax[0].axhline(y=100 * totalflux[2] / galflux + yerr[2], color='red', linestyle='--')
ax[0].set_ylim(-10, 110)


yerr = [yerr, yerr]

(_, caps, _) = ax[0].errorbar(np.arange(0.0, 3.2, 0.4), 100 * np.array(totalflux) / galflux, yerr=yerr, capsize=5, elinewidth=2, fmt='bo-')

for cap in caps:
    cap.set_color('black')
    cap.set_markeredgewidth(1.5)






hdulist = fits.open('SimulatedPSFSUBComb10-r5.fit')
data = hdulist[0].data.astype(float)
#data -= sigma_clipped_stats(data[0], sigma=3.0, iters=5)[0]
#galflux = calc_galflux(1.0, 10**(-1.2/2.5))
galflux = 10**(-1.2/2.5)
print(calc_galflux(1.0, 10**(-1.2/2.5)) / galflux)
totalflux = []
yerr = []

#for i in range(len(data)):
#	data[i] -= sigma_clipped_stats(data[i], sigma=3.0, iters=5)[1]

for i in range(len(data)):
	#totalflux = sum(sum(data[i]))
	total_flux = 0
	#data[i] -= sigma_clipped_stats(data[i], sigma=3.0, iters=5)[1]
	signal = []

	for r in range(len(data[i])):
		for c in range(len(data[i][0])):
			if distance(r, c, 50, 50 + i) <= 1.0 / 0.396 * 1:				
				total_flux += data[i, r, c]

	"""
	for r in range(len(data[i])):
		for c in range(len(data[i][0])):
			if distance(r, c, 50, 25) <= 1.0 / 0.396 * 4:
				signal.append(data[i, r, c])
	"""
	#total_flux
	totalflux.append(total_flux)
	yerr.append(signal2[i] / total_flux * 100)

	#total_flux /= len(data[0])**2

	"""
	fluxarr = []
	flux_radii = 0
	for r in range(15):
		flux_radii = photoncount(data[i], r + 1, r) / galflux
		fluxarr.append(flux_radii)
	"""


#totalflux /= totalflux[7]
print(totalflux)
#ax[1].plot(np.arange(0.0, 3.2, 0.4), 100 * np.array(totalflux) / galflux, 'bo-')
ax[1].set_title('Seeing FWHM = 1.0 arcsec', y=1.03, fontsize=20)
ax[1].xaxis.set_major_locator(majorLocator)
ax[1].xaxis.set_major_formatter(majorFormatter)
ax[1].yaxis.set_minor_locator(minorLocator)
ax[1].set_ybound(0, 100)
#ax[1].set_xmargin(0.005)
#ax[1].set_ymargin(0.005)
ax[1].autoscale_view()

#yerr = [[20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621], [20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621]]

ax[1].axhline(y=100 * totalflux[1] / galflux + yerr[1], color='red', linestyle='--')
ax[1].axhline(y=100 * totalflux[2] / galflux + yerr[2], color='red', linestyle='--')
ax[1].set_ylim(-10, 110)

yerr = [yerr, yerr]

(_, caps, _) = ax[1].errorbar(np.arange(0.0, 3.2, 0.4), 100 * np.array(totalflux) / galflux, yerr=yerr, capsize=5, elinewidth=2, fmt='bo-')

for cap in caps:
    cap.set_color('black')
    cap.set_markeredgewidth(1.5)








hdulist = fits.open('SimulatedPSFSUBComb06-r5.fit')
data = hdulist[0].data.astype(float)
#data -= sigma_clipped_stats(data[0], sigma=3.0, iters=5)[0]
galflux = 10**(-1.2/2.5)
#galflux = calc_galflux(0.6, 10**(-1.2/2.5))
print(calc_galflux(0.6, 10**(-1.2/2.5)) / galflux)
totalflux = []
yerr = []

#for i in range(len(data)):
#	data[i] -= sigma_clipped_stats(data[i], sigma=3.0, iters=5)[1]

for i in range(len(data)):
	#totalflux = sum(sum(data[i]))
	total_flux = 0
	#data[i] -= sigma_clipped_stats(data[i], sigma=3.0, iters=5)[1]
	#signal = []

	for r in range(len(data[i])):
		for c in range(len(data[i][0])):
			if distance(r, c, 50, 50 + i) <= 0.6 / 0.396 * 1:
				total_flux += data[i, r, c]

	"""
	for r in range(len(data[i])):
		for c in range(len(data[i][0])):
			if distance(r, c, 50, 25) <= 0.6 / 0.396 * 4:
				signal.append(data[i, r, c])
	"""
	
	totalflux.append(total_flux)
	yerr.append(signal2[i] / total_flux * 100)
	#total_flux /= len(data[0])**2

	"""
	fluxarr = []
	flux_radii = 0
	for r in range(15):
		flux_radii = photoncount(data[i], r + 1, r) / galflux
		fluxarr.append(flux_radii)
	"""



#totalflux /= totalflux[7]
#print(10**(1/2.5))
print(totalflux)
#ax[2].plot(np.arange(0.0, 3.2, 0.4), 100 * np.array(totalflux) / galflux, 'bo-')
ax[2].set_title('Seeing FWHM = 0.6 arcsec', y=1.03, fontsize=20)
ax[2].xaxis.set_major_locator(majorLocator)
ax[2].xaxis.set_major_formatter(majorFormatter)
ax[2].yaxis.set_minor_locator(minorLocator)

#yerr = [[20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621], [20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621, 20.5621]]

ax[2].axhline(y=100 * totalflux[1] / galflux + yerr[1], color='red', linestyle='--')
ax[2].axhline(y=100 * totalflux[2] / galflux + yerr[2], color='red', linestyle='--')
ax[2].set_ylim(-10, 110)


yerr = [yerr, yerr]


(_, caps, _) = ax[2].errorbar(np.arange(0.0, 3.2, 0.4), 100 * np.array(totalflux) / galflux, yerr=yerr, capsize=5, elinewidth=2, fmt='bo-')

for cap in caps:
    cap.set_color('black')
    cap.set_markeredgewidth(1.5)


plt.show()

print(100 * np.array(totalflux) / galflux)

"""
circularflux = 0
for r in range(1, 9):
	circularflux = photoncount(data[i], r, r - 1)
	print(circularflux)
"""

