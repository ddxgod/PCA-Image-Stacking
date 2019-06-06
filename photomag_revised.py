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

#pseudo code for SB profile and error calculation here:
def distance(x, y, x1, y1):
	return math.sqrt((x - x1)**2 + (y - y1)**2)



def photoncount(scidata, radius1, radius2):
	flux = 0
	length = 0
	#print(np.shape(scidata))

	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if distance(i, j, 50, 50) <= radius1 and distance(i, j, 50, 50) >= radius2:
				flux += scidata[i, j]
				length += 1

	return flux / length




# Find the number of points within 10 kpc

def getlength(scidata, radius1, radius2):
	length = 0
	mean, median, stddev = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if distance(i, j, 50, 50) <= radius1 and distance(i, j, 50, 50) >= radius2:
				length += 1

	return length



#load the data here:
hdulist = fits.open('MGSUBComb62-g.fit')
scidata = hdulist[0].data.astype(float)

scale = cosmo.kpc_proper_per_arcmin(0.48) * u.arcmin / u.kiloparsec * 0.396 / 60 #kpc/pixel

#get the size of image:
width = scidata.shape[0]
max_boundary = width/scale*0.5
print(max_boundary)

fig, ax = plt.subplots()


boundaries = [3, 10, 13, 16, 19, 27, 37, 51, 67, 100, 140] #kpc


#calcluate the sky error here:
sky_count = photoncount(scidata, boundaries[-1] / scale,  max_boundary/ scale)
sky_flux = sky_count / 2000 / 10**8

SBarray = []
error_array = []
outter = []
for j in range(len(boundaries) - 1):

	f_count = photoncount(scidata, boundaries[j + 1] / scale, boundaries[j] / scale)
	print("%f, %f" % (boundaries[j + 1], boundaries[j]))
	print(f_count)

	if f_count < 0:
		continue

	#outter.append(6 * (25 / 3)**((1.0 / 8)*j))
	outter.append(boundaries[j + 1] / scale * 0.396) # the unit is arcsec
	f_flux = f_count / 2000 / 10**8 #covert the count
	print(f_flux)

	error = np.sqrt(f_flux) / np.sqrt(2000 * 10**8)


	#flux per pxiel**2 to surface_brightness per arcsec**2
	mag = -2.5 * np.log10(f_flux/0.396**2)

	# start to calculate the error by assuming is the gaussain error.
	mean_f = f_flux
	sigma_f = error
	print(error)
	f_array = np.random.normal(mean_f, sigma_f, 1000)
	mag_array = -2.5 * np.log10(f_array/0.396**2)
	error_f = sigma_clipped_stats(f_array, mask_value=float('nan'))[2]
	error_mag = (-2.5 * np.log10((mean_f + error_f) / 0.396**2) - mag) 
	print(error_f)
	print(mag)
	print(error_mag)
	SBarray.append(mag)
	error_array.append(error_mag)


plt.plot(np.unique(outter), np.poly1d(np.polyfit(outter, SBarray, 3))(np.unique(outter)), 'r')
plt.plot(outter, SBarray, 'ro')

"""
yerr = [np.array(error_array), np.array(error_array)]

ax.errorbar(outter, SBarray, fmt='ro')#, yerr = yerr)

(_, caps, _) = ax.errorbar(outter, SBarray, yerr=yerr, capsize=5, elinewidth=2, fmt='ro')

for cap in caps:
	cap.set_color('black')
	cap.set_markeredgewidth(1.5)
"""




#load the data here:
hdulist = fits.open('TotComb62-g.fit')
scidata = hdulist[0].data.astype(float)

scale = cosmo.kpc_proper_per_arcmin(0.48) * u.arcmin / u.kiloparsec * 0.396 / 60 #kpc/pixel

#get the size of image:
width = scidata.shape[0]
max_boundary = width/scale*0.5
print(max_boundary)



#fig, ax = plt.subplots()

boundaries = [3, 10, 13, 16, 19, 27, 37, 51, 67, 100, 140] #kpc


#calcluate the sky error here:
sky_count = photoncount(scidata, boundaries[-1] / scale,  max_boundary/ scale)
sky_flux = sky_count# / 2000 / 10**8

SBarray = []
error_array = []
outter = []
for j in range(len(boundaries) - 1):

	f_count = photoncount(scidata, boundaries[j + 1] / scale, boundaries[j] / scale)
	print("%f, %f" % (boundaries[j + 1], boundaries[j]))
	print(f_count)

	f_flux = 0
	error = 0


	if f_count < 0:
		"""
		mag = 36
		f_flux = np.abs(f_count)
		error = np.sqrt(f_flux)
		outter.append(boundaries[j + 1] / scale * 0.396)
		"""

	else:

		#outter.append(6 * (25 / 3)**((1.0 / 8)*j))
		outter.append(boundaries[j + 1] / scale * 0.396) # the unit is arcsec
		f_flux = f_count / 2000 / 10**8 #covert the count
		print(f_flux)

		error = np.sqrt(f_count / 2000 * 2479/ 4.77) / (2479 * 10**8)
		#error /= np.sqrt(883)


		#flux per pxiel**2 to surface_brightness per arcsec**2
		mag = -2.5 * np.log10(f_flux/0.396**2)

		# start to calculate the error by assuming is the gaussain error.
		mean_f = f_flux
		sigma_f = error
		print(error)
		f_array = np.random.normal(mean_f, sigma_f, 1000)
		mag_array = -2.5 * np.log10(f_array/0.396**2)
		error_f = sigma_clipped_stats(f_array, mask_value=float('nan'))[2]# / 4
		error_mag = (-2.5 * np.log10((mean_f + error_f) / 0.396**2) - mag)
		print(error_f)
		print(mag)
		print(error_mag)
		SBarray.append(mag)
		error_array.append(error_mag)






plt.plot(np.unique(outter), np.poly1d(np.polyfit(outter, SBarray, 3))(np.unique(outter)), 'b')


yerr = [np.array(error_array), np.array(error_array)]

ax.errorbar(outter, SBarray, fmt='bo')#, yerr = yerr)



(_, caps, _) = ax.errorbar(outter, SBarray, yerr=yerr, capsize=5, elinewidth=2, fmt='bo')

for cap in caps:
	cap.set_color('black')
	cap.set_markeredgewidth(1.5)


plt.axis([1, 30, 34, 24])

majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(0.5)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)


#plt.axes().yaxis.set_tick_params(which='minor', right = 'off')
plt.xscale('log')
plt.xlabel('R [arcsec]', fontsize = 12)
plt.ylabel('Total SB [mag/arcsec$^2$]', fontsize = 12)
plt.title('Mg II z-band 0.37 $\leq z_{abs}$ < 0.55', fontsize=20)

#MG_patch = mpatches.Patch(color='blue', label='2DA SB Profile')
#NET_patch = mpatches.Patch(color='blue', label='Net SB Profile')
#plt.legend(handles=[MG_patch])
plt.show()
