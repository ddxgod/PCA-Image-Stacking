# Plot the SB curve for each band



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



# Finds the mean of all photon counts where the value is 3 sigma above the mean

def photoncount(scidata, radius1, radius2):
	flux = 0
	length = 0
	#print(np.shape(scidata))
				
	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if distance(i, j, 50, 50) <= radius1 and distance(i, j, 50, 50) >= radius2:
				flux += scidata[i][j] / 2000 / 10**8
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


def photoncount2(scidata, radius1, radius2):
	flux = 0
	length = 0
	#print(np.shape(scidata))
				
	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if scidata[i][j] > 0 and distance(i, j, 500, 500) <= radius1 and distance(i, j, 500, 500) >= radius2:
				flux += scidata[i][j]
				length += 1

	return flux / length




def errorcalc(scidata, radius1, radius2):
	array = []

	for i in range(len(scidata)):
		for j in range(len(scidata[0])):
			if distance(i, j, 50, 50) <= radius1 and distance(i, j, 50, 50) >= radius2:
				array.append(scidata[i][j])

	#array -= np.mean(array)
	return np.std(array)
	#return np.sqrt(np.mean(np.array(array)**2))

#def getdist(z):
	
"""
dirs = os.listdir('/data/marvels/billzhu/stacked/')
writer = open('SBMagnitudes2.txt', 'w')

success = 0
fail = 0

for i in range(len(dirs)):
	filename = '/data/marvels/billzhu/stacked/' + dirs[i]
	print(filename)
	hdulist = fits.open(filename)
	scidata = hdulist[0].data.astype(float)
	SBarray = []
	scale = cosmo.kpc_proper_per_arcmin(float(dirs[i].split('_')[2])) * u.arcmin / u.kiloparsec * 0.396 / 60
	print(scale)
	
	for j in range(1, 7):
		count = photoncount(scidata, 10 * j / scale)
		if count == 0:
			continue
		f = count/(10**8 * 2000)
		mag = -2.5 * np.log10(f)
		print(mag)

		# Calculate the surface brightness by adding the amount of kpc, similar to traditional mg/arcsec^2 formula
		surface_brightness = mag + 2.5 * np.log10(math.pi/4 * (10 * j)**2)
		print(surface_brightness)
		SBarray.append(surface_brightness)

	writer.write(dirs[i] + '\t' + str(mag) + '\t' + str(SBarray) + '\n')
	success += 1



writer.write('\n\n\nSuccess: ' + str(success))
writer.write('\nFail: ' + str(fail))
writer.close()        
	# How to convert counts to mag in image where no reference stars are available
	# Need to convert reference stars to same redshift, then calculate instrumental magnitude
	# Need apparent magnitude of reference star as well
	# AstroImageJ?



	# Which flux20 to use in final image?
	# Rescaling seems ok based off of flux20 data

	# Bigger cutouts means masking array?
	# Their method is janky
"""


"""
x = np.array([0, 1, 2, 3])
y = np.array([-1, 0.2, 0.9, 2.1])
A = np.vstack([x, np.ones(len(x))]).T
m, c = np.linalg.lstsq(A, y)[0]
plt.plot(x, y, 'o', label='Original data', markersize=10)
plt.plot(x, m*x + c, 'r', label='Fitted line')
plt.legend()
plt.show()
"""


"""
hdulist = fits.open('TotComb13-z.fit')
scidata = hdulist[0].data.astype(float)
mean, median, stddev = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
#scidata -= median * 8
#spline = interpolate.interp2d(np.arange(len(scidata)), np.arange(len(scidata)), scidata)
#scidata = spline(np.arange(0, len(scidata), 0.1), np.arange(0, len(scidata), 0.1))
#scidata *= 1.15
SBarray = []
error_low = []
error_high = []
outter = []
scale = cosmo.kpc_proper_per_arcmin(0.48) * u.arcmin / u.kiloparsec * 0.396 / 60
print(scale)

fig, ax = plt.subplots()

for j in range(1, 9):
	#print(5 * j / scale)
	f = photoncount(scidata, 6 * (25 / 3)**((1.0 / 8)*j) / scale, 6 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale)
		
	print("%f, %f" % (6 * (25 / 3)**((1.0 / 8)*j) / scale, 6 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale))
	outter.append(6 * (25 / 3)**((1.0 / 8)*j))
	#f /= 50

	npix = getlength(scidata, 6 * (25 / 3)**((1.0 / 8)*j) / scale, 6 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale)
	f *= 2000 * 10**8# * npix
	print(f)

	npix = getlength(scidata, 6 * (25 / 3)**((1.0 / 8)*j) / scale, 6 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale)
	print(npix)
	#error = math.sqrt((f + 70) / 4.8) / 2000 / 10**8
	error = math.sqrt(f / npix) / 2000 / 10**8
	f /= 2000 * 10**8# * npix
	#error = errorcalc(scidata, 6 * (25 / 3)**((1.0 / 8)*j) / scale, 6 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale)
	
	#print("%f, %f" % (f, error))
	mag = -2.5 * np.log10(f)
	el_mag = -2.5 * np.log10(f - error) + 2.5 * np.log10(scale**2)
	eh_mag = -2.5 * np.log10(f + error) + 2.5 * np.log10(scale**2)
	print("%f, %f" % (el_mag, eh_mag))

	# Calculate the surface brightness by adding the amount of kpc, similar to traditional mg/arcsec^2 formula
	surface_brightness = mag + 2.5 * np.log10(scale**2)
	print(surface_brightness)
	SBarray.append(surface_brightness)
	error_low.append(el_mag)
	error_high.append(eh_mag)

print(SBarray)

# fit least-squares with an intercept
#print(np.array(outter))
#print(np.array(SBarray))
#w = np.linalg.lstsq(np.vstack([outter, np.ones(8)]).T, SBarray)[0]

# plot best-fit line
#plt.plot(outter, w[0]*np.array(outter) + w[1], 'r')


error = [np.array(SBarray) - np.array(error_high), np.array(error_low) - np.array(SBarray)]
#yerr = np.multiply(error, [[6, 6, 6, 6], [6, 6, 6, 6]])
yerr = error

print(error)

ax.errorbar(outter, SBarray, yerr = yerr, fmt='bo')
(_, caps, _) = ax.errorbar(outter, SBarray, yerr=yerr, capsize=5, elinewidth=2, fmt='bo')

for cap in caps:
	cap.set_color('black')
	cap.set_markeredgewidth(1.5)

plt.plot(np.unique(outter), np.poly1d(np.polyfit(outter, SBarray, 3))(np.unique(outter)), 'b')
#plt.plot(outter, SBarray, 'bo')
#plt.plot([3.8, 10, 13, 16, 19, 27, 37, 51], [31.6, 31.5, 31.4, 31.6, 32.1, 32.3, 33, 33.25], 'go')
#plt.plot([3.8, 10, 13, 16, 19, 27, 37, 51], [29.6, 30.3, 30.7, 31.2, 31.8, 32.1, 32.8, 33.25], 'yo')

#fit = powerlaw.Fit(SBarray)
#print(fit.power_law.alpha)
#print(fit.power_law.xmin)


"""
fig, ax = plt.subplots()

hdulist = fits.open('TotComb50-g.fit')
scidata = hdulist[0].data.astype(float)
error_low = []
error_high = []
mean, median, stddev = sigma_clipped_stats(scidata, sigma=3.0, iters=5)
print(median)
#scidata -= median * 24
#spline = interpolate.interp2d(np.arange(len(scidata)), np.arange(len(scidata)), scidata)
#scidata = spline(np.arange(0, len(scidata), 0.1), np.arange(0, len(scidata), 0.1))
#scidata *= 1.15

SBarray = []
outter = []
scale = cosmo.kpc_proper_per_arcmin(0.48) * u.arcmin / u.kiloparsec * 0.396 / 60
print(scale)
#print(cosmo.kpc_proper_per_arcmin(1.03))
#print(cosmo.kpc_proper_per_arcmin(5))

#print("%s" % ('%.4E' % Decimal(photoncount(scidata, 250, 200))))
#print(-2.5 * math.log10(photoncount(scidata, 250, 200)) + 2.5 * math.log10(scale**2))

#print(-2.5 * math.log10(2/(10**8 * 2000)) + 2.5 * np.log10(scale**2))
#print(scale)
boundaries = [3, 10, 13, 16, 19, 27, 37, 51, 67, 100, 140]

for j in range(len(boundaries) - 1):
	#print(5 * j / scale)
	#f = photoncount(scidata, 6 * (25 / 3)**((1.0 / 8)*j) / scale, 6 * (25 / 3)**((1.0 / 8)*(j - 1)) / scale)
	f = photoncount(scidata, boundaries[j + 1] / scale, boundaries[j] / scale)
	print("%f, %f" % (boundaries[j + 1], boundaries[j]))
	print(f)
		
	#outter.append(6 * (25 / 3)**((1.0 / 8)*j))
	outter.append(boundaries[j + 1] / scale * 0.396)
	#error = errorcalc(scidata, boundaries[j + 1] / scale, boundaries[j] / scale) / 2000 / 10**8 / np.sqrt(getlength(scidata, boundaries[j + 1] / scale, boundaries[j] / scale)) * 3
	error = np.sqrt(f * 2000 * 10**8) / 2000 / 10**8
	print(error)
	print(f - error)
	#f /= 50

	#mag = 22.5 - 2.5 * np.log10(f)
	mag = -2.5 * np.log10(f)
	if np.isnan(mag):
		mag = 38
	
	#print(mag)

	# Calculate the surface brightness by adding the amount of kpc, similar to traditional mg/arcsec^2 formula
	surface_brightness = mag + 2.5 * math.log10(0.396**2)#np.log10(math.pi * ((60 * (25 / 3)**((1.0 / 8) * j))**2 - (60 * (25 / 3)**((1.0 / 8) * (j - 1)))**2) / 50)


	el_mag = -2.5 * np.log10(f - error) + 2.5 * np.log10(0.396**2)
	eh_mag = -2.5 * np.log10(f + error) + 2.5 * np.log10(0.396**2)


	#if np.isnan(el_mag):
	#    el_mag = 40


	print(surface_brightness)
	SBarray.append(surface_brightness)
	error_low.append(el_mag)
	error_high.append(eh_mag)

	

print(SBarray)
print(error_high)


# plot best-fit line
#plt.plot(outter, w[0]*np.array(outter) + w[1], 'r')
#plt.plot(np.unique(outter), np.poly1d(np.polyfit(outter, SBarray, 3))(np.unique(outter)), 'b')
plt.plot(outter, SBarray, 'ro')
#plt.plot(outter, error_low, 'go')
#plt.plot(outter, error_high, 'go')




yerr = [np.array(SBarray) - np.array(error_high), np.array(error_low) - np.array(SBarray)]

ax.errorbar(outter, SBarray, fmt='bo')#, yerr = yerr)



(_, caps, _) = ax.errorbar(outter, SBarray, yerr=yerr, capsize=5, elinewidth=2, fmt='bo')

for cap in caps:
	cap.set_color('black')
	cap.set_markeredgewidth(1.5)


plt.axis([1, 50, 38, 24])

majorLocator = MultipleLocator(2)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(0.5)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)


#plt.axes().yaxis.set_tick_params(which='minor', right = 'off')
plt.xscale('log')
plt.xlabel('R [arc]', fontsize = 12)
plt.ylabel('Total SB [mag/arcsec$^2$]', fontsize = 12)
plt.title('Mg II $g$-band 0.37 $\leq z_{abs}$ < 0.55', fontsize=20)

#MG_patch = mpatches.Patch(color='blue', label='2DA SB Profile')
#NET_patch = mpatches.Patch(color='blue', label='Net SB Profile')
#plt.legend(handles=[MG_patch])
plt.show()

	
