from astropy.io import fits
from astropy.table import Table
import numpy as np
import io
import math
from decimal import Decimal
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from astropy.stats import sigma_clipped_stats
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

def distance(x1, y1, x2, y2):
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)


def integrated_flux(filename):
    band = fits.open(filename)
    scidata = band[0].data.astype(float)
    total_flux = 0

    for i in range(len(scidata)):
        for j in range(len(scidata)):
            if distance(i, j, 50, 50) >= 5 and distance(i, j, 50, 50) <= 42:
                total_flux += scidata[i][j]

    print(total_flux)
    total_flux /= (10**8 * 2000)
    total_flux = -2.5 * math.log10(total_flux)
    return total_flux


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
    print(stddev * math.sqrt(count))
    return stddev * math.sqrt(count)



"""
gflux = integrated_flux('MGSUBComb22-g.fit')#TotComb13-g.fit')# - 1.23
rflux = integrated_flux('MGSUBComb22-r.fit')#'TotComb13-r.fit')# - 0.52
iflux = integrated_flux('MGSUBComb22-i.fit')#'TotComb13-i.fit')# - 0.10
zflux = integrated_flux('MGSUBComb22-z.fit')#'TotComb13-z.fit')# - 0.07

gsed = sedfitter(gflux, 4770)
rsed = sedfitter(rflux, 6231)
ised = sedfitter(iflux, 7625)
zsed = sedfitter(zflux, 9134)
gerror_low = -2.5 * math.log10(math.sqrt(10**(gflux / -2.5) * 2000 * 10**8 / 4 + errorcalc('MGSUBComb22-g.fit')) / 2000 / 10**8 + 10**(gflux / -2.5))
rerror_low = -2.5 * math.log10(math.sqrt(10**(rflux / -2.5) * 2000 * 10**8 / 4.8 + errorcalc('MGSUBComb22-r.fit')) / 2000 / 10**8 + 10**(rflux / -2.5))
ierror_low = -2.5 * math.log10(math.sqrt(10**(iflux / -2.5) * 2000 * 10**8 / 4.8 + errorcalc('MGSUBComb22-i.fit')) / 2000 / 10**8 + 10**(iflux / -2.5))
zerror_low = -2.5 * math.log10(math.sqrt(10**(zflux / -2.5) * 2000 * 10**8  / 4.8 + errorcalc('MGSUBComb22-z.fit')) / 2000 / 10**8 + 10**(zflux / -2.5))
gerror_high = -2.5 * math.log10(-math.sqrt(10**(gflux / -2.5) * 2000 * 10**8  / 4 + errorcalc('MGSUBComb22-g.fit')) / 2000 / 10**8 + 10**(gflux / -2.5))
rerror_high = -2.5 * math.log10(-math.sqrt(10**(rflux / -2.5) * 2000 * 10**8  / 4.8 + errorcalc('MGSUBComb22-r.fit')) / 2000 / 10**8 + 10**(rflux / -2.5))
ierror_high = -2.5 * math.log10(-math.sqrt(10**(iflux / -2.5) * 2000 * 10**8  / 4.8 + errorcalc('MGSUBComb22-i.fit')) / 2000 / 10**8 + 10**(iflux / -2.5))
zerror_high = -2.5 * math.log10(-math.sqrt(10**(zflux / -2.5) * 2000 * 10**8  / 4.8 + errorcalc('MGSUBComb22-z.fit')) / 2000 / 10**8 + 10**(zflux / -2.5))

print("%s: %f, %f, %f, %s" % ('g-band', gflux, gerror_low, gerror_high, '%.4E' % Decimal(gsed)))
print("%s: %f, %f, %f, %s" % ('r-band', rflux, rerror_low, rerror_high, '%.4E' % Decimal(rsed)))
print("%s: %f, %f, %f, %s" % ('i-band', iflux, ierror_low, ierror_high, '%.4E' % Decimal(ised)))
print("%s: %f, %f, %f, %s" % ('z-band', zflux, zerror_low, zerror_high, '%.4E' % Decimal(zsed)))



fig, ax = plt.subplots()
error = [[sedfitter(gerror_low, 4770) - gsed, sedfitter(rerror_low, 6231) - rsed, sedfitter(ierror_low, 7625) - ised, sedfitter(zerror_low, 9134) - zsed], [gsed - sedfitter(gerror_high, 4770), rsed - sedfitter(rerror_high, 6231), ised - sedfitter(ierror_high, 7625), zsed - sedfitter(zerror_high, 9134)]]

#error *= np.array([[10, 10, 10, 10], [10, 10, 10, 10]])
#yerr = np.multiply(error, [[6, 6, 6, 6], [6, 6, 6, 6]])
yerr = error
xerr = [[770, 731, 775, 634], [730, 619, 875, 616]]
print(error)

(_, caps, _) = ax.errorbar([4770, 6231, 7625, 9134], [gsed, rsed, ised, zsed], xerr = xerr, yerr = yerr, capsize=5, elinewidth=2, fmt='bo')

for cap in caps:
    cap.set_color('black')
    cap.set_markeredgewidth(1.5)

plt.plot([4770, 6231, 7625, 9134], [4.9 * 10**-18, 9.0 * 10**-18, 9.35 * 10**-18, 1.09 * 10**-17], 'ro')


"""
def errorcalc(filename):
    band = fits.open(filename)
    scidata = band[0].data.astype(float)
    outside = []
    count = 0
    total_flux = 0
    
    for i in range(len(scidata)):
        for j in range(len(scidata)):
            if distance(i, j, 50, 50) >= 10 and distance(i, j, 50, 50) <= 60:
                outside.append(scidata[i][j])
                count += 1

    mean, median, stddev = sigma_clipped_stats(outside, sigma=3.0, iters=5)
    return stddev * math.sqrt(count)

fig, ax = plt.subplots()

total_flux = []

hdulist = fits.open('2175SUBComb30-g.fit')
scidata = hdulist[0].data.astype(float)
nmgy_count = 0

for i in range(len(scidata)):
    for j in range(len(scidata[0])):
        if distance(i, j, 50, 50) >= 3 and distance(i, j, 50, 50) <= 30:
            nmgy_count += scidata[i][j]

print(nmgy_count)
if nmgy_count < 0:
    nmgy_count = 0.000001
gmag = -2.5 * math.log10(nmgy_count) + 22.5
print(gmag)
total_flux.append(10**(gmag/-2.5) * 3631 * 10**-23 * 29979245800/(4770 * 10**-8)**2 * 10**-8)


hdulist = fits.open('2175SUBComb30-r.fit')
scidata = hdulist[0].data.astype(float)
nmgy_count = 0


for i in range(len(scidata)):
    for j in range(len(scidata[0])):
        if distance(i, j, 50, 50) >= 3 and distance(i, j, 50, 50) <= 30:
            nmgy_count += scidata[i][j]

print(nmgy_count)
if nmgy_count < 0:
    nmgy_count = 0.000001
rmag = -2.5 * math.log10(nmgy_count) + 22.5
print(rmag)
total_flux.append(10**(rmag/-2.5) * 3631 * 10**-23 * 29979245800/(6231 * 10**-8)**2 * 10**-8)


hdulist = fits.open('2175SUBComb30-i.fit')
scidata = hdulist[0].data.astype(float)
nmgy_count = 0

for i in range(len(scidata)):
    for j in range(len(scidata[0])):
        if distance(i, j, 50, 50) >= 3 and distance(i, j, 50, 50) <= 30:
            nmgy_count += scidata[i][j]

print(nmgy_count)
if nmgy_count < 0:
    nmgy_count = 0.000001
imag = -2.5 * math.log10(nmgy_count) + 22.5
print(imag)
total_flux.append(10**(imag/-2.5) * 3631 * 10**-23 * 29979245800/(7625 * 10**-8)**2 * 10**-8)


hdulist = fits.open('2175SUBComb30-z.fit')
scidata = hdulist[0].data.astype(float)
nmgy_count = 0

for i in range(len(scidata)):
    for j in range(len(scidata[0])):
        if distance(i, j, 50, 50) >= 3 and distance(i, j, 50, 50) <= 30:
            nmgy_count += scidata[i][j]

print(nmgy_count)
if nmgy_count < 0:
    nmgy_count = 0.000001
zmag = -2.5 * math.log10(nmgy_count) + 22.5
print(zmag)
total_flux.append(10**(zmag/-2.5) * 3631 * 10**-23 * 29979245800/(9134 * 10**-8)**2 * 10**-8)

gflux = gmag
rflux = rmag
iflux = imag
zflux = zmag

gsed = total_flux[0]
rsed = total_flux[1]
ised = total_flux[2]
zsed = total_flux[3]

plt.plot([4770, 6231, 7625, 9134], [gsed, rsed, ised, zsed], 'bo')

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



#fig, ax = plt.subplots()
error = [[sedfitter(gerror_low, 4770) - gsed, sedfitter(rerror_low, 6231) - rsed, sedfitter(ierror_low, 7625) - ised, sedfitter(zerror_low, 9134) - zsed], [gsed - sedfitter(gerror_high, 4770), rsed - sedfitter(rerror_high, 6231), ised - sedfitter(ierror_high, 7625), zsed - sedfitter(zerror_high, 9134)]]
#yerr = np.multiply(error, [[6, 6, 6, 6], [6, 6, 6, 6]])
yerr = error
xerr = [[770, 731, 775, 634], [730, 619, 875, 616]]

print(error)

ax.errorbar([4770, 6231, 7625, 9134], [gsed, rsed, ised, zsed], xerr = xerr, fmt='go')

(_, caps, _) = ax.errorbar([4770, 6231, 7625, 9134], [gsed, rsed, ised, zsed], xerr = xerr, capsize=5, elinewidth=2, fmt='go')

for cap in caps:
    cap.set_color('black')
    cap.set_markeredgewidth(1.5)

"""

plt.axis([3000, 10000, 0, 2 * 10**-17])
#majorLocator = MultipleLocator(10**-17)
#majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(0.25 * 10**-17)
#ax.yaxis.set_major_locator(majorLocator)
#ax.yaxis.set_major_formatter(majorFormatter)
ax.yaxis.set_minor_locator(minorLocator)

plt.xlabel('Wavelength [$\AA$]', fontsize=12)
plt.ylabel('Flux Density [erg sec^-1 cm^-2 $\AA$^-1]', fontsize=12)
MG_patch = mpatches.Patch(color='red', label='Zibetti et. al. SED')
#Zibetti_patch = mpatches.Patch(color='blue', label='Calculated SED')
#DUST_patch = mpatches.Patch(color='green', label='Dust Absorber SED')
plt.legend(handles=[MG_patch])
ax.set_title('SED of 2DA Host Galaxies, 1.00 <= Zabs < 2.50', y=1.05, fontsize=20)

plt.show()

