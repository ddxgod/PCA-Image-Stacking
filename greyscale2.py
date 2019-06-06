from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter




def sbgrey(data, i, j, fig, ax, majorLocator, majorFormatter, minorLocator):
	sb = np.zeros((101, 101))

	for r in range(len(data[0])):
		for c in range(len(data[0])):
			if data[r, c] > 0:
				#sb[r, c] = 22.5 - 2.5 * np.log10(data[r, c]) + 2.5 * np.log10(0.396**2)
				sb[r, c] = -2.5 * np.log10(data[r, c] / 2000 / 10**8) + 2.5 * np.log10(0.396**2)
			else:
				sb[r, c] = 32



	ax[i, j].yaxis.set_major_locator(majorLocator)
	ax[i, j].yaxis.set_major_formatter(majorFormatter)
	ax[i, j].yaxis.set_minor_locator(minorLocator)
	ax[i, j].xaxis.set_major_locator(majorLocator)
	ax[i, j].xaxis.set_major_formatter(majorFormatter)
	ax[i, j].xaxis.set_minor_locator(minorLocator)
	ax[i, j].axis([0, 100, 0, 100])

	return sb, ax








fig, ax = plt.subplots(2, 2)
fig.set_size_inches(10, 10)
#plt.pcolor(cmap='grey', norm=LogNorm(vmin=32, vmax=vmax))
majorLocator = MultipleLocator(50)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(10)

vmin = 26
vmax = 32
gamma = 0.17


"""
hdulist = fits.open('TestPSF400-g.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 0, 0, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[0].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[0].set_title('$g$-band 400 Mg II Stack', y=1.0, fontsize=18)


hdulist = fits.open('2175SUBComb30-g.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 0, 1, fig, ax, majorLocator, majorFormatter, minorLocator)
im=ax[1].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[1].set_title('$g$-band 400 2DA Stack', y=1.0, fontsize=18)
"""




hdulist = fits.open('StarComb-g.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 0, 0, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[0, 0].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[0, 0].set_title('g-band PSF Star Stack', y=1.0, fontsize=14)

hdulist = fits.open('StarComb-r.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 0, 1, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[0, 1].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[0, 1].set_title('r-band PSF Star Stack', y=1.0, fontsize=14)

hdulist = fits.open('StarComb-i.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 1, 0, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[1, 0].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[1, 0].set_title('i-band PSF Star Stack', y=1.0, fontsize=14)

hdulist = fits.open('StarComb-z.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 1, 1, fig, ax, majorLocator, majorFormatter, minorLocator)
im = ax[1, 1].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[1, 1].set_title('z-band PSF Star Stack', y=1.0, fontsize=14)







"""
hdulist = fits.open('2175SUBComb30-g.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 0, 0, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[0, 0].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[0, 0].set_title('g-band 2DA Stack', y=1.0, fontsize=14)

hdulist = fits.open('2175SUBComb30-r.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 0, 1, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[0, 1].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[0, 1].set_title('r-band 2DA Stack', y=1.0, fontsize=14)

hdulist = fits.open('2175SUBComb30-i.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 1, 0, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[1, 0].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[1, 0].set_title('i-band 2DA Stack', y=1.0, fontsize=14)

hdulist = fits.open('2175SUBComb30-z.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 1, 1, fig, ax, majorLocator, majorFormatter, minorLocator)
im = ax[1, 1].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[1, 1].set_title('z-band 2DA Stack', y=1.0, fontsize=14)
"""




"""
# g-band
hdulist = fits.open('MGSUBComb50-g.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 0, 0, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[0, 0].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[0, 0].set_title('g-band Mg II Stack', y=1.0, fontsize=10)

hdulist = fits.open('REFSUBComb50-g.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 0, 1, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[0, 1].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[0, 1].set_title('g-band Reference QSO Stack', y=1.0, fontsize=10)

hdulist = fits.open('TotComb50-g.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 0, 2, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[0, 2].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[0, 2].set_title('g-band Net Stack', y=1.0, fontsize=10)




# r-band
hdulist = fits.open('MGSUBComb50-r.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 1, 0, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[1, 0].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[1, 0].set_title('r-band Mg II Stack', y=1.0, fontsize=10)

hdulist = fits.open('REFSUBComb50-r.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 1, 1, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[1, 1].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[1, 1].set_title('r-band Reference QSO Stack', y=1.0, fontsize=10)

hdulist = fits.open('TotComb50-r.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 1, 2, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[1, 2].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[1, 2].set_title('r-band Net Stack', y=1.0, fontsize=10)




# i-band
hdulist = fits.open('MGSUBComb50-i.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 2, 0, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[2, 0].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[2, 0].set_title('i-band Mg II Stack', y=1.0, fontsize=10)

hdulist = fits.open('REFSUBComb50-i.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 2, 1, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[2, 1].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[2, 1].set_title('i-band Reference QSO Stack', y=1.0, fontsize=10)

hdulist = fits.open('TotComb50-i.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 2, 2, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[2, 2].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[2, 2].set_title('i-band Net Stack', y=1.0, fontsize=10)




# z-band
hdulist = fits.open('MGSUBComb50-z.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 3, 0, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[3, 0].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[3, 0].set_title('z-band Mg II Stack', y=1.0, fontsize=10)

hdulist = fits.open('REFSUBComb50-z.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 3, 1, fig, ax, majorLocator, majorFormatter, minorLocator)
ax[3, 1].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[3, 1].set_title('z-band Reference QSO Stack', y=1.0, fontsize=10)

hdulist = fits.open('TotComb50-z.fit')
data = hdulist[0].data.astype(float)
sb, ax = sbgrey(data, 3, 2, fig, ax, majorLocator, majorFormatter, minorLocator)
im = ax[3, 2].imshow(sb, cmap='gray', norm=PowerNorm(vmin=vmin, vmax=vmax, gamma=gamma))
ax[3, 2].set_title('z-band Net Stack', y=1.0, fontsize=10)



#plt.xlabel('Pixel = 0.396 arcsec ~ 3.4 kpc', fontsize = 10)
#plt.ylabel('Pixel = 0.396 arcsec ~ 3.4 kpc', fontsize = 10)
#plt.title(str(round(0.4 * i, 1)) + ' arcsec impact parameter', y=1.0, fontsize=20)


fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.grid(False)
plt.xlabel('Pixel = 0.396 arcsec ~ 2.4 kpc', y=1.0, fontsize = 12)
plt.ylabel('Pixel = 0.396 arcsec ~ 2.4 kpc', x=1.2, fontsize = 12)
#plt.title('$g$-band 400 Stacked Images of Mg II and 2DA Systems', y=1.04, fontsize=24)
#plt.title('PSF Star Stacked Frames in $griz$ Bands', y=1.07, fontsize=24)
"""

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])


cbar = fig.colorbar(im, cax=cbar_ax, ticks=[26, 26.5, 27, 28, 30])
cbar.set_label('mag arcsec$^{-2}$', rotation=270, fontsize=16, labelpad=25)



"""
#plt.tight_layout()
plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
plt.setp([a.get_xticklabels() for a in ax[1, :]], visible=False)
plt.setp([a.get_xticklabels() for a in ax[2, :]], visible=False)
plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
plt.setp([a.get_yticklabels() for a in ax[:, 2]], visible=False)
"""

plt.show()