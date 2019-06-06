from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


#hdulist = fits.open('ExpTotComb0_6.fit')
hdulist = fits.open('SimulatedPSFSUBComb14-r5.fit')
data = hdulist[0].data.astype(float)
#data[0] -= hdulist[0].data[0] * 3/4
#data -= np.median(data)



fig, ax = plt.subplots(2, 4)
#plt.pcolor(cmap='grey', norm=LogNorm(vmin=32, vmax=100))
majorLocator = MultipleLocator(50)
majorFormatter = FormatStrFormatter('%d')
minorLocator = MultipleLocator(10)

for i in range(len(data)):
	sb = np.zeros(np.shape(data[i]))

	for r in range(len(data[0])):
		for c in range(len(data[0])):
			if data[i][r, c] > 0:
				sb[r, c] = 22.5 - 2.5 * np.log10(data[i][r, c]) + 2.5 * np.log10(0.396**2)
				#sb[r, c] = -2.5 * np.log10(data[r, c] / 2000 / 10**8) + 2.5 * np.log10(2.4**2)
			else:
				sb[r, c] = 34



	ax[i//4, i % 4].yaxis.set_major_locator(majorLocator)
	ax[i//4, i % 4].yaxis.set_major_formatter(majorFormatter)
	ax[i//4, i % 4].yaxis.set_minor_locator(minorLocator)
	ax[i//4, i % 4].xaxis.set_major_locator(majorLocator)
	ax[i//4, i % 4].xaxis.set_major_formatter(majorFormatter)
	ax[i//4, i % 4].xaxis.set_minor_locator(minorLocator)
	
	im = ax[i//4, i % 4].imshow(sb, cmap='gray', norm=PowerNorm(vmin=25, vmax=34, gamma=0.17))
	ax[i//4, i % 4].axis([0, 100, 0, 100])

	#ax[i//4, i % 4].set_xlabel('Pixel = 0.396 arcsec ~ 3.4 kpc', fontsize = 8)
	#ax[i//4, i % 4].set_ylabel('Pixel = 0.396 arcsec ~ 3.4 kpc', fontsize = 8)
	ax[i//4, i % 4].set_title(str(round(0.4 * i, 1)) + ' arcsec', y=1.03, fontsize=12)

	#plt.imshow(sb, cmap='gray', norm=PowerNorm(vmin=32, vmax=100, gamma=0.2))

#plt.xlabel('Pixel = 0.396 arcsec ~ 3.4 kpc', fontsize = 12)
#plt.ylabel('Pixel = 0.396 arcsec ~ 3.4 kpc', fontsize = 12)
#plt.title(str(round(0.4 * i, 1)) + ' arcsec impact parameter', y=1.03, fontsize=20)


fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
plt.grid(False)
plt.xlabel('Pixel = 0.396 arcsec ~ 3.4 kpc', y=1.03, fontsize = 18)
plt.ylabel('Pixel = 0.396 arcsec ~ 3.4 kpc', x=1.08, fontsize = 18)
#plt.title('Seeing FWHM = 0.6 arcsec, residue image at different impact parameters', y=1.07, fontsize=24)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])

cbar = fig.colorbar(im, cax=cbar_ax, ticks=[25, 26, 27.5, 30.0])
cbar.set_label('mag arcsec$^{-2}$', rotation=270, fontsize=16, labelpad=25)

plt.show()
plt.tight_layout()
plt.colorbar()
plt.show()