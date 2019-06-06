import numpy as np
import scipy as sp
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u



scale = cosmo.kpc_proper_per_arcmin(0.37) * u.arcmin / u.kiloparsec * 0.396 / 60
b400 = 400 / scale * 0.396
b500 = 500 / scale * 0.396
b10 = 10 / scale * 0.396
b100 = 100 / scale * 0.396
scale = cosmo.kpc_proper_per_arcmin(0.55) * u.arcmin / u.kiloparsec * 0.396 / 60
t400 = 400 / scale * 0.396
t500 = 500 / scale * 0.396



fig, ax = plt.subplots()

majorLocator = MultipleLocator(10)
majorFormatter = FormatStrFormatter('%.1f')
minorLocator = MultipleLocator(2.5)
minorFormatter = FormatStrFormatter('%f')



with open('g-bkg.txt', 'rb') as fp:
	arr = pickle.load(fp)


with open('g2-bkg.txt', 'rb') as fp2:
	arr2 = pickle.load(fp2)


with open('g-stddev.txt', 'rb') as fp3:
	arr3 = pickle.load(fp3)


print(arr)
arr = arr[2: 101]
arr3 = arr3[2 : 101]
ax.xaxis.set_major_locator(majorLocator)
ax.xaxis.set_major_formatter(majorFormatter)
ax.set_xlabel('Arcsec From Center', fontsize=12)
ax.set_ylabel('Mean Background Count', fontsize=12)
ax.set_title('g-band Mean Background Count at Increasing Arcsec Annuli', fontsize=16)
ax.plot(range(3, 102), arr)
ax.plot(range(3, 102), arr - arr3, color='black', linestyle='--')
ax.plot(range(3, 102), arr + arr3, color='black', linestyle='--')

ax.axhline(y=arr2[0], color='red', linestyle='--')
print(arr2[0])
print(arr2[2])
ax.axvspan(t400, b400, alpha=0.2, color='green')
ax.axvspan(t500, b500, alpha=0.2, color='blue')
ax.axvspan(b10, b100, alpha=0.1, color='red')



"""
with open('r-bkg.txt', 'rb') as fp:
	arr = pickle.load(fp)


with open('r2-bkg.txt', 'rb') as fp2:
	arr2 = pickle.load(fp2)


with open('r-stddev.txt', 'rb') as fp3:
	arr3 = pickle.load(fp3)


arr = arr[4: 101]
arr3 = arr3[4 : 101]
ax[1, 0].xaxis.set_major_locator(majorLocator)
ax[1, 0].xaxis.set_major_formatter(majorFormatter)
ax[1, 0].set_xlabel('Arcsec From Center', fontsize=12)
ax[1, 0].set_ylabel('Mean Background Count', fontsize=12)
ax[1, 0].set_title('r-band Mean Background Count at Increasing Arcsec Annuli', fontsize=16)
ax[1, 0].plot(range(5, 102), arr)
ax[1, 0].plot(range(5, 102), arr - arr3, color='black')
ax[1, 0].plot(range(5, 102), arr + arr3, color='black')
print(arr2[0])
print(arr2[2])

ax[1, 0].axhline(y=arr2[0], color='red', linestyle='--')
ax[1, 0].axvspan(t400, b400, alpha=0.2, color='green')
ax[1, 0].axvspan(t500, b500, alpha=0.2, color='blue')
ax[1, 0].axvspan(b10, b100, alpha=0.2, color='red')




with open('i-bkg.txt', 'rb') as fp:
	arr = pickle.load(fp)


with open('i2-bkg.txt', 'rb') as fp2:
	arr2 = pickle.load(fp2)


with open('i-stddev.txt', 'rb') as fp3:
	arr3 = pickle.load(fp3)


arr = arr[4: 101]
arr3 = arr3[4 : 101]
ax[0, 1].xaxis.set_major_locator(majorLocator)
ax[0, 1].xaxis.set_major_formatter(majorFormatter)
ax[0, 1].set_xlabel('Arcsec From Center', fontsize=12)
ax[0, 1].set_ylabel('Mean Background Count', fontsize=12)
ax[0, 1].set_title('i-band Mean Background Count at Increasing Arcsec Annuli', fontsize=16)
ax[0, 1].plot(range(5, 102), arr)
ax[0, 1].plot(range(5, 102), arr - arr3, color='black')
ax[0, 1].plot(range(5, 102), arr + arr3, color='black')
print(arr2[0])
print(arr2[2])

ax[0, 1].axhline(y=arr2[0], color='red', linestyle='--')
ax[0, 1].axvspan(t400, b400, alpha=0.2, color='green')
ax[0, 1].axvspan(t500, b500, alpha=0.2, color='blue')
ax[0, 1].axvspan(b10, b100, alpha=0.2, color='red')




with open('z-bkg.txt', 'rb') as fp:
	arr = pickle.load(fp)


with open('z2-bkg.txt', 'rb') as fp2:
	arr2 = pickle.load(fp2)


with open('z-stddev.txt', 'rb') as fp3:
	arr3 = pickle.load(fp3)


arr = arr[4: 101]
arr3 = arr3[4: 101]
ax[1, 1].xaxis.set_major_locator(majorLocator)
ax[1, 1].xaxis.set_major_formatter(majorFormatter)
ax[1, 1].set_xlabel('Arcsec From Center', fontsize=12)
ax[1, 1].set_ylabel('Mean Background Count', fontsize=12)
ax[1, 1].set_title('z-band Mean Background Count at Increasing Arcsec Annuli', fontsize=16)
ax[1, 1].plot(range(5, 102), arr)
ax[1, 1].plot(range(5, 102), arr - arr3, color='black')
ax[1, 1].plot(range(5, 102), arr + arr3, color='black')
print(arr2[0])
print(arr2[2])

ax[1, 1].axhline(y=arr2[0], color='red', linestyle='--')
ax[1, 1].axvspan(t400, b400, alpha=0.2, color='green')
ax[1, 1].axvspan(t500, b500, alpha=0.2, color='blue')
ax[1, 1].axvspan(b10, b100, alpha=0.2, color='red')
"""

print("%f, %f, %f, %f" % (b400, t400, b500, t500))


fig.set_figheight(6)
fig.set_figwidth(8)
ax.set_xscale('log')
#ax[0, 1].set_xscale('log')
#ax[1, 0].set_xscale('log')
#ax[1, 1].set_xscale('log')
#plt.yscale('log')
#plt.tight_layout()
#plt.setp([a.get_xticklabels() for a in ax[0, :]], visible=False)
#plt.setp([a.get_xticklabels() for a in ax[1, :]], visible=False)
#plt.setp([a.get_yticklabels() for a in ax[:, 1]], visible=False)
plt.show()