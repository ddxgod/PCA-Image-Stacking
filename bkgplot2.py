import numpy as np
import scipy as sp
import pickle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
from astropy.table import Table


fig, ax = plt.subplots()
majorLocator = MultipleLocator(1000)
#majorFormatter = FormatStrFormatter('%.1f')
minorLocator = MultipleLocator(500)
minorFormatter = FormatStrFormatter('%f')


with open('z3-bkg.txt', 'rb') as fp:
	arr = pickle.load(fp)

ax.hist(arr, bins=100,  facecolor='blue', alpha=0.3)
ax.set_title('z-band Mean Residue', fontsize=20)
ax.yaxis.set_major_locator(majorLocator)
ax.yaxis.set_minor_locator(minorLocator)


"""
with open('r3-bkg.txt', 'rb') as fp:
	arr = pickle.load(fp)

ax[1].hist(arr, bins=100,  facecolor='blue', alpha=0.3)
ax[1].set_title('r-band Mean Residue', fontsize=16)
#ax[0, 1].yaxis.set_major_locator(majorLocator)
#ax[0, 1].yaxis.set_minor_locator(minorLocator)



with open('i3-bkg.txt', 'rb') as fp:
	arr = pickle.load(fp)

ax[2].hist(arr, bins=100,  facecolor='blue', alpha=0.3)
ax[2].set_title('i-band Mean Residue', fontsize=16)

#ax[1, 0].yaxis.set_major_locator(majorLocator)
#ax[1, 0].yaxis.set_minor_locator(minorLocator)



with open('z3-bkg.txt', 'rb') as fp:
	arr = pickle.load(fp)

ax[3].hist(arr, bins=100,  facecolor='blue', alpha=0.3)
ax[3].set_title('z-band Mean Residue', fontsize=16)
#ax[1, 1].yaxis.set_major_locator(majorLocator)
#ax[1, 1].yaxis.set_minor_locator(minorLocator)
"""





#fig.add_subplot(111, frameon=False)
#plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
#fig.(0.3, 0, 'Mean residue left over 400 - 500 kpc post backround subtraction', fontsize=16)
#(0, 0, 'Number of PSF stars', fontsize=16, rotation='vertical')
#for axe in ax.flat:
#ax.set(xlabel='Mean residue in counts post backround subtraction', ylabel='Number of PSF stars')
plt.xlabel('Mean residue in counts post backround subtraction', fontsize=12)
plt.ylabel('Number of PSF stars', fontsize=12)


#fig.tight_layout()

fig.subplots_adjust(hspace=0.8)

plt.show()
#fig.savefig('Mg II bkg residue.png')

