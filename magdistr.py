import string
import astropy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import math

from astropy.table import Table

reader = open('Full Data.txt', 'r')

gmag = []
rmag = []
imag = []

for i in range(1, 100001):
    line = reader.readline().split()
    gmag.append(float(line[6]))
    rmag.append(float(line[8]))
    imag.append(float(line[10]))


tableDUST = Table.read('final_catalog_full.fit', hdu=1)
drmag = []

for i in range(len(tableDUST)):
    drmag.append(tableDUST['R_BAND'][i])


"""
plt.plot(range(1, 100001), gmag, 'go')
plt.show()
plt.pause(4)
plt.close()

plt.plot(range(1, 100001), rmag, 'ro')
plt.show()
plt.pause(4)
plt.close()

plt.plot(range(1, 100001), imag, 'bo')
plt.show()
plt.pause(4)
plt.close()
"""

#print(gmag)

fig, ax = plt.subplots()


bins = np.arange(16, 23, 0.25)
arr = plt.hist(drmag, bins = bins, alpha = 0.25, color = 'r')
for i in range(len(bins) - 1):
    if arr[0][i] > 0:
        plt.text(arr[1][i],arr[0][i],str(int(arr[0][i])))

minorLocator = MultipleLocator(500)
ax.yaxis.set_minor_locator(minorLocator)
minorLocator = MultipleLocator(0.25)
ax.xaxis.set_minor_locator(minorLocator)
plt.xlabel('Apparent Magnitude')
plt.ylabel('Number of Absorber QSOs')
plt.title('Mg II Absorption QSO Apparent Magnitude Distribution in r-band')
#MG_patch = mpatches.Patch(color='red', label='')
#plt.legend(handles=[MG_patch])
plt.show()


