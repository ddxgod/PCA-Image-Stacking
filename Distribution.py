from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

table_list = fits.open('Trimmed_SDSS_DR7_107.fits')
table_readMG = Table.read('QSObased_Trimmed_SDSS_DR7_107.fits')
table_readDUST = Table.read('final_catalog_full.fit')

redshift_MG = []
redshift_DUST = []
average = []

fit, ax = plt.subplots(1, 2)


for i in range(len(table_readMG)):
    redshift_MG.append(table_readMG['ZABS'][i][0])
    if table_readMG['ZABS'][i][0] >= 0.37 and table_readMG['ZABS'][i][0] < 0.55 and table_readMG['REW_MGII_2796'][i][0] > 0.8 and table_readMG['NABS'][i] == 1:
        average.append(table_readMG['ZABS'][i][0])

average = np.mean(average)
print(average)

average = []
for i in range(len(table_readDUST)):
    redshift_DUST.append(table_readDUST['ZABS'][i])
    average.append(table_readDUST['ZABS'][i])

average = np.mean(average)
print(average)


bins_MG = np.arange(0.2, 2.5, 0.1)
bins_DUST = np.arange(1.0, 2.5, 0.1)
arr_MG = ax[0].hist(redshift_MG, bins = bins_MG, alpha = 0.1, color = 'b')
arr_DUST = ax[1].hist(redshift_DUST, bins = bins_DUST, alpha = 0.1, color = 'g')


"""
for i in range(len(bins) - 1):    
    if arr_MG[0][i] > 0:
        plt.text(arr_MG[1][i], arr_MG[0][i], str(int(arr_MG[0][i])))

    #if arr_DUST[0][i] > 0:
    #    plt.text(arr_DUST[1][i], arr_DUST[0][i], str(int(arr_DUST[0][i])))
"""

ax[0].axvline(x=0.37, color='r', linestyle='--')
ax[0].axvline(x=0.55, color='r', linestyle='--')
#MG_patch = mpatches.Patch(color='blue', label='Mg II Absorber')
#DUST_patch = mpatches.Patch(color='green', label='2175 $\AA$ Dust Absorber')

ax[0].set_xticks(np.arange(0.2, 2.5, 0.1))
ax[0].set_yticks(np.arange(0, 3600, 400))
ax[0].set_xlabel('Redshift', fontsize=14)
ax[0].set_ylabel('Number of Absorber QSOs', fontsize=14)
ax[0].set_title('Redshift Distribution of Mg II Absorbers', fontsize=20)


ax[1].set_xticks(np.arange(1.0, 2.6, 0.1))
ax[1].set_yticks(np.arange(0, 50, 10))
ax[1].set_xlabel('Redshift', fontsize=14)
ax[1].set_ylabel('Number of Absorber QSOs', fontsize=14)
ax[1].set_title('Redshift Distribution of 2175 $\AA$ Dust Absorbers', fontsize=20)
plt.show()
