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
from scipy.constants import c, pi


def lambda_flambda_to_fnu(wavelength, flambda):
    """
    Convert a Fλ vs λ spectrum to Fν vs λ

    Parameters
    ----------
    wavelength: list-like of floats
        The wavelengths in nm.
    flambda: list-like of floats
        Fλ flux density in W/m²/nm (or Lλ luminosity density in W/nm).

    Returns
    -------
    fnu: array of floats
        The Fν flux density in mJy (or the Lν luminosity density in
        1.e-29 W/Hz).

    """
    wavelength = np.array(wavelength, dtype=float)
    flambda = np.array(flambda, dtype=float)

    # Factor 1e+29 is to switch from W/m²/Hz to mJy
    # Factor 1e-9 is to switch from nm to m (only one because the other nm
    # wavelength goes with the Fλ in W/m²/nm).
    fnu = 1e+29 * 1e-9 * flambda * wavelength * wavelength / c

    return fnu





#fluxtable = Table.open('sb07_burst_0.1_z_0.008_chab.model', hdu=1)
table = Table.read('MGIIAbsorbers_best_model.fits', hdu=1)
#fluxdata = file[0].data

#for i in range(len(fluxdata)):
#	fluxdata[i] *= 


fig, ax = plt.subplots()
plt.plot(10 * table['wavelength'][400 : 1120], table['Fnu'][400 : 1120], 'black')

#nebular_continuum_old = lambda_flambda_to_fnu(table['wavelength'], table['nebular.lines_old'])
#plt.plot(10 * table['wavelength'][237 : 1115], nebular_continuum_old[237 : 1115])
#plt.plot([486, 628, 770, 913], [0.01028078, 0.01721973, 0.0178659554, 0.0328115124], 'ro')

xerr = [[770, 731, 775, 634], [730, 619, 875, 616]]
error2 = [[0.0008196941, 0.001300651, 0.001591861, 0.004950240], [0.0008196941, 0.001300651, 0.001591861, 0.004950240]]



# K-corrected: [0.01028078, 0.01721973, 0.0178659554, 0.0328115124], [0.00141492, 0.00166133, 0.00172368, 0.0055615314]

(_, caps, _) = ax.errorbar([4770, 6231, 7625, 9134], [0.004112112, 0.010353394, 0.015508566, 0.027546101], xerr = xerr, yerr = error2, capsize=5, elinewidth=2, fmt='bo')

for cap in caps:
    cap.set_color('blue')
    cap.set_markeredgewidth(1.5)




#plt.axes().yaxis.set_tick_params(which='minor', right = 'off')
plt.xscale('log')
plt.xlabel('Wavelength [$\AA$]', fontsize = 12)
plt.ylabel('Flux Density [mJy]', fontsize = 12)
plt.title('Mg II Absorber Host Galaxy SED 0.37 $\leq$ $z_{abs}$ < 0.55', fontsize=20, y=1.04)

MG_patch = mpatches.Patch(color='blue', label='SDSS Photometric Fluxes')
NET_patch = mpatches.Patch(color='black', label='CIGALE SED Fitting')
plt.legend(handles=[MG_patch, NET_patch])


plt.show()