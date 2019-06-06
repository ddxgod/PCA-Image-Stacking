from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import scipy as sp
import os
import linecache
import math
import string
from scipy import interpolate
from astropy.cosmology import WMAP9 as cosmo
from astropy import units as u
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from decimal import Decimal


ngc2865 = open('NGC 2865 Spectrum Data.txt', 'r')
wavelength = []
data = []

ngc2865.readline()
ngc2865.readline()

for line in ngc2865:
	wavelength.append(float(line.split())[0])
	counts = float(line.split()[1])
	mag = -2.5 * math.log10(counts)

print(data)
