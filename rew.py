from astropy.table import Table
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import numpy as np

da = Table.read('final_catalog_full.fit')
print(np.mean(da['EW2796']))
print(sigma_clipped_stats(da['EW2796'], sigma=3.0, iters=5)[0])

bins = np.arange(0, 5, 0.1)

arr = plt.hist(np.array(da['EW2796']).T, bins = bins, alpha = 0.1, color = 'r')
print(np.std(da['EW2796']))
plt.show()