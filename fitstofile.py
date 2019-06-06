from astropy.table import Table
import numpy as np

writer = open('DR12 QSO.txt', 'w')
dust_table = Table.read('DR12Q.fits', hdu=1)

for i in range(len(dust_table)):
    print(i)
    writer.write(' '.join(str(s) for s in dust_table[i]) + '\n')

writer.close()
