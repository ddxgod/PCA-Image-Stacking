import ezgal
import numpy as np

model = ezgal.model('www.baryons.org/ezgal/models/bc03_burst_0.1_z_0.008_salp.model')
#model = ezgal.model('www.baryons.org/ezgal/models/bc03_ssp_z_0.02_chab.model')

# Desired formation redshift
zf = 1.4
# Fetch an array of redshifts out to given formation redshift
zs = np.arange(0.35, zf, 0.005)
#print(np.shape(zs))
model.set_normalization('sloan_g', 0.00000000000000001, -22.4)
thresh_mags = model.get_apparent_mags(zf, filters=['sloan_u', 'sloan_g', 'sloan_r', 'sloan_i', 'sloan_z'], zs=zs)

thefile = open('Ezgal mag.txt', 'w')
for item in thresh_mags:
  thefile.write("%s\n" % item)


print(thresh_mags)