from astropy.io import fits
from shutil import copyfile
import linecache
import os

gdirs = os.listdir('/data/marvels/billzhu/Reference PSF Subtract/0.37 - 0.55/g/')
objdirs = os.listdir('/data/marvels/billzhu/Reference Obj/0.37 - 0.55/')
print(objdirs)

for i in range(len(gdirs)):
    index = int(gdirs[i].split('-')[0])
    print(index)
    if str(index) + '.fit' not in objdirs:
        copyfile('/data/marvels/jerryxu/dr7/raw_catalog/' + str(index) + '.fit', '/data/marvels/billzhu/Reference Obj/0.37 - 0.55/' + str(index) + '.fit')
