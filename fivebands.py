from astropy.io import fits
import os
from shutil import copyfile

gdirs = os.listdir('/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/g/')

for i in range(len(gdirs)):
    index = int(gdirs[i].split('-')[0])
    mg = int(gdirs[i].split('-')[1])
    print(index)
    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(index) + '-r.fit', '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/r/' + str(index) + '-' + str(mg) + '-r.fit')
    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(index) + '-i.fit', '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/i/' + str(index) + '-' + str(mg) + '-i.fit')
    copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(index) + '-z.fit', '/data/marvels/billzhu/Reference Dataset/0.37 - 0.55/z/' + str(index) + '-' + str(mg) + '-z.fit')
    #copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(index) + '-u.fit', '/data/marvels/billzhu/MG II Dataset/0.37 - 0.55/u/' + str(index) + '-u.fit')
