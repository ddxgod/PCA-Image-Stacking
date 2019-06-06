from shutil import copyfile
from astropy.io import fits
from astropy.table import Table
import os
import os.path
import numpy as np
import scipy as sp
import math
import string
import linecache


copyfile('/data/marvels/jerryxu/dr7/catalog/' + str(9816) + '-r.fit', '/data/marvels/billzhu/' + str(9816) + '-r.fit')
