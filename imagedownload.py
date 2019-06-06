import string
import urllib
from urllib.request import urlretrieve
import linecache

reader = open('Full Data.txt', 'r')
#for i in range(1, 11):
data = linecache.getline('Full Data.txt', 9816).split()

if len(data[51]) < 2:
    data[51] = "00" + data[51]
if len(data[51]) < 3:
    data[51] = "0" + data[51]

link = format('http://das.sdss.org/imaging/%s/%s/corr/%s/' % (data[44], data[49], data[50]))

if len(data[44]) < 3:
    data[44] = '00' + data[44]
if len(data[44]) < 4:
    data[44] = '0' + data[44]

link = format("%sfpC-00%s-g%s-0%s.fit.gz" % (link, data[44], data[50], data[51]))
print(link)

try:
    urlretrieve(link, 'c:/Research Project/' + str(9816) + ".fit")
except:
    print("Unable to download")
