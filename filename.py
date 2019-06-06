from astropy.io import fits
import linecache


writefile = open('Reference QSO Data.txt', 'w')
indexfile = open('Reference QSO.txt', 'r')

count = len(open('Reference QSO.txt').readlines(  )) - 1
writefile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("run", "imgMJD", "specMJD", "plate", "fiber", "rerun", "column", "frame", "object"))

total = indexfile.readline()

for i in range(count):
	print(i)
	data = indexfile.readline()
	ref = int(data.split()[0])
	mg = int(data.split()[1])
	line_data = linecache.getline('Full Data.txt', ref).split()
	writefile.write("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n" % (int(line_data[44]), int(line_data[45]), int(line_data[46]), int(line_data[47]), int(line_data[48]), int(line_data[49]), int(line_data[50]), int(line_data[51]), int(line_data[52])))

writefile.close()