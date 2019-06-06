import numpy as np
import astropy
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
import multiprocessing
from multiprocessing import Pool
import os
import matplotlib.pyplot as plt
import pickle
import math


# Reference Background g-band: (5.860117404482102, 6.0, 1.6851090617072355)
# Reference Background r-band: (5.893401986391389, 6.0, 1.804174203392319)
# Reference Background i-band: (6.042582374762306, 6.0, 2.003988762834255)
# Reference Background z-band: (4.6067475093086445, 4.0, 1.6801096788726562)




def begin(index, tot_mean, tot_median, tot_stddev, tot_400):
	i = int(index.split('-')[0])
	mgi = int(index.split('-')[1])
	color = index.split('-')[2].split('_')[0]
	number = index.split('-')[2].split('_')[1]

	bkg_table = Table.read('/data/marvels/billzhu/Reference Background/' + color + '/' + index + '.fits', hdu=1)
	median_arr = bkg_table['median']

	for j in range(len(bkg_table['median']) - 1):
		if bkg_table['median'][j] == bkg_table['median'][j + 1]:
			bkg_limit.append(j)
			print(counter)
			break


	tot_mean += bkg_table['mean']
	tot_median += bkg_table['median']
	tot_stddev += bkg_table['stddev']
	columns = bkg_table.columns
	#print(columns[3])
	tot_400 += columns[3][:3]





def begin2(index, res_arr):
	i = int(index.split('-')[0])
	mgi = int(index.split('-')[1])
	color = index.split('-')[2].split('_')[0]
	number = index.split('-')[2].split('_')[1]

	bkg_table = Table.read('/data/marvels/billzhu/Reference Background/' + color + '/' + index + '.fits', hdu=1)
	column_data = bkg_table.columns
	names = bkg_table.colnames
	lower_lim = math.ceil(float(names[3].split('-')[0]))
	higher_lim = math.ceil(float(names[3].split('-')[1].split()[0]))
	res_arr.append(np.mean(column_data[0][lower_lim - 1: higher_lim]) - column_data[3][0])
	print(counter)





if __name__ == '__main__':
	bkg_limit = []

	gdirs = os.listdir('/data/marvels/billzhu/Reference Background/g/')
	rdirs = os.listdir('/data/marvels/billzhu/Reference Background/r/')
	idirs = os.listdir('/data/marvels/billzhu/Reference Background/i/')
	zdirs = os.listdir('/data/marvels/billzhu/Reference Background/z/')

	rangelist = []


	counter = 0
	tot_mean = np.zeros(200)
	tot_median = np.zeros(200)
	tot_stddev = np.zeros(200)
	tot_400 = np.zeros(3)

	sigma_mean = 0

	res_arr = []




	"""
	for d in gdirs:
		index = d.split('.')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		#begin(index, tot_mean, tot_median, tot_stddev, tot_400)
		begin2(index, res_arr)
		counter += 1
	
	
	
	for d in rdirs:
		index = d.split('.')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		#begin(index, tot_mean, tot_median, tot_stddev, tot_400)
		begin2(index, res_arr)
		counter += 1
	
	
	
	for d in idirs:
		index = d.split('.')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		begin(index, tot_mean, tot_median, tot_stddev, tot_400)
		begin2(index, res_arr)
		counter += 1
	
	"""
	
	for d in zdirs:
		index = d.split('.')[0]
		#if str(index) + '_SUB.fit' not in check_dirs:
		#begin(index, tot_mean, tot_median, tot_stddev, tot_400)
		begin2(index, res_arr)
		counter += 1
	
	"""

	
	print(len(rangelist))
	#try:
	pool = Pool(multiprocessing.cpu_count())
	pool.map(begin, rangelist)
	

	tot_mean /= counter
	tot_median /= counter
	tot_stddev /= counter
	tot_400 /= counter

	print(bkg_limit)
	print(len(bkg_limit))
	print(sigma_clipped_stats(bkg_limit, sigma=3.0, iters=5))

	with open('z-bkg.txt', 'wb') as fp:
		pickle.dump(tot_mean, fp)
		#pickle.dump(tot_median, fp)
		#pickle.dump(tot_stddev, fp)


	with open('z2-bkg.txt', 'wb') as fp:
		pickle.dump(tot_400, fp)


	with open('z-stddev.txt', 'wb') as fp:
		pickle.dump(tot_stddev, fp)
	
	"""


	with open('z3-bkg.txt', 'wb') as fp:
		pickle.dump(res_arr, fp)



