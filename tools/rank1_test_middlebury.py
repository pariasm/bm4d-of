# test the low-rank covariance matrix for the 2nd step

import subprocess
import shutil
import time
from math import floor
import os
import glob

sequences = [
'Army',
'Walking',
'DogDance',
'Mequon',
'Evergreen']
#'Army',
#'Backyard',
#'Basketball']#,
#'Beanbags',
#'DogDance',
#'Dumptruck',
#'Evergreen',
#'Grove',
#'Grove2',
#'Grove3',
#'Hydrangea',
#'Mequon',
#'MiniCooper',
#'RubberWhale',
#'Schefflera',
#'Urban',
#'Urban2',
#'Urban3',
#'Walking',
#'Wooden',
#'Yosemite'

base_dir = '/home/pariasm/ipol/local/video_nlbayes/input/'

sigma_str  = ['40','20','10']
rank1      = ['4','8','12','16','20']
npatches1  = ['80','120','160','300','600']
method = 'vnlbayes'

for np1 in npatches1:
	for r1 in rank1:
		for sigma in sigma_str:
			for iseq in range(0,len(sequences)):
				seq = sequences[iseq]
				print "Processing " + seq + " with r1 = " + r1 + " and np1 " + np1

				# time execution
				start_time = time.time()

				# launch process
				proc = subprocess.Popen(['bin/' + method,
						'-i', base_dir + seq + '/i%04d.png',
						'-f', '0', '-l', '7',
						'-sigma', sigma,
						'-px1', '5', '-px2', '5',
						'-pt1', '4', '-pt2', '4',
						'-wt1', '2', '-wt2', '2',
						'-wx2', '37',
						'-np1', np1, '-r1', r1,
						'-np2', '160', '-r2', '16'],
						stdout = subprocess.PIPE)

				stdout_val = proc.communicate()[0]
				print stdout_val

				# time execution
				end_time = time.time()
				run_time = floor(10*(end_time - start_time))/10

				# append execution time to measures.txt
				fi = open("measures.txt", "ab")
				fi.write('-time = ' + str(run_time) + '\n')
				fi.close()

				# create folder for results
				results_folder = '../results/' + method + '/table_rank1_1ps5x5x4_2ps5x5x4_2wx37/' + seq + '_s' + sigma + '_r1' + r1 + '_np1' + np1
				if not os.path.exists(results_folder):
					os.makedirs(results_folder)

				# move measures to folder
				shutil.move('measures.txt', results_folder + '/measures')

				# move images to folder
				results_prefixes = ['nisy', 'deno', 'bsic', 'diff', 'wei1','wei2']
				for prefix in results_prefixes:
					files = glob.glob(prefix + '*.png')
					for f in files:
						shutil.move(f, results_folder + '/' + f)



