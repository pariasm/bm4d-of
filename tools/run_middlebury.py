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

patch_temp = ['4','3','2']
sigma_str  = ['40','20','10']
method = 'vnlbayes'

for pt in patch_temp:
	for sigma in sigma_str:
		for iseq in range(0,len(sequences)):
			seq = sequences[iseq]
			print "Processing " + seq + " with pt = " + pt

			# time execution
			start_time = time.time()

			# launch process
			proc = subprocess.Popen(['bin/' + method,
					'-i', base_dir + seq + '/i%04d.png',
					'-f', '0', '-l', '7',
					'-sigma', sigma,
					'-px1', '5', '-px2', '5',
					'-pt1', pt , '-pt2', pt ,
					'-wt1', '2', '-wt2', '2',
					'-wx2', '37'],
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
			results_folder = '../results/' + method + '/table_1ps5_2ps5_2wx37/' + seq + '_' + 's' + sigma + '_' + 'pt' + pt
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



