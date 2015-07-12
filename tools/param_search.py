import subprocess
import shutil
import time
from math import floor
import os
import glob

sequences = [
'Army']
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

base_dir = '/data/pablo/ipol/local/video_nlbayes/input/'

sigma_str  = ['10', '20', '40']
method = 'vnlbayes'

for isigma in range(0,len(sigma_str)):
	sigma = sigma_str[isigma]
	for seq in sequences:
		for xrad in spacial_radii2[isigma]:
			for psz in patch_size2[isigma]:

				print "Processing " + seq \
					+ ": sigma = " + sigma \
					+ "  xrad = " + xrad \
					+ "  psz  = " + psz
	
				# time execution
				start_time = time.time()
	
				# launch process
				proc = subprocess.Popen(['bin/' + method,
					'-i', base_dir + seq + '/i%04d.png',
					'-f', '0', '-l', '7',
					'-sigma', sigma,
					'-wx2', xrad,
					'-ps2', psz,
					'-wt1', '2', '-wt2', '2',
					'-flat-area1', '0', '-flat-area2', '0'],
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
				results_folder = '../results/' \
						+ method + '/prms_step2_offset/' + seq \
						+ '_' + 's' + sigma \
						+ '_' + 'tr' + '2' \
						+ '_' + 'xr' + xrad \
						+ '_' + 'ps' + psz
				if not os.path.exists(results_folder):
					os.makedirs(results_folder)
	
				# move measures to folder
				shutil.move('measures.txt',
						results_folder + '/measures')
	
				# move images to folder
				results_prefixes = ['nisy', 'deno', 'bsic', 'diff']
				for prefix in results_prefixes:
					files = glob.glob(prefix + '*.png')
					for f in files:
						shutil.move(f, results_folder + '/' + f)



