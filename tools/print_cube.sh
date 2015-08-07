#! /bin/bash
# Prints a table with the results of a sequence given as input
#
# Tables are as follows:
#
# px \ wx
#    13 25 37 49
# 3
# 5
# 7
# 9
#
# Usage
# ./print_table.sh SEQUENCE SIGMA [MEASURE]
#
# SEQUENCE can be Army, DogDance, Evergreen, Mequon, Walking
# SIGMA    can be 10, 20, 40
# MEASURE  can be [PSNR_final], PSNR_basic, RMSE_final, RMSE_basic, time

# default value for MEASURE argument
measure=${3:-PSNR_final}

# declare arrays to read tables
nps=('080' '160' '320' '640')
ranks=('08' '16' '32' '64')
pxs=('3' '5' '7' '9')

for px in ${pxs[@]};
do
	echo "% table for patch size $px"
	for r in ${ranks[@]};
	do
		for ((i=0; i < ${#nps[@]}; i++));
		do
			# read specified measure from measures
			t=`cat $1_s$2_1px${px}_1r${r}_1np${nps[$i]}/measures | grep $measure | sed "s/-$measure = //"` 

			# round to 2 decimals and store in vals array
			val[$i]=`echo "scale=2; (10^2*$t + 0.5) / 10^2" | bc`
		done
		
		# printf array separated by tabs into vals
		printf -v vals "%s\t" "${val[@]}"

		# remove last character (it's a tab)
		vals=${vals%?}

		# printf vals with a newline end character to stdout
		echo "$vals"
	done
done

