#!/bin/bash
# Tune the algorithm's parameters

# noise levels
sigmas=(10 20 40)

# number of trials
ntrials=2000

# test sequences
seqs=(\
derf-hd/aspen \
derf-hd/old_town_cross \
derf-hd/rush_field_cuts \
derf-hd/rush_hour \
)

# seq folder
sf='/home/pariasm/remote/lime/denoising/data/'
sff='/home/pariasm/remote/lime/denoising/projects/nldct/results/'

output=${1:-"trials"}

export OMP_NUM_THREADS=1

# we assume that the binaries are in the same folder as the script
BIN=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
#echo $BIN

for ((i=0; i < $ntrials; i++))
do
	# pick a noise level
	r=$(awk -v M=3 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*M)}')
	s=${sigmas[$r]}

	# randomly draw parameters
	n1=$(awk -v M=7   -v S=1  -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M - S + 1) + S)}')
	n2=$(awk -v M=7   -v S=1  -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M - S + 1) + S)}')
	b1=$(awk -v M=4.0 -v S=.5 -v s=$RANDOM 'BEGIN{srand(s); print rand()*(M - S) + S}')
	b2=1
	t1=$(awk -v M=100. -v S=0. -v s=$RANDOM 'BEGIN{srand(s); print rand()*(M - S) + S}')
	t2=$(awk -v M=60.  -v S=0. -v s=$RANDOM 'BEGIN{srand(s); print rand()*(M - S) + S}')

	n1=$(plambda -c "2 $n1 pow")
	n2=$(plambda -c "2 $n2 pow")

	# format as string
	s=$(printf "%02d" $s)
	n1=$(printf "%03d" $n1)
	n2=$(printf "%03d" $n2)
	b1=$(printf "%3.1f" $b1)
	b2=$(printf "%3.1f" $b2)
	t1=$(printf "%04.1f" $t1)
	t2=$(printf "%04.1f" $t2)

	trialfolder="$output/s$s-$n1-$n2-$b1-$b2-$t1-$t2"

	# run denoising on all sequences and compute average psnr
	mse_bsic=0
	mse_deno=0
	nseqs=${#seqs[@]}
	ff=231
	lf=250
	if [ ! -d $trialfolder ]
	then
		mkdir -p $trialfolder
		for seq in ${seqs[@]}
		do
			prms="-px2 10 -pt2 2 -wx2 21 -wt2 4 -np2 $n2 -b2 $b2 -t2 $t2"
			prms="-px1 10 -pt1 2 -wx1 21 -wt1 4 -np1 $n1 -b1 $b1 -t1 $t1 $prms"
			flow="-bof ${sff}${seq}/s$s/%03d_b.flo"
			flow="-fof ${sff}${seq}/s$s/%03d_f.flo $flow"
			echo   "$BIN/nldct -i ${sf}${seq}/%03d.png -f $ff -l $lf $flow -sigma $s -verbose 0 $prms"
			mse=($($BIN/nldct -i ${sf}${seq}/%03d.png -f $ff -l $lf $flow -sigma $s -verbose 0 $prms))
			echo ${mse[0]} ${mse[1]}
			mse_deno=$(echo "$mse_deno + ${mse[0]}/$nseqs" | bc -l)
			mse_bsic=$(echo "$mse_bsic + ${mse[1]}/$nseqs" | bc -l)
		done
	fi
	
	# store table entry with params and average psnr
	printf "$s $n1 $b1 $n2 $b2 $t1 $t2 %7.4f %7.4f\n" $mse_bsic $mse_deno >> $output/table
done
