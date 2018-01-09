#!/bin/bash
# Tune the algorithm's parameters

# noise levels
sigmas=(40)

# number of trials
ntrials=200

# test sequences
seqs=(\
derf-hd/park_joy \
derf-hd/speed_bag \
derf-hd/station2 \
derf-hd/sunflower \
derf-hd/tractor \
)

# seq folder
sf='/home/pariasm/remote/lime/denoising/data/'

output=${1:-"trials"}

export OMP_NUM_THREADS=1

# we assume that the binaries are in the same folder as the script
BIN=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
echo $BIN

for ((i=0; i < $ntrials; i++))
do
	# randomly draw parameters
	d=$(awk -v M=20   -v S=0    -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M - S) + S)}')
	t=$(awk -v M=6000 -v S=3000 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M - S) + S)}')
	b=$(awk -v M=4    -v S=2    -v s=$RANDOM 'BEGIN{srand(s); print rand()*(M - S) + S}')

	# format as string
	d=$(printf "%02d" $d)
	t=$(printf "%04d" $t)
	b=$(printf "%3.1f" $b)

	trialfolder=$(printf "$output/dsub%stau%sbt%s\n" $d $t $b)

	# run denoising on all sequences and compute average psnr
	mpsnr=0
	nseqs=${#seqs[@]}
	ff=70
	lf=85
	if [ ! -d $trialfolder ]
	then
		mkdir -p $trialfolder
		for seq in ${seqs[@]}
		do
#			echo  "$BIN/nldct -i ${sf}${seq}/%03d.png -f $ff -l $lf -sigma 40 -px2 0 -verbose 0 -dsub1 $d -tau1 $t -b1 $b"
			psnr=$($BIN/nldct -i ${sf}${seq}/%03d.png -f $ff -l $lf -sigma 40 -px2 0 -verbose 0 -dsub1 $d -tau1 $t -b1 $b)
			mpsnr=$(echo "$mpsnr + $psnr/$nseqs" | bc -l)
		done
	fi
	
	# store table entry with params and average psnr
	printf "$d $t $b %7.4f\n" $mpsnr >> $output/table
done
