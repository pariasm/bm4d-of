#!/bin/bash
# Tune the algorithm's parameters

# noise levels
sigmas=(10 20 40)

# number of trials
ntrials=200

# test sequences
seqs=(\
gray_2.4_alley.png \
gray_2.4_computer.png \
gray_2.4_gardens.png \
gray_2.4_girl_crop.png \
gray_2.4_man2_crop.png \
gray_2.4_street1_crop.png \
gray_2.4_traffic.png \
gray_2.4_yard.png \
)
#seqs=(\
#gray_2.4_alley.png \
#gray_2.4_building1.png \
#gray_2.4_computer.png \
#gray_2.4_flowers1.png \
#gray_2.4_gardens.png \
#gray_2.4_girl_crop.png \
#gray_2.4_hallway.png \
#gray_2.4_street1_crop.png \
#gray_2.4_traffic.png \
#gray_2.4_trees.png \
#gray_2.4_yard.png \
#)

# seq folder
sf='/home/pariasm/denoising/data/images/tune/gray/'

output=${1:-"trials"}

export OMP_NUM_THREADS=1

# we assume that the binaries are in the same folder as the script
BIN=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
echo $BIN

for ((i=0; i < $ntrials; i++))
do
	# noise level
	r=$(awk -v M=3 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*M)}')
	s=${sigmas[$r]}

	# randomly draw parameters
	t1=$(awk -v M=4000 -v S=2000 -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M - S + 1) + S)}')
	t2=$(awk -v M=600  -v S=200  -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M - S + 1) + S)}')
	n1=$(awk -v M=6    -v S=1    -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M - S + 1) + S)}')
	n2=$(awk -v M=6    -v S=1    -v s=$RANDOM 'BEGIN{srand(s); print int(rand()*(M - S + 1) + S)}')
	b1=$(awk -v M=3.0  -v S=0.8  -v s=$RANDOM 'BEGIN{srand(s); print rand()*(M - S) + S}')
	b2=$(awk -v M=2.0  -v S=0.8  -v s=$RANDOM 'BEGIN{srand(s); print rand()*(M - S) + S}')
	n1=$(plambda -c "2 $n1 pow")
	n2=$(plambda -c "2 $n2 pow")
#	w1=
#	w2=
#	bm1=
#	bm2=

	# format as string
	n1=$(printf "%02.0f"  $n1)
	n2=$(printf "%02.0f"  $n2)
	b1=$(printf "%3.1f" $b1)
	b2=$(printf "%3.1f" $b2)
	t1=$(printf "%04d"  $t1)
	t2=$(printf "%04d"  $t2)

	trial=$(printf "$output/%d-n%s-b%s-t%s-n%s-b%s-t%s" $s $n1 $b1 $t1 $n2 $b2 $t2)

	# run denoising on all sequences and compute average psnr
	mse_bsic=0
	mse_deno=0
	nseqs=${#seqs[@]}
	if [ ! -d $trial ]
	then
		echo $trial
		mkdir -p $trial
		for seq in ${seqs[@]}
		do
#			echo  "$BIN/nldct -i $sf$seq -sigma $s -order-inv 1 -verbose 0 -np1 $n1 -b1 $b1 -t1 $t1 -np2 $n2 -b2 $b2 -t2 $t2"
			mse=($($BIN/nldct -i $sf$seq -sigma $s -order-inv 1 -verbose 0 -np1 $n1 -b1 $b1 -t1 $t1 -np2 $n2 -b2 $b2 -t2 $t2))
			echo ${mse[0]} ${mse[1]}
			mse_deno=$(echo "$mse_deno + ${mse[0]}/$nseqs" | bc -l)
			mse_bsic=$(echo "$mse_bsic + ${mse[1]}/$nseqs" | bc -l)
		done
	fi
	
	# store table entry with params and average psnr
	printf "$s $n1 $b1 $t1 $n2 $b2 $t2 %7.4f %7.4f\n" $mse_bsic $mse_deno >> $output/table
done
