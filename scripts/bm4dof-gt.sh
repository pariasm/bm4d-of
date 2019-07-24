#!/bin/bash
# Evals vnlm using ground truth

SEQ=$1 # sequence path
FFR=$2 # first frame
LFR=$3 # last frame
SIG=$4 # noise standard dev.
OUT=$5 # output folder
PRM=$6 # denoiser parameters

mkdir -p $OUT

# we assume that the binaries are in the same folder as the script
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )

# error checking {{{1
for i in $(seq $FFR $LFR);
do
	file=$(printf $SEQ $i)
	if [ ! -f $file ]
	then
		echo ERROR: $file not found
		exit 1
	fi
done

# add noise {{{1
$DIR/awgn $SIG $SEQ $OUT/"%03d.tif" $FFR $LFR

# denoise {{{1
$DIR/bm4dof.sh $OUT/"%03d.tif" $FFR $LFR $SIG $OUT $PRM

# compute psnr {{{1
$DIR/psnr $SEQ $OUT/deno_%03d.tif $FFR $LFR $OUT/measures-deno 1>/dev/null
$DIR/psnr $SEQ $OUT/bsic_%03d.tif $FFR $LFR $OUT/measures-bsic 1>/dev/null


# vim:set foldmethod=marker:
