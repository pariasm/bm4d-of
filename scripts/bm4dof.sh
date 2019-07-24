#!/bin/bash
# Helper script for the video NL-Bayes denoiser. Given
# a noisy video it computes the tvl1 optical flow and then
# calls the denoiser.

SEQ=$1 # noisy sequence path, in printf format, e.g. /my/sequence/frame-%02d.tif
FFR=$2 # first frame
LFR=$3 # last frame
SIG=$4 # noise standard dev.
OUT=$5 # output folder
PRM=$6 # denoiser parameters (between quotes, e.g. "-px1 5 -pt1 4 -verbose 0" 

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

# compute optical flow {{{1
FLOWOPTS=""
if [ $LFR -gt $FFR ]; then
	$DIR/tvl1flow-seq.sh $SEQ $FFR $LFR $OUT
	FLOWOPTS="-fof $OUT/tvl1_%03d_f.flo -bof $OUT/tvl1_%03d_b.flo"
fi

# run denoising {{{1

# run first step
echo "$DIR/bm4d-of \
 -i $SEQ -f $FFR -l $LFR -sigma $SIG $FLOWOPTS \
 -bsic $OUT/bsic_%03d.tif -px2 0 $PRM"

$DIR/bm4d-of \
 -i $SEQ -f $FFR -l $LFR -sigma $SIG $FLOWOPTS \
 -bsic $OUT/bsic_%03d.tif -px2 0 $PRM

# run second step
echo "$DIR/bm4d-of \
 -i $SEQ -f $FFR -l $LFR -sigma $SIG -b $OUT/bsic_%03d.tif $FLOWOPTS \
 -deno $OUT/deno_%03d.tif -px1 0 $PRM"

$DIR/bm4d-of \
 -i $SEQ -f $FFR -l $LFR -sigma $SIG -b $OUT/bsic_%03d.tif $FLOWOPTS \
 -deno $OUT/deno_%03d.tif -px1 0 $PRM

# vim:set foldmethod=marker:
