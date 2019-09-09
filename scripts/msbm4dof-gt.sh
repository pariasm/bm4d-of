#!/bin/bash
if [ $# -lt 5 ]; then
	echo 'Multiscale VNLB video denoiser'
	echo ''
	echo ' msbm4dof-gt.sh'
	echo '    /path/to/frame-%03d.png # printf path to frames'
	echo '    first-frame             # first frame'
	echo '    last-frame              # last frame'
	echo '    sigma                   # noise'
	echo '    out-folder              # output folder'
	echo '    bm4dof-params             # (optional) parameters for VNLB'
	echo '    pyr-recomp              # (optional) recomposition ratio, e.g. 0.7'
	echo '    pyr-levels              # (optional) number of levels, e.g. 3'
	echo '    pyr-factor              # (optional) downscaling factor, e.g. 2'
	echo ''
	echo ' Example:'
	echo ' msbm4dof-gt.sh /path/to/frame-%03d.png 1 50 20 outs "-verbose 0" 3 2 0.7'
	exit 1
fi

SEQ=$1            # input sequence in printf format
FFR=$2            # first frame
LFR=$3            # last frame
SIG=$4            # noise std. dev
OUT=$5            # output folder
PRM=${6:-""}      # msbm4dof params
PYR_REC=${7:-0.7} # recomposition ratio
PYR_LVL=${8:--1}  # number of scales
PYR_DWN=${9:-2}   # downsampling factor

# error checking
for i in $(seq $FFR $LFR);
do
	file=$(printf $SEQ $i)
	if [ ! -f $file ]
	then
		echo ERROR: $file not found
		exit 1
	fi
done

# determine number of levels based on the image size
if [ $PYR_LVL -eq -1 ];
then
	PIXELS=$(imprintf "%N" $(printf $SEQ $FFR))
	printf -v PIXELS "%.f" "$PIXELS"
	echo $PIXELS
	if   [ ${PIXELS} -lt  500000 ]; then PYR_LVL=1
	elif [ ${PIXELS} -lt 2000000 ]; then PYR_LVL=2
	elif [ ${PIXELS} -lt 8000000 ]; then PYR_LVL=3
	else                                 PYR_LVL=4
	fi
fi

# add noise
$DIR/awgn $SIG $SEQ $OUT/"%03d.tif" $FFR $LFR

# denoise
$DIR/msbm4dof.sh "$OUT/ms${l}/%03d.tif" $FFR $LFR $LSIG $OUT/ms$l\
	"$PRM" $PYR_REC $PYR_LVL $PYR_DWN

# compute psnr
$DIR/psnr $SEQ $OUT/msdeno_%03d.tif $FFR $LFR $OUT/measures-msdeno 1>/dev/null
$DIR/psnr $SEQ $OUT/msbsic_%03d.tif $FFR $LFR $OUT/measures-msbsic 1>/dev/null
$DIR/psnr $SEQ $OUT/ms0/deno_%03d.tif $FFR $LFR $OUT/measures-deno 1>/dev/null
$DIR/psnr $SEQ $OUT/ms0/bsic_%03d.tif $FFR $LFR $OUT/measures-bsic 1>/dev/null
