#!/bin/bash
if [ $# -lt 5 ]; then
	echo 'Multiscale VNLB video denoiser'
	echo ''
	echo ' vnlb-ms.sh'
	echo '    /path/to/frame-%03d.png # printf path to frames'
	echo '    first-frame             # first frame'
	echo '    last-frame              # last frame'
	echo '    sigma                   # noise'
	echo '    out-folder              # output folder'
	echo '    vnlb-params             # (optional) parameters for VNLB'
	echo '    pyr-recomp              # (optional) recomposition ratio, e.g. 0.7'
	echo '    pyr-levels              # (optional) number of levels, e.g. 3'
	echo '    pyr-factor              # (optional) downscaling factor, e.g. 2'
	echo ''
	echo ' Example:'
	echo ' vnlb-ms.sh /path/to/frame-%03d.png 1 50 20 outs "-verbose 0" 3 2 0.7'
	exit 1
fi

SEQ=$1            # input sequence in printf format
FFR=$2            # first frame
LFR=$3            # last frame
SIG=$4            # noise std. dev
OUT=$5            # output folder
PRM=${6:-""}      # denoiser params
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

echo "Scales: $PYR_LVL"

# we assume that the binaries are in the same folder as the script
DIR=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
DIRMS=$DIR"/ms-lanczos3"

# create folders for scales
for ((l=PYR_LVL-1; l>=0; --l)); do mkdir -p $OUT/ms$l; done

for ((i=FFR; i<=LFR; ++i))
do
	$DIRMS/lanczos3_decompose.m $(printf $SEQ $i) "$OUT/ms"\
		$PYR_LVL $(printf "/%03d.tif" $i)
done

for ((l=PYR_LVL-1; l>=0; --l))
do
	LSIG=$(bc <<< "scale=2; $SIG / ${PYR_DWN}^$l")
	$DIR/bm4dof.sh "$OUT/ms${l}/%03d.tif" $FFR $LFR $LSIG $OUT/ms$l "$PRM"
done

for ((i=FFR; i<=LFR; ++i))
do
	$DIRMS/lanczos3_recompose.m $(printf "$OUT/msbsic_%03d.tif" $i) "$OUT/ms" $PYR_LVL $(printf "/bsic_%03d.tif" $i) $PYR_REC
	$DIRMS/lanczos3_recompose.m $(printf "$OUT/msdeno_%03d.tif" $i) "$OUT/ms" $PYR_LVL $(printf "/deno_%03d.tif" $i) $PYR_REC
done
