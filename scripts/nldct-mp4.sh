#!/bin/bash

#this script is used to execute the program from the IPOL demo
bindir=$1
vidin=$2
denoout=$3
nisyout=$4
sigma=$5
px=$6
pt=$7
wx=$8
wt=$9
np=${10}
pxlfmt=${11}

#extract info from video
info=`avprobe -v error -show_streams  $vidin`
info="${info#*codec_type=video}"
echo

width=`echo ${info#*width=}| cut -d' ' -f 1` 
height=`echo ${info#*height=}| cut -d' ' -f 1` 
framerate=`echo ${info#*avg_frame_rate=}| cut -d' ' -f 1 | sed 's/\/1//'`
nframes=`echo ${info#*nb_frames=}| cut -d' ' -f 1`
size=${width}x${height}

if [ "$framerate" == "0/0" ] ; then
  echo "Error reading the frame rate of the video (default 30)"
  framerate="30"
fi

echo input parameters and video info
echo " "sigma: $sigma
echo " "psz: $px x $pt
echo " "wsz: $wx x $wt
echo " "np:  $np
echo " "input video size: $width x $height x $nframes @ $framerate fps
echo

echo extract frames from input video
if [ $pxlfmt == "gray" ]
then
	echo avconv -v error -i $vidin -f image2 -vf format=gray i%04d.png
	time avconv -v error -i $vidin -f image2 -vf format=gray i%04d.png

	for i in $(seq 1 $nframes)
	do
		plambda i$(printf %04d $i).png "x[0]" -o i$(printf %04d $i).png
	done
else
	echo avconv -v error -i $vidin -f image2 i%04d.png
	time avconv -v error -i $vidin -f image2 i%04d.png
fi
echo

echo run nldct denoising
export OMP_NUM_THREADS=8 # set max number of threads
echo $bindir/nldct -i i%04d.png -f 1 -l $nframes -bsic d%04d.png -nisy n%04d.png \
	  -sigma $sigma -px2 0 -px1 $px -pt1 $pt -wx1 $wx -wt1 $wt -np1 $np
time $bindir/nldct -i i%04d.png -f 1 -l $nframes -bsic d%04d.png -nisy n%04d.png \
	  -sigma $sigma -px2 0 -px1 $px -pt1 $pt -wx1 $wx -wt1 $wt -np1 $np
echo

echo save output video as lossless mp4
echo avconv -y -v error -framerate $framerate -f image2 -i d%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $denoout
time avconv -y -v error -framerate $framerate -f image2 -i d%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $denoout

echo avconv -y -v error -framerate $framerate -f image2 -i n%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $nisyout
time avconv -y -v error -framerate $framerate -f image2 -i n%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $nisyout

if [ $pxlfmt == "gray" ]
then
	# overwrite input video (if it was RGB, it will be overwritten by
	# a grayscale video)
	echo avconv -y -v error -framerate $framerate -f image2 -i i%04d.png \
	     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $vidin
	time avconv -y -v error -framerate $framerate -f image2 -i i%04d.png \
	     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $vidin
fi
echo

nkeep=19
echo keep $nkeep frames at the beginning of the sequence
f1=1
f2=$((nkeep > nframes ? nframes : nkeep))

echo remove pngs
for ((i=$nkeep+1; i<=$nframes; i++))
do
	rm -R {i,d,n}$(printf "%04d" $i).png
done
