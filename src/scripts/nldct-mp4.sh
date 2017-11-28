#!/bin/bash

#this script is used to execute the program from the IPOL demo
bindir=$1
vidin=$2
vidout=$3
sigma=$4
px=$5
pt=$6
wx=$7
wt=$8
np=$9

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
echo avconv -v error -i $vidin -f image2 i%04d.png
time avconv -v error -i $vidin -f image2 i%04d.png
echo

echo run nldct denoising
echo $bindir/nldct -i i%04d.png -f 1 -l $nframes -bsic d%04d.png -sigma $sigma \
     -px2 0 -px1 $px -pt1 $pt -wx1 $wx -wt1 $wt -np1 $np
time $bindir/nldct -i i%04d.png -f 1 -l $nframes -bsic d%04d.png -sigma $sigma \
     -px2 0 -px1 $px -pt1 $pt -wx1 $wx -wt1 $wt -np1 $np
echo

echo save output video as lossless mp4
echo avconv -v error -framerate $framerate -f image2 -i d%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $vidout
time avconv -v error -framerate $framerate -f image2 -i d%04d.png \
     -r $framerate -c:v libx264 -preset ultrafast -crf 0 $vidout
echo

echo remove pngs
rm -R {i,d}????.png
