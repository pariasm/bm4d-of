/*
 * Copyright (c) 2015, Pablo Arias <pariasm@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

#include <string>
#include <sstream>

#include "Utilities.h"
#include "LibVideo.h"

/**
 * @file   testLibVideo.cpp
 * @brief  Executable file to test libVideo
 *
 *
 *
 * @author PABLO ARIAS  <pariasm@gmail.com>
 **/



void print_video_size(const std::string& name, const Video_f32& vid)
{
	printf("%s", name.c_str());
	printf("\tFrames: %d - Channels: %d - Height: %d - Width %d\n", 
			vid.frames, vid.channels, vid.height, vid.width);
}

int main(int argc, char **argv)
{
	//! Check if there is the right call for the algorithm
	if (argc != 6)
	{
		fprintf(stdout, "Usage: %s path-to-frames first-frame last-frame frame-step sigma\n", argv[0]);
		return EXIT_FAILURE;
	}

	//! Get command line inputs
	const char* i_video_path =    argv[1] ;
	const int i_firstFrame = atoi(argv[2]);
	const int i_lastFrame  = atoi(argv[3]);
	const int i_frameStep  = atoi(argv[4]);
	const int i_sigma      = atoi(argv[5]);

	//! Create empty video
	Video_f32 vid1;
	print_video_size("empty video", vid1);

	//! Save video / should do nothing!
	vid1.saveVideo("/tmp/vid1_empty_%02d.png",i_firstFrame, i_frameStep);

	//! Load a video through loadVideo 
	vid1.loadVideo(i_video_path, i_firstFrame, i_lastFrame, i_frameStep);
	print_video_size("loaded video1", vid1);


	/*/! Accessing pixels through coordinates
	for (int f =  0; f < vid1.frames; f += 2)
	for (int y = 10; y < vid1.height - 10; y++)
	for (int x = 10; x < vid1.width  - 10; x++)
		vid1(x,y,f,0) = 250;

	//! Accessing pixels through indices
	for (int i = 0; i < vid1.whcf; i += 17) vid1(i) = 0;//*/

	//! Save video
	vid1.saveVideo("/tmp/vid1_%02d.png", i_firstFrame, i_frameStep);
//	vid1.saveVideoAscii("/tmp/vid1"    , i_firstFrame, i_frameStep);//*/

	//! Subdivide into smaller videos with border
	unsigned nvids = 8;
	unsigned border = 10;
	std::vector<Video_f32> subvids1;
	VideoUtils::subDivide(vid1, subvids1, border, nvids); 

	for (int n = 0; n < subvids1.size(); n++)
	{
		char name[1024];
		sprintf(name, "/tmp/vid1_sub%02d_%%02d.png", n);
		subvids1[n].saveVideo(name, i_firstFrame, i_frameStep);
	}

	//! Join again into large video, removing border
	Video_f32 vid3(vid1.size());
	VideoUtils::subBuild (subvids1, vid3, border); 

	float psnr_13, rmse_13;
	VideoUtils::computePSNR(vid1, vid3, psnr_13, rmse_13);
	printf("After joining the subvideos... RMSE: %f - PSNR: %f\n", rmse_13, psnr_13);//*/


	/*/! Pad video by symmetry
	Video_f32 vid1_sym;
	VideoUtils::addBorder(vid1, vid1_sym, 4, true);
	print_video_size("video with border added by symmetrizing", vid1_sym);
	vid1_sym.saveVideo("/tmp/vid1_sym_%02d.png", i_firstFrame, i_frameStep);//*/

	/*/! Change colorspace
	VideoUtils::transformColorSpace(vid1, true);
	vid1.saveVideoAscii("/tmp/vid1_yuv", i_firstFrame, i_frameStep);

	VideoUtils::transformColorSpace(vid1, false);//*/

	/*/! Add noise
	Video_f32 vid1_noise;
	VideoUtils::addNoise(vid1, vid1_noise, i_sigma, true);

	//! Save video
	vid1_noise.saveVideo("/tmp/vid1_noise_%02d.png",i_firstFrame, i_frameStep);//*/

	/*/! Compute PSNR
	float psnr, rmse;
	VideoUtils::computePSNR(vid1, vid1_noise, psnr, rmse);
	printf("Computed RMSE: %f - PSNR: %f\n", rmse, psnr);//*/

	/*/! Compute difference
	Video_f32 diff;
	VideoUtils::computeDiff(vid1, vid1_noise, diff, i_sigma);
	print_video_size("difference", diff);

	//! Save difference video
	diff.saveVideo("/tmp/diff_%02d.png",i_firstFrame, i_frameStep);//*/

	/*/! Load a video through constructor
	Video_f32 vid2(i_video_path, i_firstFrame, i_lastFrame, i_frameStep);
	print_video_size("loaded video2", vid2);

	//! Save video
	vid2.saveVideo("/tmp/vid2_%02d.png",i_firstFrame, i_frameStep);

	//! Create an empty copy
	Video_f32 vid3(vid1);
	print_video_size("copied video3", vid3);

	//! Save
	vid2.saveVideo("/tmp/vid3_%02d.png",i_firstFrame, i_frameStep);*/

	return EXIT_SUCCESS;
}
