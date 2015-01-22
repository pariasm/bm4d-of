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
			vid.nFrames, vid.nChannels, vid.height, vid.width);
}

int main(int argc, char **argv)
{
	//! Check if there is the right call for the algorithm
	if (argc != 5)
	{
		fprintf(stdout, "Usage: %s path-to-frames first-frame last-frame frame-step\n", argv[0]);
		return EXIT_FAILURE;
	}

	//! Get CMD line inputs
	const char* i_video_path =    argv[1] ;
	const int i_firstFrame = atoi(argv[2]);
	const int i_lastFrame  = atoi(argv[3]);
	const int i_frameStep  = atoi(argv[4]);

	//! Create empty video
	Video_f32 vid1;
	print_video_size("empty video", vid1);

	//! Save video
	vid1.saveVideo("/tmp/vid1_empty_%02d.png",i_firstFrame, i_frameStep);

	//! Load a video through loadVideo 
	if (vid1.loadVideo(i_video_path,
			i_firstFrame, i_lastFrame, i_frameStep) == EXIT_FAILURE)
	{
		fprintf(stderr, "Exiting. Failed to load video.\n");
		return EXIT_FAILURE;
	}
	print_video_size("loaded video1", vid1);

	// accessing pixels through coordinates
	for (int f =  0; f < vid1.nFrames; f += 2)
	for (int y = 10; y < vid1.height - 10; y++)
	for (int x = 10; x < vid1.width  - 10; x++)
		vid1(x,y,f,0) = 250;

	// accessing pixels through indices
	for (int i = 0; i < vid1.whcf; i += 17) vid1(i) = 0;

	//! Save video
	vid1.saveVideo("/tmp/vid1_%02d.png",i_firstFrame, i_frameStep);

	//! Load a video through constructor
	Video_f32 vid2(i_video_path, i_firstFrame, i_lastFrame, i_frameStep);
	print_video_size("loaded video2", vid2);

	//! Save video
	vid2.saveVideo("/tmp/vid2_%02d.png",i_firstFrame, i_frameStep);

	//! Create an empty copy
	Video_f32 vid3(vid1);
	print_video_size("copied video3", vid3);

	//! Save
	vid2.saveVideo("/tmp/vid3_%02d.png",i_firstFrame, i_frameStep);

	return EXIT_SUCCESS;
}
