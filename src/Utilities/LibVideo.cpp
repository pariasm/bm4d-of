/*
 * Copyright (c) 2013, Pablo Arias <pariasm@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file LibVideo.cpp
 * @brief Video container with basic functionalities
 *
 * @author Pablo Arias <pablo.arias@gmail.com>
 **/

#include "LibVideo.h"
#include "LibImages.h"
#include <iostream>
#include <string>
#include <cstdio> // sprintf
#include <cstdlib> // EXIT_FAILURE/SUCCESS
#include <omp.h>

Video_f32::Video_f32(void)
	: width(0)
	, height(0)
	, frames(0)
	, channels(0)
	, wh(0)
	, whc(0)
	, whcf(0)
	, data(0)
{ }

Video_f32::Video_f32(const Video_f32& i_in)
	: width(i_in.width)
	, height(i_in.height)
	, frames(i_in.frames)
	, channels(i_in.channels)
	, wh(i_in.wh)
	, whc(i_in.whc)
	, whcf(i_in.whcf)
	, data(i_in.data)
{ }

Video_f32::Video_f32(
	const std::string i_pathToFiles
,	unsigned i_firstFrame
,	unsigned i_lastFrame
,	unsigned i_frameStep
)
	: width(0)
	, height(0)
	, frames(0)
	, channels(0)
	, wh(0)
	, whc(0)
	, whcf(0)
	, data(0)
{
	int loaded_ok = loadVideo(i_pathToFiles,
			                    i_firstFrame, i_lastFrame, i_frameStep);
	// TODO what to do if not loaded_ok?
}

Video_f32::Video_f32(
	unsigned i_width
,	unsigned i_height
,	unsigned i_frames
,	unsigned i_channels
)
	: width(i_width)
	, height(i_height)
	, frames(i_frames)
	, channels(i_channels)
	, wh(width * height)
	, whc(width * height * channels)
	, whcf(width * height * channels * frames)
	, data(whcf)
{ }

Video_f32::Video_f32(const VideoSize& i_size)
	: width(i_size.width)
	, height(i_size.height)
	, frames(i_size.frames)
	, channels(i_size.channels)
	, wh(width * height)
	, whc(width * height * channels)
	, whcf(width * height * channels * frames)
	, data(whcf)
{ }

void Video_f32::clear(void)
{
	width = 0;
	height = 0;
	frames = 0;
	channels = 0;
	wh = 0;
	whc = 0;
	whcf = 0;
	data.clear();
}

void Video_f32::resize(
	unsigned i_width
,	unsigned i_height
,	unsigned i_frames
,	unsigned i_channels
){
	if (width != i_width ||
	    height != i_height ||
	    frames != i_frames ||
	    channels != i_channels ) 
	{
		clear();

		width = i_width;
		height = i_height;
		frames = i_frames;
		channels = i_channels;
		wh = width * height;
		whc = wh * channels;
		whcf = whc * frames;
		data.resize(whcf);
	}
}

void Video_f32::resize(const VideoSize& i_size)
{
	resize(i_size.width, i_size.height, i_size.frames, i_size.channels);
}


int Video_f32::loadVideo(
    const std::string i_pathToFiles
,   unsigned i_firstFrame
,   unsigned i_lastFrame
,   unsigned i_frameStep
){
	clear();

	//! open first frame and allocate memory
	std::vector<float>::iterator p_data;
	{
		char filename[1024];
		std::sprintf(filename, i_pathToFiles.c_str(), i_firstFrame);

		ImageSize frameSize;
		std::vector<float> frame;
		if (loadImage(filename, frame, frameSize) == EXIT_FAILURE)
			return EXIT_FAILURE;

		//! set size
		width    = frameSize.width;
		height   = frameSize.height;
		channels = frameSize.nChannels;
		frames   = (i_lastFrame - i_firstFrame + 1)/i_frameStep;
		wh       = width * height;
		whc      = width * height * channels;
		whcf     = width * height * channels * frames;
		
		//! allocate
		data.resize(whcf);

		//! copy data in first frame
		p_data = std::copy(frame.begin(), frame.end(), data.begin());
	}
	
	//! load rest of frames
	for (unsigned f = i_firstFrame + i_frameStep; f <= i_lastFrame; f += i_frameStep)
	{
		char filename[1024];
		std::sprintf(filename, i_pathToFiles.c_str(), f);

		ImageSize frameSize;
		std::vector<float> frame;

		if (loadImage(filename, frame, frameSize) == EXIT_FAILURE)
		{
			std::cerr << "Error @ Video_f32::loadVideo: Frame "
			          << f << " cannot be loaded" << std::endl;
			return EXIT_FAILURE;
		}

		if (frameSize.width != width ||
		    frameSize.height != height ||
		    frameSize.nChannels != channels)
		{
			std::cerr << "Error @ Video_f32::loadVideo: Size of frame "
			          << f << " does not match video size" << std::endl;
			return EXIT_FAILURE;
		}

		//! copy data in first frame // FIXME: avoidable memcpy
		p_data = std::copy(frame.begin(), frame.end(), p_data);
	}

	return EXIT_SUCCESS;
}

int Video_f32::saveVideo(
	const std::string i_pathToFiles
,	unsigned i_firstFrame
,	unsigned i_frameStep
,	float i_pmin
,	float i_pmax
) const {

	//! size structure for frames
	ImageSize frameSize;
	frameSize.width     = width;
	frameSize.height    = height;
	frameSize.nChannels = channels;
	frameSize.wh        = width * height;
	frameSize.whc       = width * height * channels;

	unsigned lastFrame = i_firstFrame + frames * i_frameStep;
	std::vector<float> frame(whc);
	std::vector<float>::const_iterator p_data = data.begin();
	for (unsigned f = i_firstFrame; f < lastFrame; f += i_frameStep, p_data += whc)
	{
		//! copy data to frame // FIXME: avoidable memcpy
		std::copy(p_data, p_data + whc, frame.begin());

		//! file name
		char filename[1024];
		std::sprintf(filename, i_pathToFiles.c_str(), f);
		
		if (saveImage(filename, frame, frameSize, i_pmin, i_pmax) == EXIT_FAILURE)
		{
			std::cerr << "Error @ Video_f32::saveVideo: Frame "
			          << f << " cannot be saved" << std::endl;
			return EXIT_FAILURE;
		}
	}

	return EXIT_SUCCESS;
}

int Video_f32::saveVideoAscii(
	const std::string i_prefix
,	unsigned i_firstFrame
,	unsigned i_frameStep
) const {
	char channel[32], frame[32];

	int C = channels;
	int F = frames;
	int W = width;
	int H = height;

	for (size_t c = 0; c < C; ++c)
	for (size_t f = 0; f < F; ++f)
	{
		if (C > 1) std::sprintf(channel,"_ch%d",(int)c); else channel[0] = 0;
		if (F > 1) std::sprintf(frame  ,"_%03d",(int)f); else frame[0] = 0;

		std::string filename = i_prefix + channel + frame + ".asc";
		std::FILE *const nfile = std::fopen(filename.c_str(),"w");

		unsigned idx = index(0,0,f,c);
		for (size_t y = 0; y < H; ++y)
		{
			for (size_t x = 0; x < W; ++x, ++idx)
				std::fprintf(nfile,"%.16g " , (double)data[idx]);
			std::fprintf(nfile, "\n");
		}

		std::fclose(nfile);
	}
	return EXIT_SUCCESS;
}

