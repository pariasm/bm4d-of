/*
 * Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
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
 * @file lib_videot.cpp
 * @brief Spetializations of template functions defined in lib_videot.hpp
 *
 * @author Pablo Arias <pariasm@gmail.com>
 **/

#include "lib_videot.hpp"

//#include <fftw3.h>
#include "lib_image.h"

#include <cstdio>
#include <cstdlib> // EXIT_FAILURE


template <> 
void Video<float>::loadVideo(
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
			throw std::runtime_error("Video<T>::loadVideo: loading of " + 
					std::string(filename) + " failed");

		//! set size
		sz.width    = frameSize.width;
		sz.height   = frameSize.height;
		sz.channels = frameSize.nChannels;
		sz.frames   = (i_lastFrame - i_firstFrame + 1)/i_frameStep;
		sz.update_fields();
		
		//! allocate
		data.resize(sz.whcf);

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

		loadImage(filename, frame, frameSize);
		if (loadImage(filename, frame, frameSize) == EXIT_FAILURE)
			throw std::runtime_error("Video<T>::loadVideo: loading of " + 
					std::string(filename) + " failed");

		if (frameSize.width != sz.width ||
		    frameSize.height != sz.height ||
		    frameSize.nChannels != sz.channels)
			throw std::runtime_error("Video<T>::loadVideo: size of " + 
					std::string(filename) + " does not match video size");

		//! copy data in first frame // FIXME: avoidable memcpy
		p_data = std::copy(frame.begin(), frame.end(), p_data);
	}

	return;
}
	
template <> 
void Video<float>::saveVideo(
	const std::string i_pathToFiles
,	unsigned i_firstFrame
,	unsigned i_frameStep
,	float i_pmin
,	float i_pmax
) const {

	//! size structure for frames
	ImageSize frameSize;
	frameSize.width     = sz.width;
	frameSize.height    = sz.height;
	frameSize.nChannels = sz.channels;
	frameSize.wh        = sz.wh;
	frameSize.whc       = sz.whc;

	unsigned lastFrame = i_firstFrame + sz.frames * i_frameStep;
	std::vector<float> frame(sz.whc);
	std::vector<float>::const_iterator p_data = data.begin();
	for (unsigned f = i_firstFrame; f < lastFrame; f += i_frameStep, p_data += sz.whc)
	{
		//! copy data to frame // FIXME: avoidable memcpy
		std::copy(p_data, p_data + sz.whc, frame.begin());

		//! file name
		char filename[1024];
		std::sprintf(filename, i_pathToFiles.c_str(), f);
		
		if (saveImage(filename, frame, frameSize, i_pmin, i_pmax) == EXIT_FAILURE)
			throw std::runtime_error("Video<T>::saveVideo: writing of " + 
					std::string(filename) + " failed");
	}

	return;
}

template <> 
void Video<float>::saveVideoAscii(
	const std::string i_prefix
,	unsigned i_firstFrame
,	unsigned i_frameStep
) const {

	char channel[32], frame[32];

	int C = sz.channels;
	int F = sz.frames;
	int W = sz.width;
	int H = sz.height;

	for (size_t c = 0; c < C; ++c)
	for (size_t f = 0; f < F; ++f)
	{
		if (C > 1) std::sprintf(channel,"_ch%d",(int)c); else channel[0] = 0;
		if (F > 1) std::sprintf(frame  ,"_%03d",(int)f); else frame[0] = 0;

		std::string filename = i_prefix + channel + frame + ".asc";
		std::FILE *const nfile = std::fopen(filename.c_str(),"w");

		unsigned idx = sz.index(0,0,f,c);
		for (size_t y = 0; y < H; ++y)
		{
			for (size_t x = 0; x < W; ++x, ++idx)
				std::fprintf(nfile,"%.16g " , (double)data[idx]);
			std::fprintf(nfile, "\n");
		}

		std::fclose(nfile);
	}
	return;
}

