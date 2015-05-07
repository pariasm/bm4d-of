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
#include <cassert>

#include <string>
#include <sstream>

#include "Utilities.h"
#include "LibVideoT.hpp"


/**
 * @file   testLibVideoT.cpp
 * @brief  Executable file to test libVideoT
 *
 *
 *
 * @author PABLO ARIAS  <pariasm@gmail.com>
 **/

#if 0
/**
 * @brief Structure containing size informations of a video.
 *
 * @param width    : width of the image;
 * @param height   : height of the image;
 * @param channels : number of channels in the image;
 * @param frames   : number of frames in the video;
 * @param wh       : equal to width * height. Provided for convenience;
 * @param whc      : equal to width * height * channels. Provided for convenience.
 * @param whcf     : equal to width * height * frames * channels. Provided for convenience.
 * @param whf      : equal to width * height * frames. Provided for convenience.
 **/
struct VideoSize
{
	unsigned width;
	unsigned height;
	unsigned frames;
	unsigned channels;
	unsigned wh;
	unsigned whc;
	unsigned whcf;
	unsigned whf;

	inline bool operator == (const VideoSize& sz)
	{
		return (width    == sz.width     &&
		        height   == sz.height    &&
		        channels == sz.channels  &&
		        frames   == sz.frames    );
	}

	inline bool operator != (const VideoSize& sz)
	{ 
		return !operator==(sz);
	}

	//! Updates products of dimensions
	inline void update_fields(void)
	{
		wh = width * height;
		whc = wh * channels;
		whcf = whc * frames;
		whf  = wh  * frames;
	}

	//! Returns index
	inline signed index(unsigned x, unsigned y, unsigned t, unsigned c) const
	{
		assert(x < width && y < height && t < frames && c < channels);
		return t*whc + c*wh + y*width + x;
	}

	//! Returns index assuming the video has one channel
	inline unsigned index(unsigned x, unsigned y, unsigned t) const
	{
		assert(x < width && y < height && t < frames);
		return t*wh + y*width + x;
	}

	//! Compute coordinates from index
	inline
	void coords(unsigned idx, unsigned& x, unsigned& y, unsigned& t, unsigned& c) const
	{
		assert(idx < whcf);
		t = (idx      ) / whc;
		c = (idx % whc) / wh ;
		y = (idx % wh ) / width;
		x = (idx % width  );
	}

	//! Compute coordinates from index assuming the video has one channel
	inline
	void coords(unsigned idx, unsigned& x, unsigned& y, unsigned& t) const
	{
		assert(idx < whf);
		t = (idx      ) / wh;
		y = (idx % wh ) / width;
		x = (idx % width  );
	}
};

/**
 * @brief A video class template with very basic functionalities.
 *
 * @param width    : width of the video;
 * @param height   : height of the video;
 * @param channels : number of channels in the video;
 * @param frames   : number of frames in the video;
 * @param wh       : equal to width * height. Provided for convenience;
 * @param whc      : equal to width * height * channels. Provided for convenience.
 * @param whcf     : equal to width * height * frames * channels. Provided for convenience.
 * @param data     : pointer to an std::vector<float> containing the data
 **/
template <class T>
class Video
{
	public:

		//! Size
		unsigned width;
		unsigned height;
		unsigned frames;
		unsigned channels;
		unsigned wh;
		unsigned whc;
		unsigned whcf;
		unsigned whf;

		//! Data
		std::vector<T> data;

		//! Constructors
		Video(void) { }; //< empty
		Video(const Video& i_in); //< copy
		Video(const std::string i_pathToFiles,
		          unsigned i_firstFrame, unsigned i_lastFrame, unsigned i_frameStep = 1); //< from filename
		Video(unsigned i_width, unsigned i_height, unsigned i_frames, unsigned i_channels = 1);  //< alloc
		Video(unsigned i_width, unsigned i_height, unsigned i_frames, unsigned i_channels, T val);  //< init
		Video(const VideoSize& i_size);  //< alloc
		Video(const VideoSize& i_size, T val);  //< init

		//! Destructor
		~Video(void) { };

		void clear(void);
		void resize(unsigned i_width, unsigned i_height, unsigned frames, unsigned i_channels = 1);
		void resize(const VideoSize& i_size);

		//! Read/write pixel access ~ inline for efficiency
		T& operator () (unsigned idx); //< from coordinates
		T& operator () (unsigned x, unsigned y, unsigned t, unsigned c = 0); //< from coordinates

		//! Read only pixel access ~ inline for efficiency
		T operator () (unsigned idx) const; //< from index
		T operator () (unsigned x, unsigned y, unsigned t, unsigned c = 0) const; //< from coordinates

		//! Pixel access with special boundary conditions
		T& getPixelSymmetric(int x, int y, int t, unsigned c = 0);
		T  getPixelSymmetric(int x, int y, int t, unsigned c = 0) const;
		
		//! I/O
		void loadVideo(const std::string i_pathToFiles, 
		               unsigned i_firstFrame, unsigned i_lastFrame, unsigned i_frameStep = 1);
		void saveVideo(const std::string i_pathToFiles, 
		               unsigned i_firstFrame, unsigned i_frameStep = 1,
		               T i_pmin = 0, T i_pmax = 255) const;
		void saveVideoAscii(const std::string i_prefix, 
		                    unsigned i_firstFrame, unsigned i_frameStep = 1) const;

		//! Utilities
		unsigned index(unsigned x, unsigned y, unsigned t, unsigned c) const;
		void coords(unsigned index, unsigned& x, unsigned& y, unsigned& t, unsigned& c) const;

		VideoSize size(void) const
		{
			VideoSize sz = {width, height, frames, channels, wh, whc, whcf, whf};
			return sz; 
		}
};
#endif


void print_video_size(const std::string& name, const VideoSize& sz)
{
	printf("%s", name.c_str());
	printf("\tFrames: %d - Channels: %d - Height: %d - Width %d\n", 
			sz.frames, sz.channels, sz.height, sz.width);
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
	Video<float> vid1;
	print_video_size("empty video", vid1.sz);

	//! Save video / should do nothing!
	vid1.saveVideo("/tmp/vid1_empty_%02d.png",i_firstFrame, i_frameStep);

	//! Load a video through loadVideo 
	vid1.loadVideo(i_video_path, i_firstFrame, i_lastFrame, i_frameStep);
	print_video_size("loaded video1", vid1.sz);


	/*/! Accessing pixels through coordinates
	for (int f =  0; f < vid1.sz.frames; f += 2)
	for (int y = 10; y < vid1.sz.height - 10; y++)
	for (int x = 10; x < vid1.sz.width  - 10; x++)
		vid1(x,y,f,0) = 250;

	//! Accessing pixels through indices
	for (int i = 0; i < vid1.sz.whcf; i += 17) vid1(i) = 0;//*/

	/*/! Save video
	vid1.saveVideo("/tmp/vid1_%02d.png", i_firstFrame, i_frameStep);
//	vid1.saveVideoAscii("/tmp/vid1"    , i_firstFrame, i_frameStep);//*/

	/*/! Subdivide into smaller videos with border
	unsigned nvids = 8;
	unsigned border = 10;
	std::vector<Video<float> > subvids1;
	VideoUtils::subDivide(vid1, subvids1, border, nvids); 

	for (int n = 0; n < subvids1.size(); n++)
	{
		char name[1024];
		sprintf(name, "/tmp/vid1_sub%02d_%%02d.png", n);
		subvids1[n].saveVideo(name, i_firstFrame, i_frameStep);
	}

	//! Join again into large video, removing border
	Video<float> vid3(vid1.sz);
	VideoUtils::subBuild (subvids1, vid3, border); 

	float psnr_13, rmse_13;
	VideoUtils::computePSNR(vid1, vid3, psnr_13, rmse_13);
	printf("After joining the subvideos... RMSE: %f - PSNR: %f\n", rmse_13, psnr_13);//*/

	//! Subdivide into smaller videos with border
	int nvids = 8;
	int border = 10;
	std::vector<Video<float> > subvids1;
	std::vector<VideoUtils::CropPosition > crops; 
	VideoUtils::subDivideTight(vid1, subvids1, crops, border, nvids); 

	for (int n = 0; n < subvids1.size(); n++)
	{
		char name[1024];
		sprintf(name, "/tmp/vid1_sub%02d_%%02d.png", n);
		subvids1[n].saveVideo(name, i_firstFrame, i_frameStep);
	}

	//! Join again into large video, removing border
	Video<float> vid3(vid1.sz);
	VideoUtils::subBuildTight (subvids1, vid3, border); 

	float psnr_13, rmse_13;
	VideoUtils::computePSNR(vid1, vid3, psnr_13, rmse_13);
	printf("After joining the subvideos... RMSE: %f - PSNR: %f\n", rmse_13, psnr_13);//*/

	{
		vid3.saveVideo("/tmp/vid1_built_%02d.png", i_firstFrame, i_frameStep);
	}


	/*/! Pad video by symmetry
	Video<float> vid1_sym;
	VideoUtils::addBorder(vid1, vid1_sym, 4, true);
	print_video_size("video with border added by symmetrizing", vid1_sym.sz);
	vid1_sym.saveVideo("/tmp/vid1_sym_%02d.png", i_firstFrame, i_frameStep);//*/

	/*/! Change colorspace
	VideoUtils::transformColorSpace(vid1, true);
	vid1.saveVideoAscii("/tmp/vid1_yuv", i_firstFrame, i_frameStep);

	VideoUtils::transformColorSpace(vid1, false);//*/

	//! Add noise
	Video<float> vid1_noise;
	VideoUtils::addNoise(vid1, vid1_noise, i_sigma, true);

	//! Save video
	vid1_noise.saveVideo("/tmp/vid1_noise_%02d.png",i_firstFrame, i_frameStep);

	//! Compute PSNR
	float psnr, rmse;
	VideoUtils::computePSNR(vid1, vid1_noise, psnr, rmse);
	printf("Computed RMSE: %f - PSNR: %f\n", rmse, psnr);//*/

	/*/! Compute difference
	Video<float> diff;
	VideoUtils::computeDiff(vid1, vid1_noise, diff, i_sigma);
	print_video_size("difference", diff.sz);

	//! Save difference video
	diff.saveVideo("/tmp/diff_%02d.png",i_firstFrame, i_frameStep);//*/

	/*/! Load a video through constructor
	Video<float> vid2(i_video_path, i_firstFrame, i_lastFrame, i_frameStep);
	print_video_size("loaded video2", vid2.sz);

	//! Save video
	vid2.saveVideo("/tmp/vid2_%02d.png",i_firstFrame, i_frameStep);

	//! Create an empty copy
	Video<float> vid3(vid1);
	print_video_size("copied video3", vid3.sz);

	//! Save
	vid2.saveVideo("/tmp/vid3_%02d.png",i_firstFrame, i_frameStep);//*/

	return EXIT_SUCCESS;
}
