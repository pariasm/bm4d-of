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
#ifndef LIB_VIDEO_H_INCLUDED
#define LIB_VIDEO_H_INCLUDED

#include <vector>
#include <string>
#include <fftw3.h>
#include <cassert>

/**
 * @brief Structure containing size informations of a video.
 *
 * @param width     : width of the image;
 * @param height    : height of the image;
 * @param nChannels : number of channels in the image;
 * @param nFrames   : number of frames in the video;
 * @param wh        : equal to width * height. Provided for convenience;
 * @param whc       : equal to width * height * nChannels. Provided for convenience.
 * @param whcf      : equal to width * height * nFrames * nChannels. Provided for convenience.
 **/
struct VideoSize
{
	unsigned width;
	unsigned height;
	unsigned nChannels;
	unsigned nFrames;
	unsigned wh;
	unsigned whc;
	unsigned whcf;
};

/**
 * @brief A float video class with very basic functionalities.
 *
 * @param width     : width of the video;
 * @param height    : height of the video;
 * @param nChannels : number of channels in the video;
 * @param nFrames   : number of frames in the video;
 * @param wh        : equal to width * height. Provided for convenience;
 * @param whc       : equal to width * height * nChannels. Provided for convenience.
 * @param whcf      : equal to width * height * nFrames * nChannels. Provided for convenience.
 * @param data      : pointer to an std::vector<float> containing the data
 **/
class Video_f32
{
	public:

		//! Size
		unsigned width;
		unsigned height;
		unsigned nFrames;
		unsigned nChannels;
		unsigned wh;
		unsigned whc;
		unsigned whcf;

		//! Data
		std::vector<float> data;

		//! Constructors
		Video_f32(void); //< empty
		Video_f32(const Video_f32& i_in); //< copy
		Video_f32(const std::string i_pathToFiles,
		          unsigned i_firstFrame, unsigned i_lastFrame, unsigned i_frameStep = 1); //< from filename
		Video_f32(unsigned i_width, unsigned i_height, unsigned nFrames, unsigned i_nChannels = 1);  //< alloc

		//! Destructor
		~Video_f32(void) { };

		void clear(void);

		//! Read/write pixel access ~ inplace for efficiency
		float& operator () (unsigned idx); //< from coordinates
		float& operator () (unsigned x, unsigned y, unsigned t, unsigned c = 0); //< from coordinates

		//! Read only pixel access ~ inplace for efficiency
		float operator () (unsigned idx) const; //< from index
		float operator () (unsigned x, unsigned y, unsigned t, unsigned c = 0) const; //< from coordinates
		
		//! I/O
		int loadVideo(const std::string i_pathToFiles, 
		              unsigned i_firstFrame, unsigned i_lastFrame, unsigned i_frameStep = 1);
		int saveVideo(const std::string i_pathToFiles, 
		              unsigned i_firstFrame, unsigned i_frameStep = 1,
		              float i_pmin = 0, float i_pmax = 255) const;

		//! Utilities
		unsigned index(unsigned x, unsigned y, unsigned t, unsigned c = 0) const;
		void coords(unsigned index, unsigned& x, unsigned& y, unsigned& t, unsigned& c) const;
};

inline float& Video_f32::operator () (unsigned idx) 
{
	assert(idx < whcf);
	return data[idx];
}

inline float& Video_f32::operator() (
	unsigned x
,	unsigned y
,	unsigned t
,	unsigned c
){
	return data[index(x,y,t,c)];
}

inline float Video_f32::operator () (unsigned idx) const
{
	assert(idx < whcf);
	return data[idx];
}

inline float Video_f32::operator() (
	unsigned x
,	unsigned y
,	unsigned t
,	unsigned c
) const {
	return data[index(x,y,t,c)];
}

inline unsigned Video_f32::index(
	unsigned x
,	unsigned y
,	unsigned t
,	unsigned c
) const {
	assert(x < width && y < height && t < nFrames && c < nChannels);
	return t*whc + c*wh + y*width + x;
}

inline void Video_f32::coords(
	unsigned idx
,	unsigned& x
,	unsigned& y
,	unsigned& t
,	unsigned& c
) const {
	assert(idx < whcf);
	t = (idx      ) / whc;
	c = (idx % whc) / wh ;
	y = (idx % wh ) / width;
	x = (idx % width  );
//	c = (idx - t*whc)  / wh;
//	y = (idx - t*whc - c*wh)  / width;
//	x = (idx - t*whc - c*wh - y*width) ;
}

#if 0
//! Utilities for video
namespace VideoUtils
{

/**
 * @brief add noise to video.
 *
 * @param i_vid : original noise-free image;
 * @param o_vidNoisy = vid + noise;
 * @param p_sigma : standard deviation of the noise.
 *
 * @return none.
 **/
void addNoise(
    Video_f32 const& i_vid
,   Video_f32 &o_vidNoisy
,   const float p_sigma
,   const bool p_verbose = false
);

/**
 * @brief Compute PSNR and RMSE between i_vid1 and i_vid2
 *
 * @param i_vid1 : video 1;
 * @param i_vid2 : video 2;
 * @param o_psnr  : will contain the PSNR;
 * @param o_rmse  : will contain the RMSE;
 * @param p_videoName: name of the video;
 * @param p_verbose: if true, print values of PSNR and RMSE.
 *
 * @return EXIT_FAILURE if both videos haven't the same size.
 **/
int computePsnr(
    Video_f32 const& i_vid1
,   Video_f32 const& i_vid2
,   float &o_psnr
,   float &o_rmse
,   const char* p_videoName
,   const bool p_verbose = false
);

/**
 * @brief Compute a difference image between i_vid1 and i_vid2.
 *
 * @param i_vid1: reference image;
 * @param i_vid2: image to compare;
 * @param o_vidDiff: will contain the difference;
 * @param p_sigma: standard deviation of the noise;
 * @param p_min, p_max : range of data (usually [0, 255]);
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if i_im1 and i_im2 don't have the same size.
 **/
int computeDiff(
    Video_f32 const& i_vid1
,   Video_f32 const& i_vid2
,   Video_f32 &o_vidDiff
,   const float p_sigma
,   const float p_min = 0.f
,   const float p_max = 255.f
,   const bool p_verbose = false
);

/**
 * @brief Add boundary by symmetry.
 *
 * @param i_vid : video to symmetrize;
 * @param o_vidSym : will contain i_vid with symmetrized boundaries;
 *
 * @return none.
 **/
int addBoundary(
	std::vector<float> const& i_vid
,	std::vector<float> &o_vidSym
);

/**
 * @brief Remove boundaries added with addBoundary
 *
 * @param o_vid : will contain the inner image;
 * @param i_vidSym : contains i_vid with symmetrized boundaries;
 *
 * @return none.
 **/
int removeBoundary(
	std::vector<float> &o_vid
,	std::vector<float> const& i_vidSym
);

/**
 * @brief Add boundaries by symmetry
 *
 * @param io_vid : original image;
 * @param io_vidSym : contain io_im symmetrized;
 *
 * @return none.
 **/
void symmetrizeImage(
	std::vector<float> const& i_vid1
,	std::vector<float> &o_vid2
,	const unsigned p_borderSize
,	const bool p_isForward
);

/**
 * @brief Transform the color space of an video, from RGB to YUV, or vice-versa.
 *
 * @param io_vid: image on which the transform will be applied;
 * @param p_isForward: if true, go from RGB to YUV, otherwise go from YUV to RGB.
 *
 * @return none.
 **/
void transformColorSpace(
	std::vector<float> &io_vid
,	const bool p_isForward
);

/**
 * @brief Subdivide a video into small sub-videos
 *
 * @param i_video : image to subdivide;
 * @param o_videoSub : will contain all sub-videos;
 * @param p_videoSizeSub : size of sub-videos;
 * @param p_N : boundary around sub-videos;
 * @param p_nb : number of sub-videos wanted. Need to be a power of 2.
 *
 * @return EXIT_FAILURE in case of problems.
 **/
int subDivide(
	Video_f32 const& i_im
,	std::vector<Video_f32> &o_imSub
,	VideoSize &p_imSizeSub
,	const unsigned p_N
,	const unsigned p_nb
);

/**
 * @brief Reconstruct an video from its small sub-videos
 *
 * @param o_vid : image to reconstruct;
 * @param i_vidSub : will contain all sub-images;
 * @param p_N : boundary around sub-videos.
 *
 * @return EXIT_FAILURE in case of problems.
 **/
int subBuild(
	Video_f32 &o_vid
,	std::vector<Video_f32> const& i_vidSub
,	const unsigned p_N
);
#endif

#endif // LIB_VIDEO_H_INCLUDED
