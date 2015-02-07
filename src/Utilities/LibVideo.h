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
#include <climits>

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
 * @brief A float video class with very basic functionalities.
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
class Video_f32
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
		std::vector<float> data;

		//! Constructors
		Video_f32(void); //< empty
		Video_f32(const Video_f32& i_in); //< copy
		Video_f32(const std::string i_pathToFiles,
		          unsigned i_firstFrame, unsigned i_lastFrame, unsigned i_frameStep = 1); //< from filename
		Video_f32(unsigned i_width, unsigned i_height, unsigned i_frames, unsigned i_channels = 1);  //< alloc
		Video_f32(unsigned i_width, unsigned i_height, unsigned i_frames, unsigned i_channels, float val);  //< init
		Video_f32(const VideoSize& i_size);  //< alloc
		Video_f32(const VideoSize& i_size, float val);  //< init

		//! Destructor
		~Video_f32(void) { };

		void clear(void);
		void resize(unsigned i_width, unsigned i_height, unsigned frames, unsigned i_channels = 1);
		void resize(const VideoSize& i_size);

		//! Read/write pixel access ~ inline for efficiency
		float& operator () (unsigned idx); //< from coordinates
		float& operator () (unsigned x, unsigned y, unsigned t, unsigned c = 0); //< from coordinates

		//! Read only pixel access ~ inline for efficiency
		float operator () (unsigned idx) const; //< from index
		float operator () (unsigned x, unsigned y, unsigned t, unsigned c = 0) const; //< from coordinates

		//! Pixel access with special boundary conditions
		float& getPixelSymmetric(int x, int y, int t, unsigned c = 0);
		float  getPixelSymmetric(int x, int y, int t, unsigned c = 0) const;
		
		//! I/O
		int loadVideo(const std::string i_pathToFiles, 
		              unsigned i_firstFrame, unsigned i_lastFrame, unsigned i_frameStep = 1);
		int saveVideo(const std::string i_pathToFiles, 
		              unsigned i_firstFrame, unsigned i_frameStep = 1,
		              float i_pmin = 0, float i_pmax = 255) const;
		int saveVideoAscii(const std::string i_prefix, 
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

inline float& Video_f32::getPixelSymmetric(
	int x
,	int y
,	int t
,	unsigned c
) {
	// NOTE: assumes that -width+1 < x < 2*(width -1)
	assert(-(int)width   < x && x < 2*(int)width -1&&
	       -(int)height  < y && y < 2*(int)height-1&&
	       -(int)frames  < t && t < 2*(int)frames-1);
	// symmetrize
	x = (x < 0) ? -x : (x >= (int)width  ) ? 2*(int)width  - 2 - x : x ;
	y = (y < 0) ? -y : (y >= (int)height ) ? 2*(int)height - 2 - y : y ;
	t = (t < 0) ? -t : (t >= (int)frames ) ? 2*(int)frames - 2 - t : t ;

	return data[index(x,y,t,c)];
}

inline float Video_f32::getPixelSymmetric(
	int x
,	int y
,	int t
,	unsigned c
) const {
	// NOTE: assumes that -width+1 < x < 2*(width -1)
	assert(-(int)width   < x && x < 2*(int)width  - 1 &&
	       -(int)height  < y && y < 2*(int)height - 1 &&
	       -(int)frames  < t && t < 2*(int)frames - 1 );
	// symmetrize
	x = (x < 0) ? -x : (x >= (int)width  ) ? 2*(int)width  - 2 - x : x ;
	y = (y < 0) ? -y : (y >= (int)height ) ? 2*(int)height - 2 - y : y ;
	t = (t < 0) ? -t : (t >= (int)frames ) ? 2*(int)frames - 2 - t : t ;

	return data[index(x,y,t,c)];
}

inline unsigned Video_f32::index(
	unsigned x
,	unsigned y
,	unsigned t
,	unsigned c
) const {
	assert(x < width && y < height && t < frames && c < channels);
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
	int computePSNR(
	    Video_f32 const& i_vid1
	,   Video_f32 const& i_vid2
	,   float &o_psnr
	,   float &o_rmse
//	,   const char* p_videoName
//	,   const bool p_verbose = false
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
	);
	
	/**
	 * @brief Add boundaries by symmetry
	 *
	 * @param io_vid : original image;
	 * @param io_vidSym : contain io_im symmetrized;
	 *
	 * @return none.
	 **/
	void addBorder(
		Video_f32 const& i_vid1
	,	Video_f32 &o_vid2
	,	const unsigned p_borderSize
	,	const bool p_isForward
	);

	/**
	 * @brief 'Generalized' croping of a video (cropped video may be larger than original).
	 *
	 * @param i_vid1 : original video;
	 * @param o_vid2 : output video, already allocated to desired size;
	 * @param i_origin2 : vid1 coordinates of vid2 origin. Origin coordinates
	 * larger than corresponding vid1 dimension are redefined to center the crop
	 * in that dimension.
	 *
	 * @return none.
	 **/
	void crop(
		Video_f32 const &i_vid1
	,	Video_f32 &o_vid2
	,	int p_origin_t = INT_MAX
	,	int p_origin_x = INT_MAX
	,	int p_origin_y = INT_MAX
	);

	/**
	 * @brief 'Generalized' croping of a video (cropped video may be larger than original).
	 *
	 * @param i_vid1 : original video;
	 * @param o_vid2 : output video, already allocated to desired size;
	 * @param i_origin2 : vid1 coordinates of vid2 origin. Origin coordinates
	 * larger than corresponding vid1 dimension are redefined to center the crop
	 * in that dimension.
	 *
	 * @return none.
	 **/
	void crop(
		Video_f32 const &i_vid1
	,	Video_f32 &o_vid2
	,	const int * const p_origin
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
		Video_f32 &io_vid
	,	const bool p_isForward
	);
	
	/**
	 * @brief Subdivide a video into small sub-videos
	 *
	 * @param i_video : image to subdivide;
	 * @param o_videoSub : will contain all sub-videos;
	 * @param p_N : boundary around sub-videos;
	 * @param p_nb : number of sub-videos wanted. Need to be a power of 2.
	 *
	 * @return EXIT_FAILURE in case of problems.
	 **/
	int subDivide(
		Video_f32 const& i_vid
	,	std::vector<Video_f32> &o_vidSub
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
	 	std::vector<Video_f32> const& i_vidSub
	,	Video_f32 &o_vid
	,	const unsigned p_N
	);
}

#endif // LIB_VIDEO_H_INCLUDED
