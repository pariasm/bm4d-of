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
#ifndef LIB_VIDEOT_HPP_INCLUDED
#define LIB_VIDEOT_HPP_INCLUDED

#include <vector>
#include <string>
#include <fftw3.h>
#include <cassert>
#include <climits>
#include <stdexcept>

#include <math.h>
#include "mt19937ar.h"

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

	//! Constuctros
	VideoSize(void)
		: width(0), height(0), frames(0), channels(0)
	{
		update_fields();
	}

	VideoSize(unsigned w, unsigned h, unsigned f, unsigned c)
		: width(w), height(h), frames(f), channels(c)
	{
		update_fields();
	}

	//! Comparison operators
	inline bool operator == (const VideoSize& sz) const
	{
		return (width    == sz.width     &&
		        height   == sz.height    &&
		        channels == sz.channels  &&
		        frames   == sz.frames    );
	}

	inline bool operator != (const VideoSize& sz) const
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
	inline unsigned index(unsigned x, unsigned y, unsigned t, unsigned c) const
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
 * @param sz       : VideoSize structure with size of the video;
 * @param data     : pointer to an std::vector<T> containing the data
 **/
template <class T>
class Video
{
	public:

		//! Size
		VideoSize sz;

		//! Data
		std::vector<T> data;

		//! Constructors
		Video(void); //< empty
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
};

// Implementations

template <class T> 
Video<T>::Video(void)
	: sz(), data(0)
{ }
	
template <class T> 
Video<T>::Video(const Video& i_in)
	: sz(i_in.sz), data(i_in.data)
{ }

template <class T> 
Video<T>::Video(
	const std::string i_pathToFiles
,	unsigned i_firstFrame
,	unsigned i_lastFrame
,	unsigned i_frameStep
) : sz(), data(0)
{
	loadVideo(i_pathToFiles, i_firstFrame, i_lastFrame, i_frameStep);
}
	
template <class T> 
Video<T>::Video(
	unsigned i_width
,	unsigned i_height
,	unsigned i_frames
,	unsigned i_channels
)
	: sz(i_width, i_height, i_frames, i_channels)
	, data(sz.whcf)
{ }

template <class T> 
Video<T>::Video(
	unsigned i_width
,	unsigned i_height
,	unsigned i_frames
,	unsigned i_channels
,	T val
)
	: sz(i_width, i_height, i_frames, i_channels)
	, data(sz.whcf, val)
{ }
	
template <class T> 
Video<T>::Video(const VideoSize& i_size, T val)
	: sz(i_size)
	, data(sz.whcf, val)
{ }
	
template <class T> 
Video<T>::Video(const VideoSize& i_size)
	: sz(i_size)
	, data(sz.whcf)
{ }
	
template <class T> 
void Video<T>::clear(void)
{
	sz.width = 0;
	sz.height = 0;
	sz.frames = 0;
	sz.channels = 0;
	sz.update_fields();
	data.clear();
}

template <class T> 
void Video<T>::resize(const VideoSize& i_size)
{
	if (sz != i_size)
	{
		clear();
		sz = i_size;
		data.resize(sz.whcf);
	}
}

template <class T> 
void Video<T>::resize(
	unsigned i_width
,	unsigned i_height
,	unsigned i_frames
,	unsigned i_channels
){
	resize(VideoSize(i_width, i_height, i_frames, i_channels));
}
	
template <class T> 
void Video<T>::loadVideo(
    const std::string i_pathToFiles
,   unsigned i_firstFrame
,   unsigned i_lastFrame
,   unsigned i_frameStep
){
	throw std::runtime_error("Video<T>::loadVideo(...) is only implemented "
			"for T = float");
}

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
	
template <class T> 
void Video<T>::saveVideo(
	const std::string i_pathToFiles
,	unsigned i_firstFrame
,	unsigned i_frameStep
,	T i_pmin
,	T i_pmax
) const {
	throw std::runtime_error("Video<T>::saveVideo(...) is only implemented "
			"for T = float");
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

template <class T> 
void Video<T>::saveVideoAscii(
	const std::string i_prefix
,	unsigned i_firstFrame
,	unsigned i_frameStep
) const {
	throw std::runtime_error("Video<T>::saveVideoAscii(...) is only implemented "
			"for T = float");
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

template <class T>
inline T& Video<T>::operator () (unsigned idx) 
{
	assert(idx < sz.whcf);
	return data[idx];
}

template <class T>
inline T& Video<T>::operator() (
	unsigned x
,	unsigned y
,	unsigned t
,	unsigned c
){
	return data[sz.index(x,y,t,c)];
}

template <class T>
inline T Video<T>::operator () (unsigned idx) const
{
	assert(idx < whcf);
	return data[idx];
}

template <class T>
inline T Video<T>::operator() (
	unsigned x
,	unsigned y
,	unsigned t
,	unsigned c
) const {
	return data[sz.index(x,y,t,c)];
}

template <class T>
inline T& Video<T>::getPixelSymmetric(
	int x
,	int y
,	int t
,	unsigned c
) {
	// NOTE: assumes that -width+1 < x < 2*(width -1)
	assert(-(int)sz.width   < x && x < 2*(int)sz.width -1&&
	       -(int)sz.height  < y && y < 2*(int)sz.height-1&&
	       -(int)sz.frames  < t && t < 2*(int)sz.frames-1);
	// symmetrize
	x = (x < 0) ? -x : (x >= (int)sz.width  ) ? 2*(int)sz.width  - 2 - x : x ;
	y = (y < 0) ? -y : (y >= (int)sz.height ) ? 2*(int)sz.height - 2 - y : y ;
	t = (t < 0) ? -t : (t >= (int)sz.frames ) ? 2*(int)sz.frames - 2 - t : t ;

	return data[sz.index(x,y,t,c)];
}

template <class T>
inline T Video<T>::getPixelSymmetric(
	int x
,	int y
,	int t
,	unsigned c
) const {
	// NOTE: assumes that -width+1 < x < 2*(width -1)
	assert(-(int)sz.width   < x && x < 2*(int)sz.width  - 1 &&
	       -(int)sz.height  < y && y < 2*(int)sz.height - 1 &&
	       -(int)sz.frames  < t && t < 2*(int)sz.frames - 1 );
	// symmetrize
	x = (x < 0) ? -x : (x >= (int)sz.width  ) ? 2*(int)sz.width  - 2 - x : x ;
	y = (y < 0) ? -y : (y >= (int)sz.height ) ? 2*(int)sz.height - 2 - y : y ;
	t = (t < 0) ? -t : (t >= (int)sz.frames ) ? 2*(int)sz.frames - 2 - t : t ;

	return data[sz.index(x,y,t,c)];
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
	template <class T>
	void addNoise(
	    Video<T> const& i_vid
	,   Video<T> &o_vidNoisy
	,   const float p_sigma
	,   const bool p_verbose = false
	){
		if (p_verbose) printf("Add noise [sigma = %f] ... ", p_sigma);

		//! Initialization
		o_vidNoisy = i_vid;
//		mt_init_genrand((unsigned long int) time (NULL) +
//		                (unsigned long int) getpid());
		mt_init_genrand(0); printf("Warning: random generator seed is 0 ");

		//! Add noise
		for (unsigned k = 0; k < i_vid.sz.whcf; k++)
		{
			const double a = mt_genrand_res53();
			const double b = mt_genrand_res53();
			o_vidNoisy(k) += (T) (p_sigma *
				sqrtl(-2.0l * log(a)) * cos(2.0l * M_PI * b));
		}

		if (p_verbose) printf("done.\n");
	}
	
	/**
	 * @brief Compute PSNR and RMSE between i_vid1 and i_vid2
	 *
	 * @param i_vid1 : video 1;
	 * @param i_vid2 : video 2;
	 * @param o_psnr  : will contain the PSNR;
	 * @param o_rmse  : will contain the RMSE;
	 *
	 * @return none.
	 **/
	template <class T>
	void computePSNR(
	    Video<T> const& i_vid1
	,   Video<T> const& i_vid2
	,   float &o_psnr
	,   float &o_rmse
	){
		if (i_vid1.sz != i_vid2.sz)
			throw std::runtime_error("VideoUtils::computePSNR: videos have different sizes");

		float sum = 0.f;
		for (unsigned k = 0; k < i_vid1.sz.whcf; k++)
			sum += ((float)i_vid1(k) - (float)i_vid2(k)) *
			       ((float)i_vid1(k) - (float)i_vid2(k));

		o_rmse = sqrtf(sum / (float) i_vid1.sz.whcf);
		o_psnr = 20.f * log10f(255.f / o_rmse);

		return;
	}
	
	/**
	 * @brief Compute a difference image between i_vid1 and i_vid2.
	 *
	 * @param i_vid1: reference image;
	 * @param i_vid2: image to compare;
	 * @param o_vidDiff: will contain the difference;
	 * @param p_sigma: standard deviation of the noise;
	 * @param p_min, p_max : range of data (usually [0, 255]);
	 *
	 * @return none.
	 **/
	template <class T>
	void computeDiff(
	    Video<T> const& i_vid1
	,   Video<T> const& i_vid2
	,   Video<float> &o_vidDiff
	,   const float p_sigma
	,   const T p_min = 0.f
	,   const T p_max = 255.f
	){
		if (i_vid1.sz != i_vid2.sz)
			throw std::runtime_error("VideoUtils::computeDiff: videos have different sizes");

		o_vidDiff.resize(i_vid1.sz);
		for (unsigned k = 0; k < i_vid1.sz.whcf; k++)
		{
			float value =  ((float)i_vid1(k) - (float)i_vid2(k) + p_sigma) * p_max / (2.f * p_sigma);
			o_vidDiff(k) = clip(value, p_min, p_max);
		}

		return;
	}

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
	template <class T>
	void crop(
		Video<T> const &i_vid1
	,	Video<T> &o_vid2
	,	int p_origin_t = INT_MAX
	,	int p_origin_x = INT_MAX
	,	int p_origin_y = INT_MAX
	){
		assert(o_vid2.sz.channels == i_vid1.sz.channels);

		//! Redefine invalid origin coordinates to default (centered crop)
		if (p_origin_t > (int)i_vid1.sz.frames) p_origin_t = ((int)i_vid1.sz.frames - (int)o_vid2.sz.frames)/2;
		if (p_origin_x > (int)i_vid1.sz.width ) p_origin_x = ((int)i_vid1.sz.width  - (int)o_vid2.sz.width )/2;
		if (p_origin_y > (int)i_vid1.sz.height) p_origin_y = ((int)i_vid1.sz.height - (int)o_vid2.sz.height)/2;

		// TODO: more efficient implementation
		for (int      f = 0; f < o_vid2.sz.frames  ; f++)
		for (unsigned c = 0; c < o_vid2.sz.channels; c++)
		for (int      y = 0; y < o_vid2.sz.height  ; y++)
		for (int      x = 0; x < o_vid2.sz.width   ; x++)
			o_vid2(x,y,f,c) = 
				i_vid1.getPixelSymmetric(x + p_origin_x, y + p_origin_y, f + p_origin_t, c);
	}

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
	template <class T>
	void crop(
		Video<T> const &i_vid1
	,	Video<T> &o_vid2
	,	const int * const p_origin
	){
		crop(i_vid1, o_vid2, p_origin[2], p_origin[0], p_origin[1]);
	}
	
	
	/**
	 * @brief Add boundaries by symmetry
	 *
	 * @param io_vid : original image;
	 * @param io_vidSym : contain io_im symmetrized;
	 *
	 * @return none.
	 **/
	template <class T>
	void addBorder(
		Video<T> const& i_vid1
	,	Video<T> &o_vid2
	,	const unsigned p_borderSize
	,	const bool p_isForward
	){
		//! Sizes
		const unsigned w2 = i_vid1.sz.width  + (p_isForward ? 2*p_borderSize : - 2*p_borderSize);
		const unsigned h2 = i_vid1.sz.height + (p_isForward ? 2*p_borderSize : - 2*p_borderSize);
		const unsigned f2 = i_vid1.sz.frames + (p_isForward ? 2*p_borderSize : - 2*p_borderSize);
		const unsigned ch = i_vid1.sz.channels;

		//! Position of vid2 origin in vid1 coordinates
		const int tx = p_isForward ? -p_borderSize : p_borderSize;
		const int ty = p_isForward ? -p_borderSize : p_borderSize;
		const int tf = p_isForward ? -p_borderSize : p_borderSize;

		//! Resize output image, if necessary
		o_vid2.resize(w2, h2, f2, ch);

		//! Call generalized crop
		crop(i_vid1, o_vid2, tx, ty, tf);
	}

	/**
	 * @brief Transform the color space of an video, from RGB to YUV, or vice-versa.
	 *
	 * @param io_vid: image on which the transform will be applied;
	 * @param p_isForward: if true, go from RGB to YUV, otherwise go from YUV to RGB.
	 *
	 * @return none.
	 **/
	template <class T>
	void transformColorSpace(
		Video<T> &io_vid
	,	const bool p_isForward
	){
		//! If the image as only one channel, do nothing
		if (io_vid.sz.channels == 1) return;

		//! Initialization
		const unsigned width  = io_vid.sz.width;
		const unsigned height = io_vid.sz.height;
		const unsigned chnls  = io_vid.sz.channels;
		const unsigned wh     = io_vid.sz.wh;

		for (int f = 0; f < io_vid.sz.frames; f++)
			if (p_isForward) //< RGB to YUV
			{
				if (chnls == 3)
				{
					const unsigned red   = f * io_vid.sz.whc;
					const unsigned green = red + wh;
					const unsigned blue  = green + wh;
					const float a = 1.f / sqrtf(3.f);
					const float b = 1.f / sqrtf(2.f);
					const float c = 2.f * a * sqrtf(2.f);

					float yuv[3];
					for (unsigned k = 0; k < wh; k++)
					{
						yuv[0] = a * ((float)io_vid(k + red) + (float)io_vid(k + green) + (float)io_vid(k + blue));
						yuv[1] = b * ((float)io_vid(k + red) - (float)io_vid(k + blue));
						yuv[2] = c * (0.25f * (float)io_vid(k + red ) - 0.5f * (float)io_vid(k + green)
								      + 0.25f * (float)io_vid(k + blue));

						io_vid(k + red  ) = (T)yuv[0];
						io_vid(k + green) = (T)yuv[1];
						io_vid(k + blue ) = (T)yuv[2];
					}
				}
				else //< chnls == 4
				{
					const unsigned Gr = f * io_vid.sz.whc;
					const unsigned R  = Gr + wh;
					const unsigned B  = R  + wh;
					const unsigned Gb = B  + wh;
					const float a = 0.5f;
					const float b = 1.f / sqrtf(2.f);

					float wtf[4];
					for (unsigned k = 0; k < wh; k++)
					{
						wtf[0] = a * ( (float)io_vid(k + Gr) + (float)io_vid(k + R ) +
						               (float)io_vid(k + B ) + (float)io_vid(k + Gb));
						wtf[1] = b * ( (float)io_vid(k + R ) - (float)io_vid(k + B ));
						wtf[2] = a * (-(float)io_vid(k + Gr) + (float)io_vid(k + R ) +
						               (float)io_vid(k + B ) - (float)io_vid(k + Gb));
						wtf[3] = b * (-(float)io_vid(k + Gr) + (float)io_vid(k + Gb));
						io_vid(k + Gr) = (T)wtf[0];
						io_vid(k + R ) = (T)wtf[1];
						io_vid(k + B ) = (T)wtf[2];
						io_vid(k + Gb) = (T)wtf[3];
					}
				}
			}
			else //< YUV to RGB
			{
				if (chnls == 3)
				{
					const unsigned red   = f * io_vid.sz.whc;
					const unsigned green = red + wh;
					const unsigned blue  = green + wh;
					const float a = 1.f / sqrtf(3.f);
					const float b = 1.f / sqrtf(2.f);
					const float c = a / b;

					float rgb[3];
					for (unsigned k = 0; k < wh; k++)
					{
						rgb[0] = a * (float)io_vid(k + red) + b * (float)io_vid(k + green)
						                             + c * 0.5f * (float)io_vid(k + blue );
						rgb[1] = a * (float)io_vid(k + red) - c * (float)io_vid(k + blue );
						rgb[2] = a * (float)io_vid(k + red) - b * (float)io_vid(k + green)
						                      + c * 0.5f * (float)io_vid(k + blue );
						io_vid(k + red  ) = (T)rgb[0];
						io_vid(k + green) = (T)rgb[1];
						io_vid(k + blue ) = (T)rgb[2];
					}
				}
				else //! chnls == 4
				{
					const unsigned Gr = f * io_vid.sz.whc;
					const unsigned R  = Gr + wh;
					const unsigned B  = R  + wh;
					const unsigned Gb = B  + wh;
					const float a = 0.5f;
					const float b = 1.f / sqrtf(2.f);

					float wtf[4];
					for (unsigned k = 0; k < wh; k++)
					{
						wtf[0] = a * (float)io_vid(k + Gr) - a * (float)io_vid(k + B) - b * (float)io_vid(k + Gb);
						wtf[1] = a * (float)io_vid(k + Gr) + b * (float)io_vid(k + R) + a * (float)io_vid(k + B );
						wtf[2] = a * (float)io_vid(k + Gr) - b * (float)io_vid(k + R) + a * (float)io_vid(k + B );
						wtf[3] = a * (float)io_vid(k + Gr) - a * (float)io_vid(k + B) + b * (float)io_vid(k + Gb);
						io_vid(k + Gr) = (T)wtf[0];
						io_vid(k + R ) = (T)wtf[1];
						io_vid(k + B ) = (T)wtf[2];
						io_vid(k + Gb) = (T)wtf[3];
					}
				}
			}
	}
	
	/**
	 * @brief Subdivide a video into small sub-videos
	 *
	 * @param i_video : image to subdivide;
	 * @param o_videoSub : will contain all sub-videos;
	 * @param p_N : boundary around sub-videos;
	 * @param p_nb : number of sub-videos wanted. Need to be a power of 2.
	 *
	 * @return none.
	 **/
	template <class T>
	void subDivide(
		Video<T> const& i_vid
	,	std::vector<Video<T> > &o_vidSub
	,	const unsigned p_N
	,	const unsigned p_nb
	){
		/* FIXME current version splits the video only spatially. 
		 *       The reason is to mantain consistency with Marc's
		 *       code. For its proper extension to video, we need
		 *       to determine how to split the video in space and
		 *       time. */
		
		//! Determine number of sub-images
		unsigned nW, nH;
		determineFactor(p_nb, nW, nH);
		const unsigned wTmp = ceil(float(i_vid.sz.width ) / float(nW)); // sizes w/out 
		const unsigned hTmp = ceil(float(i_vid.sz.height) / float(nH)); //     borders

		//! Obtain sub-images
		VideoSize imSubSize;
		imSubSize.width    = wTmp + 2 * p_N; // each sub-image has border
		imSubSize.height   = hTmp + 2 * p_N;
		imSubSize.frames   = i_vid.sz.frames; // NOTE: same frames as original
		imSubSize.channels = i_vid.sz.channels;
		imSubSize.update_fields();

		o_vidSub.resize(p_nb);
		for (unsigned p = 0, n = 0; p < nH; p++)
		for (unsigned q = 0;        q < nW; q++, n++)
		{
			o_vidSub[n].resize(imSubSize);

			// The origin is shifted -p_N to account for the subimage border
			int origin[3] = {q * wTmp - p_N, p * hTmp - p_N, 0};

			// Crop using symmetric boundary conditions
			VideoUtils::crop(i_vid, o_vidSub[n], origin);
		}

		return;
	}
	
	/**
	 * @brief Reconstruct an video from its small sub-videos
	 *
	 * @param o_vid : image to reconstruct;
	 * @param i_vidSub : will contain all sub-images;
	 * @param p_N : boundary around sub-videos.
	 *
	 * @return none.
	 **/
	template <class T>
	void subBuild(
	 	std::vector<Video<T> > const& i_vidSub
	,	Video<T> &o_vid
	,	const unsigned p_N
	){
		/* FIXME current version builds a video that has been split
		 *       only spatially by subDivide. 
		 *       The reason is to mantain consistency with Marc's
		 *       code. For its proper extension to video, we need
		 *       to determine how to split the video in space and
		 *       time. */

		assert(i_vidSub.size());
		assert(i_vidSub[0].sz.whcf);
		assert(o_vid.sz.whcf);
		assert(o_vid.sz.frames   == i_vidSub[0].sz.frames  );
		assert(o_vid.sz.channels == i_vidSub[0].sz.channels);

		//! Determine width and height composition
		unsigned nW, nH;
		determineFactor(i_vidSub.size(), nW, nH);
		const unsigned hTmp = i_vidSub[0].sz.height - 2 * p_N;
		const unsigned wTmp = i_vidSub[0].sz.width  - 2 * p_N;

		//! Obtain inner image (containing boundaries)
		// TODO pending decision for video
		for (unsigned py = 0, n = 0; py < nH*hTmp; py += hTmp)
		for (unsigned px = 0       ; px < nW*wTmp; px += wTmp, n++)
		{
			/* Diagram for a 1D image with W = 8, covered
			 * with 2 sub images of w = 5, with border 2.
			 * Symmetrized pixels are indicated with an s.
			 *
			 * ori         0  1  2  3  4  5  6  7
			 * sub1 s0 s1  2  3  4  5  6 s7 s8
			 * sub2                s0 s1  2  3  4 s5 s6 s7 s8
			 * 
			 * px         0*w            1*w       <-- don't exceed 1*w + W-1*w
			 *
			 * Notation: [px,py] coords on big image of sub-image top-left point 
			 *           [qx,qy] point on big image
			 *           [sx,sy] corresponding point on sub-image
			 */
			unsigned wmax = std::min(wTmp, o_vid.sz.width  - px) + p_N;
			unsigned hmax = std::min(hTmp, o_vid.sz.height - py) + p_N;

			for (unsigned f = 0; f < o_vid.sz.frames  ; f++)
			for (unsigned c = 0; c < o_vid.sz.channels; c++)
			for (unsigned sy = p_N, qy = py; sy < hmax; sy++, qy++)
			for (unsigned sx = p_N, qx = px; sx < wmax; sx++, qx++)
				o_vid(qx, qy, f, c) = i_vidSub[n](sx, sy, f, c);
		}

		return;
	}
}

#endif // LIB_VIDEOT_HPP_INCLUDED
