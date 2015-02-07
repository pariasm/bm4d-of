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
 * @file LibVideo-utils.cpp
 * @brief Utilities extending the functionalities of the video container 
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

#include <math.h>
#include "Utilities.h"
#include "mt19937ar.h"

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
	,   const bool p_verbose
	){
		if (p_verbose) printf("Add noise [sigma = %f] ... ", p_sigma);

		//! Initialization
		o_vidNoisy = i_vid;
//		mt_init_genrand((unsigned long int) time (NULL) +
//		                (unsigned long int) getpid());
		mt_init_genrand(0); printf("Warning: random generator seed is 0 ");

		//! Add noise
		for (unsigned k = 0; k < i_vid.whcf; k++)
		{
			const double a = mt_genrand_res53();
			const double b = mt_genrand_res53();
			o_vidNoisy(k) += p_sigma *
				(float) (sqrtl(-2.0l * log(a)) * cos(2.0l * M_PI * b));
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
//	,   const bool p_verbose
	){
		if (i_vid1.size() != i_vid2.size())
		{
			fprintf(stderr, "Error @ LibVideo::computePsnr: "
					          "videos have different sizes:\n");
			fprintf(stderr, "\ti_vid1 %dx%dx%d - %d ch\n", i_vid1.width, i_vid1.height, i_vid1.frames, i_vid1.channels);
			fprintf(stderr, "\ti_vid2 %dx%dx%d - %d ch\n", i_vid2.width, i_vid2.height, i_vid2.frames, i_vid2.channels);
			return EXIT_FAILURE;
		}

		float sum = 0.f;
		for (unsigned k = 0; k < i_vid1.whcf; k++)
			sum += (i_vid1(k) - i_vid2(k)) * (i_vid1(k) - i_vid2(k));

		o_rmse = sqrtf(sum / (float) i_vid1.whcf);
		o_psnr = 20.f * log10f(255.f / o_rmse);

//		if (p_verbose)
//		{
//			std::cout << p_imageName << endl;
//			std::cout << "PSNR = " << o_psnr << endl;
//			std::cout << "RMSE = " << o_rmse << endl;
//		}

		return EXIT_SUCCESS;

	}
	
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
	,   const float p_min
	,   const float p_max
	){
		if (i_vid1.size() != i_vid2.size())
		{
			fprintf(stderr, "Error @ VideoUtils::computeDiff: "
			                "videos have different sizes\n");
			fprintf(stderr, "\ti_vid1 %dx%dx%d - %d ch\n", i_vid1.width, i_vid1.height, i_vid1.frames, i_vid1.channels);
			fprintf(stderr, "\ti_vid2 %dx%dx%d - %d ch\n", i_vid2.width, i_vid2.height, i_vid2.frames, i_vid2.channels);
			return EXIT_FAILURE;
		}

		o_vidDiff.resize(i_vid1.size());
		for (unsigned k = 0; k < i_vid1.whcf; k++)
		{
			float value =  (i_vid1(k) - i_vid2(k) + p_sigma) * p_max / (2.f * p_sigma);
			o_vidDiff(k) = clip(value, p_min, p_max);
		}

		return EXIT_SUCCESS;
	}
	
	/**
	 * @brief Add borders by symmetrizing
	 *
	 * @param i_vid1 : original video;
	 * @param o_vid2 : output with border added or removed;
	 * @param p_borderSize: size of added border;
	 * @param p_isForward: add border if true, else remove border;
	 * 
	 * @return none.
	 **/
	void addBorder(
		Video_f32 const &i_vid1
	,	Video_f32 &o_vid2
	,	const unsigned p_borderSize
	,	const bool p_isForward
	){
		//! Sizes
		const unsigned w2 = i_vid1.width  + (p_isForward ? 2*p_borderSize : - 2*p_borderSize);
		const unsigned h2 = i_vid1.height + (p_isForward ? 2*p_borderSize : - 2*p_borderSize);
		const unsigned f2 = i_vid1.frames + (p_isForward ? 2*p_borderSize : - 2*p_borderSize);
		const unsigned ch = i_vid1.channels;

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
	,	int p_origin_t
	,	int p_origin_x
	,	int p_origin_y
	){
		assert(o_vid2.channels == i_vid1.channels);

		//! Redefine invalid origin coordinates to default (centered crop)
		if (p_origin_t > (int)i_vid1.frames) p_origin_t = ((int)i_vid1.frames - (int)o_vid2.frames)/2;
		if (p_origin_x > (int)i_vid1.width ) p_origin_x = ((int)i_vid1.width  - (int)o_vid2.width )/2;
		if (p_origin_y > (int)i_vid1.height) p_origin_y = ((int)i_vid1.height - (int)o_vid2.height)/2;

		// TODO: more efficient implementation
		for (int      f = 0; f < o_vid2.frames  ; f++)
		for (unsigned c = 0; c < o_vid2.channels; c++)
		for (int      y = 0; y < o_vid2.height  ; y++)
		for (int      x = 0; x < o_vid2.width   ; x++)
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
	void crop(
		Video_f32 const &i_vid1
	,	Video_f32 &o_vid2
	,	const int * const p_origin
	){
		crop(i_vid1, o_vid2, p_origin[2], p_origin[0], p_origin[1]);
	}

	/**
	 * @brief Transform the color space of an video, from RGB to YUV, or vice-versa.
	 *
	 * @param io_vid: image on which the transform will be applied;
	 * @param p_isForward: if true, go from RGB to YUV, otherwise go from YUV to RGB.
	 *
	 * @return none.
	 **/
	void transformColorSpace(
		Video_f32& io_vid
	,	const bool p_isForward
	){
		//! If the image as only one channel, do nothing
		if (io_vid.channels == 1) return;

		//! Initialization
		const unsigned width  = io_vid.width;
		const unsigned height = io_vid.height;
		const unsigned chnls  = io_vid.channels;
		const unsigned wh     = io_vid.wh;

		for (int f = 0; f < io_vid.frames; f++)
			if (p_isForward) //< RGB to YUV
			{
				if (chnls == 3)
				{
					const unsigned red   = f * io_vid.whc;
					const unsigned green = red + wh;
					const unsigned blue  = green + wh;
					const float a = 1.f / sqrtf(3.f);
					const float b = 1.f / sqrtf(2.f);
					const float c = 2.f * a * sqrtf(2.f);

					float yuv[3];
					for (unsigned k = 0; k < wh; k++)
					{
						yuv[0] = a * (io_vid(k + red) + io_vid(k + green) + io_vid(k + blue));
						yuv[1] = b * (io_vid(k + red) - io_vid(k + blue));
						yuv[2] = c * (0.25f * io_vid(k + red ) - 0.5f * io_vid(k + green)
								      + 0.25f * io_vid(k + blue));

						io_vid(k + red  ) = yuv[0];
						io_vid(k + green) = yuv[1];
						io_vid(k + blue ) = yuv[2];
					}
				}
				else //< chnls == 4
				{
					const unsigned Gr = f * io_vid.whc;
					const unsigned R  = Gr + wh;
					const unsigned B  = R  + wh;
					const unsigned Gb = B  + wh;
					const float a = 0.5f;
					const float b = 1.f / sqrtf(2.f);

					float wtf[4];
					for (unsigned k = 0; k < wh; k++)
					{
						wtf[0] = a * ( io_vid(k + Gr) + io_vid(k + R ) +
						               io_vid(k + B ) + io_vid(k + Gb));
						wtf[1] = b * ( io_vid(k + R ) - io_vid(k + B ));
						wtf[2] = a * (-io_vid(k + Gr) + io_vid(k + R ) +
						               io_vid(k + B ) - io_vid(k + Gb));
						wtf[3] = b * (-io_vid(k + Gr) + io_vid(k + Gb));
						io_vid(k + Gr) = wtf[0];
						io_vid(k + R ) = wtf[1];
						io_vid(k + B ) = wtf[2];
						io_vid(k + Gb) = wtf[3];
					}
				}
			}
			else //< YUV to RGB
			{
				if (chnls == 3)
				{
					const unsigned red   = f * io_vid.whc;
					const unsigned green = red + wh;
					const unsigned blue  = green + wh;
					const float a = 1.f / sqrtf(3.f);
					const float b = 1.f / sqrtf(2.f);
					const float c = a / b;

					float rgb[3];
					for (unsigned k = 0; k < wh; k++)
					{
						rgb[0] = a * io_vid(k + red) + b * io_vid(k + green)
						                      + c * 0.5f * io_vid(k + blue );
						rgb[1] = a * io_vid(k + red) - c * io_vid(k + blue );
						rgb[2] = a * io_vid(k + red) - b * io_vid(k + green)
						                      + c * 0.5f * io_vid(k + blue );
						io_vid(k + red  ) = rgb[0];
						io_vid(k + green) = rgb[1];
						io_vid(k + blue ) = rgb[2];
					}
				}
				else //! chnls == 4
				{
					const unsigned Gr = f * io_vid.whc;
					const unsigned R  = Gr + wh;
					const unsigned B  = R  + wh;
					const unsigned Gb = B  + wh;
					const float a = 0.5f;
					const float b = 1.f / sqrtf(2.f);

					float wtf[4];
					for (unsigned k = 0; k < wh; k++)
					{
						wtf[0] = a * io_vid(k + Gr) - a * io_vid(k + B) - b * io_vid(k + Gb);
						wtf[1] = a * io_vid(k + Gr) + b * io_vid(k + R) + a * io_vid(k + B );
						wtf[2] = a * io_vid(k + Gr) - b * io_vid(k + R) + a * io_vid(k + B );
						wtf[3] = a * io_vid(k + Gr) - a * io_vid(k + B) + b * io_vid(k + Gb);
						io_vid(k + Gr) = wtf[0];
						io_vid(k + R ) = wtf[1];
						io_vid(k + B ) = wtf[2];
						io_vid(k + Gb) = wtf[3];
					}
				}
			}
	}
	
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
	// TODO pending decision for video
	int subDivide(
		Video_f32 const& i_im
	,	std::vector<Video_f32> &o_imSub
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
		const unsigned wTmp = ceil(float(i_im.width ) / float(nW)); // sizes w/out 
		const unsigned hTmp = ceil(float(i_im.height) / float(nH)); //     borders

		//! Obtain sub-images
		VideoSize imSubSize;
		imSubSize.width    = wTmp + 2 * p_N; // each sub-image has border
		imSubSize.height   = hTmp + 2 * p_N;
		imSubSize.frames   = i_im.frames; // NOTE: same frames as original
		imSubSize.channels = i_im.channels;
		imSubSize.update_fields();

		o_imSub.resize(p_nb);
		for (unsigned p = 0, n = 0; p < nH; p++)
		for (unsigned q = 0;        q < nW; q++, n++)
		{
			o_imSub[n].resize(imSubSize);

			// The origin is shifted -p_N to account for the subimage border
			int origin[3] = {q * wTmp - p_N, p * hTmp - p_N, 0};

			// Crop using symmetric boundary conditions
			VideoUtils::crop(i_im, o_imSub[n], origin);
		}
		return EXIT_SUCCESS;
	}

	/**
	 * @brief Reconstruct an video from its small sub-videos
	 *
	 * @param o_vid : image to reconstruct;
	 * @param i_vidSub : will contain all sub-images;
	 * @param p_N : boundary around sub-videos.
	 *
	 * @return EXIT_FAILURE in case of problems.
	 **/
	// TODO pending decision for video
	int subBuild(
	 	std::vector<Video_f32> const& i_vidSub
	,	Video_f32 &o_vid
	,	const unsigned p_N
	){
		/* FIXME current version builds a video that has been split
		 *       only spatially by subDivide. 
		 *       The reason is to mantain consistency with Marc's
		 *       code. For its proper extension to video, we need
		 *       to determine how to split the video in space and
		 *       time. */

		assert(i_vidSub.size());
		assert(i_vidSub[0].whcf);
		assert(o_vid.whcf);
		assert(o_vid.frames   == i_vidSub[0].frames  );
		assert(o_vid.channels == i_vidSub[0].channels);

		//! Determine width and height composition
		unsigned nW, nH;
		determineFactor(i_vidSub.size(), nW, nH);
		const unsigned hTmp = i_vidSub[0].height - 2 * p_N;
		const unsigned wTmp = i_vidSub[0].width  - 2 * p_N;

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
			unsigned wmax = std::min(wTmp, o_vid.width  - px) + p_N;
			unsigned hmax = std::min(hTmp, o_vid.height - py) + p_N;

			for (unsigned f = 0; f < o_vid.frames  ; f++)
			for (unsigned c = 0; c < o_vid.channels; c++)
			for (unsigned sy = p_N, qy = py; sy < hmax; sy++, qy++)
			for (unsigned sx = p_N, qx = px; sx < wmax; sx++, qx++)
				o_vid(qx, qy, f, c) = i_vidSub[n](sx, sy, f, c);
		}

		return EXIT_SUCCESS;
	}
}
