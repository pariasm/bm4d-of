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
			fprintf(stderr, "\ti_vid1 %dx%dx%d - %d ch\n", i_vid1.width, i_vid1.height, i_vid1.nFrames, i_vid1.nChannels);
			fprintf(stderr, "\ti_vid2 %dx%dx%d - %d ch\n", i_vid2.width, i_vid2.height, i_vid2.nFrames, i_vid2.nChannels);
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
	,   const bool p_verbose
	){
		if (i_vid1.size() != i_vid2.size())
		{
			fprintf(stderr, "Error @ VideoUtils::computeDiff: "
			                "videos have different sizes\n");
			fprintf(stderr, "\ti_vid1 %dx%dx%d - %d ch\n", i_vid1.width, i_vid1.height, i_vid1.nFrames, i_vid1.nChannels);
			fprintf(stderr, "\ti_vid2 %dx%dx%d - %d ch\n", i_vid2.width, i_vid2.height, i_vid2.nFrames, i_vid2.nChannels);
			return EXIT_FAILURE;
		}

		if (p_verbose) printf("Compute difference... ");
		o_vidDiff.resize(i_vid1.size());
		for (unsigned k = 0; k < i_vid1.whcf; k++)
		{
			float value =  (i_vid1(k) - i_vid2(k) + p_sigma) * p_max / (2.f * p_sigma);
			o_vidDiff(k) = clip(value, p_min, p_max);
		}

		if (p_verbose) printf("done.");
		return EXIT_SUCCESS;
	}
	

	/**
	 * @brief Add boundary by symmetry.
	 *
	 * @param i_vid : video to symmetrize;
	 * @param o_vidSym : will contain i_vid with symmetrized boundaries;
	 *
	 * @return none.
	 **
	int addBoundary(
		std::vector<float> const& i_vid
	,	std::vector<float> &o_vidSym
	){
		//! Parameters declarations
		const unsigned width  = p_vid.width;
		const unsigned height = p_vid.height;
		const unsigned chnls  = p_vid.nChannels;
		const unsigned h      = p_vidSym.height;
		const unsigned w      = p_vidSym.width;

		if (w < width || h < height)
		{
			fprintf(stderr, "Error @ VideoUtils::addBoundary: symmetrized"
			                "video has to be larger than original\n");
			return EXIT_FAILURE;
		}

		//! Resize output image if nenessary
		if (o_imSym.size() != chnls * h * w) o_imSym.resize(chnls * w * h);

		//! Declaration
		for (unsigned c = 0; c < chnls; c++)
		{
			const unsigned dc1 = c * width * height;
			const unsigned dc2 = c * w * h;

			//! Center of the image
			for (unsigned i = 0; i < height; i++)
			for (unsigned j = 0; j < width ; j++)
				o_imSym[dc2 + i * w + j] = i_im[dc1 + i * width + j];

			//! Right
			for (unsigned i = 0    ; i < height; i++)
			for (unsigned j = width; j < w     ; j++)
				o_imSym[dc2 + i * w + j] = o_imSym[dc2 + i * w + 2 * width - j - 1];

			//! Bottom
			for (unsigned i = height; i < h; i++)
			for (unsigned j = 0     ; j < w; j++)
				o_imSym[dc2 + i * w + j] = o_imSym[dc2 + (2 * height - i - 1) * w + j];
		}

		return EXIT_SUCCESS;
	}*/

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
	){
	}
	
	/**
	 * @brief Add boundaries by symmetry
	 *
	 * @param io_vid : original image;
	 * @param io_vidSym : contain io_im symmetrized;
	 *
	 * @return none.
	 **/
	// TODO: make more general: we will need different sizes for different
	//       borders. More general inputs could be:
	//				 - vid_in   : input video
	//				 - vid_out  : output video
	//				 - size_out : size of output video
	//				 - offset   : position of vid2 origin in vid1 coordinates
	void symmetrizeVideo(
		Video_f32 const &i_vid1
	,	Video_f32 &o_vid2
	,	const unsigned p_borderSize
	,	const bool p_isForward
	){
		//! Sizes
		const unsigned w1 = i_vid1.width;
		const unsigned h1 = i_vid1.height;
		const unsigned f1 = i_vid1.nFrames;
		const unsigned w2 = w1 + (p_isForward ? 2*p_borderSize : - 2*p_borderSize);
		const unsigned h2 = h1 + (p_isForward ? 2*p_borderSize : - 2*p_borderSize);
		const unsigned f2 = f1 + (p_isForward ? 2*p_borderSize : - 2*p_borderSize);
		const unsigned ch = i_vid1.nChannels;

		//! Position of vid2 origin in vid1 coordinates
		const int tx = p_isForward ? -p_borderSize : p_borderSize;
		const int ty = p_isForward ? -p_borderSize : p_borderSize;
		const int tf = p_isForward ? -p_borderSize : p_borderSize;

		//! Resize output image, if necessary
		o_vid2.resize(w2, h2, f2, ch);

		// TODO: more efficient implementation
		for (int      f = 0; f < f2; f++)
		for (unsigned c = 0; c < ch; c++)
		for (int      y = 0; y < h2; y++)
		for (int      x = 0; x < w2; x++)
			o_vid2(x,y,f,c) = 
				i_vid1.getPixelSymmetric(x + tx,y + ty,f + tf,c);

		return;
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
		if (io_vid.nChannels == 1) return;

		//! Initialization
		const unsigned width  = io_vid.width;
		const unsigned height = io_vid.height;
		const unsigned chnls  = io_vid.nChannels;
		const unsigned wh     = io_vid.wh;

		for (int f = 0; f < io_vid.nFrames; f++)
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
	int subDivide(
		Video_f32 const& i_im
	,	std::vector<Video_f32> &o_imSub
	,	VideoSize &p_imSizeSub
	,	const unsigned p_N
	,	const unsigned p_nb
	){
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
	int subBuild(
		Video_f32 &o_vid
	,	std::vector<Video_f32> const& i_vidSub
	,	const unsigned p_N
	){
	}
}
