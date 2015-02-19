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
#include "VideoNLMeans.h"
#include "cmd_option.h"

#include <algorithm>

/**
 * @file   testVideoNLMeans.cpp
 * @brief  Executable file to test VideoNLMeans
 *
 *
 *
 * @author PABLO ARIAS  <pariasm@gmail.com>
 **/

void print_patch_from_group(
	std::vector<std::vector<float> > &group3d
,	VideoNLM::nlbParams prms
,	unsigned idx
){
	const unsigned sP    = prms.sizePatch;
	const unsigned nSimP = prms.nSimilarPatches;

	assert(idx < nSimP);


	for (unsigned hy = 0, k = 0; hy < sP; hy++)
	{
		for (unsigned hx = 0; hx < sP; hx++, k++)
		{
			printf("% 3f,% 3f,% 3f   ", group3d[0][k*nSimP + idx], 
			                            group3d[1][k*nSimP + idx],
			                            group3d[2][k*nSimP + idx]);
		}
		printf("\n");
	}
}

void compare_patch_groups(
	std::vector<std::vector<float> > &group3d
,	std::vector<std::pair<float, unsigned> > &patch_dists
,	VideoNLM::nlbParams prms
,	VideoSize sz
,	unsigned idx
){
	assert(idx < patch_dists.size());

	const unsigned w  = sz.width;
	const unsigned wh = sz.wh;

	unsigned px, py, pt, pc;
	sz.coords(patch_dists[0].second, px, py, pt, pc);

	printf("patch %d from stack, anchored at: [%d, %d, %d]:\n", 0, px, py, pt);
	print_patch_from_group(group3d, prms, 0);

	sz.coords(patch_dists[idx].second, px, py, pt, pc);
	printf("patch %d from stack, anchored at: [%d, %d, %d]:\n", idx, px, py, pt);
	print_patch_from_group(group3d, prms, idx);
	printf("\ndistance between them: %f\n\n\n", patch_dists[idx].first);
}

namespace VideoNLM
{
	std::vector<std::pair<float, unsigned> > estimateSimilarPatchesStep1_debug(
		Video<float> const& i_im
	,	std::vector<std::vector<float> > &o_group3d
	,	std::vector<unsigned> &o_index
	,	const unsigned pidx
	,	const nlbParams &p_params
	){
		//! Initialization
		const unsigned sWx   = p_params.sizeSearchWindow;
		const unsigned sWt_f = p_params.sizeSearchTimeRangeFwd;
		const unsigned sWt_b = p_params.sizeSearchTimeRangeBwd;
		const unsigned sP    = p_params.sizePatch;
		const unsigned nSimP = p_params.nSimilarPatches;
	
		//! Coordinates and index of top-left-back pixel of search box
		unsigned px, py, pt, pc;
		i_im.sz.coords(pidx, px, py, pt, pc);
	
		const unsigned rangex[2] = {px - (sWx - 1)/2, px + (sWx - 1)/2};
		const unsigned rangey[2] = {py - (sWx - 1)/2, py + (sWx - 1)/2};
		const unsigned ranget[2] = {std::max((int)pt - (int)sWt_b, 0),
		                            std::min((int)pt + (int)sWt_f, (int)i_im.sz.frames-1)};

		std::vector<std::pair<float, unsigned> > distance(sWx * sWx * 
		                                            (ranget[1] - ranget[0] + 1));
	
		//! Compute distance between patches in search range
		for (unsigned qt = ranget[0], dt = 0; qt <= ranget[1]; qt++, dt++)
		for (unsigned qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
		for (unsigned qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
		{
			//! Squared L2 distance
			float diff = 0.f;
			for (unsigned hy = 0; hy < sP; hy++)
			for (unsigned hx = 0; hx < sP; hx++)
			{
				const float tmp = i_im(px + hx, py + hy, pt)
				                - i_im(qx + hx, qy + hy, qt);
				diff += tmp * tmp;
			}
	
			//! Save distance and corresponding patch index
			distance[dt * sWx*sWx + dy * sWx + dx] = 
				std::make_pair(diff, i_im.sz.index(qx, qy, qt, 0));
		}
	
		//! Keep only the N2 best similar patches
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
		                  distance.end(), comparaisonFirst);
	
		//! Register position of patches
		for (unsigned n = 0; n < nSimP; n++) o_index[n] = distance[n].second;
	
		//! Stack selected patches into the 3D group
		const unsigned w  = i_im.sz.width;
		const unsigned wh = i_im.sz.wh;
		for (unsigned c  = 0; c < i_im.sz.channels; c++)
		for (unsigned hy = 0, k = 0; hy < sP; hy++)
		for (unsigned hx = 0;        hx < sP; hx++)
		for (unsigned n  = 0; n < nSimP; n++, k++)
			o_group3d[c][k] = i_im.data[o_index[n] + hy * w + hx + c * wh];

		return distance;
	
		/* 00  pixels from all patches
		 * 01  pixels from all patches
		 * ...
		 * 0sp pixels from all patches
		 */
	}

	void computeNlMeansEstimateStep1_debug(
		std::vector<std::vector<float> > &io_group3d
//	,	matParams &i_mat
	,	unsigned &io_nInverseFailed
	,	nlbParams const& p_params
	){
		//! Parameters
		const unsigned chnls = io_group3d.size();
		const unsigned nSimP = p_params.nSimilarPatches;
		const unsigned sP2   = p_params.sizePatch * p_params.sizePatch;
		const float factor   = 1.f/(2.f*p_params.beta*p_params.beta);

		std::vector<std::vector<float> > tmp_group3d = io_group3d; 

//		const unsigned sP = p_params.sizePatch;
//		for (unsigned i = nSimP -3; i < nSimP; i++)
//		{
//			for (unsigned hy = 0, k = 0; hy < sP; hy++)
//			{
//				for (unsigned hx = 0; hx < sP; hx++, k++)
//					printf("% 3f   ", io_group3d[0][k*nSimP + i]); 
//				printf("\n");
//			}
//			printf("\n");
//		}

		//! Vector of weight normalization factors
		std::vector<float> total_weights(nSimP, 1.f);
		//! We initialize with 1 since for the case j = i, weight = 1

		//! For each patch in patch group
		for (unsigned i = 0; i < nSimP; i++)
		{
			//! Accumulate contributions from other patches in group
			for (unsigned j = 0; j < i+1; j++)
				printf("%8.0f ", 0.f);

			//! Accumulate contributions from other patches in group
			for (unsigned j = i+1; j < nSimP; j++)
			{
				//! Compute distance between Y components of patches
				float dist = 0.f, dif;
				for (unsigned h = 0; h < sP2; h++)
					dist += (dif = io_group3d[0][h*nSimP + i] - io_group3d[0][h*nSimP + j]) * dif ;

				//! Unnormalized exponential weight
				const float weight = exp(- factor * dist / (float)sP2 );
				total_weights[i] += weight;
				total_weights[j] += weight;

				//! Accumulate nl-mean
				for (unsigned c = 0; c < chnls; c++)
				for (unsigned h = 0; h < sP2  ; h++)
				{
					tmp_group3d[c][h*nSimP + i] += weight * io_group3d[c][h*nSimP + j];
					tmp_group3d[c][h*nSimP + j] += weight * io_group3d[c][h*nSimP + i];
				}

				printf("%8.4f ", tmp_group3d[0][i]);
			}

			printf("%f \n", total_weights[i]);

			//! Normalize nl-mean while copying to output stack
			const float factor = 1.f/total_weights[i];
			for (unsigned c = 0; c < chnls; c++)
			for (unsigned h = 0; h < sP2  ; h++)
				io_group3d[c][h*nSimP + i] = tmp_group3d[c][h*nSimP + i] * factor;
		}
		io_nInverseFailed = 0;
	}

	void computeNlMeansEstimateStep1_debug_ref(
		std::vector<std::vector<float> > &io_group3d
//	,	matParams &i_mat
	,	unsigned &io_nInverseFailed
	,	nlbParams const& p_params
	){
		//! Parameters
		const unsigned chnls = io_group3d.size();
		const unsigned nSimP = p_params.nSimilarPatches;
		const unsigned sP2   = p_params.sizePatch * p_params.sizePatch;
		const float factor   = 1.f/(2.f*p_params.beta*p_params.beta);

		std::vector<std::vector<float> > tmp_group3d(io_group3d.size(),
		                     std::vector<float>( io_group3d[0].size(), 0.f)); 

//		const unsigned sP = p_params.sizePatch;
//		for (unsigned i = nSimP -3; i < nSimP; i++)
//		{
//			for (unsigned hy = 0, k = 0; hy < sP; hy++)
//			{
//				for (unsigned hx = 0; hx < sP; hx++, k++)
//					printf("% 3f   ", io_group3d[0][k*nSimP + i]); 
//				printf("\n");
//			}
//			printf("\n");
//		}

		//! Vector of weight normalization factors
		std::vector<float> total_weights(nSimP, 0.f);

		//! For each patch in patch group
		for (unsigned i = 0; i < nSimP; i++)
		{
			//! Accumulate contributions from other patches in group
			for (unsigned j = 0; j < nSimP; j++)
			{
				//! Compute distance between Y components of patches
				float dist = 0.f, dif;
				for (unsigned h = 0; h < sP2; h++)
					dist += (dif = io_group3d[0][h*nSimP + i] - io_group3d[0][h*nSimP + j]) * dif ;

				//! Unnormalized exponential weight
				const float weight = exp(- factor * dist / (float)sP2);
				total_weights[i] += weight;

				//! Accumulate nl-mean
				for (unsigned c = 0; c < chnls; c++)
				for (unsigned h = 0; h < sP2  ; h++)
					tmp_group3d[c][h*nSimP + i] += weight * io_group3d[c][h*nSimP + j];

				printf("%8.4f ", tmp_group3d[0][i]);
			}

			printf("%f \n", total_weights[i]);
		}

		//! Normalize nl-mean while copying to output stack
		for (unsigned i = 0; i < nSimP; i++)
		{
			const float factor = 1.f/total_weights[i];
			for (unsigned c = 0; c < chnls; c++)
			for (unsigned h = 0; h < sP2  ; h++)
				io_group3d[c][h*nSimP + i] = tmp_group3d[c][h*nSimP + i] * factor;
		}

		io_nInverseFailed = 0;
	}
}


void print_video_size(const std::string& name, const Video<float>& vid)
{
	printf("%s", name.c_str());
	printf("\tFrames: %d - Channels: %d - Height: %d - Width %d\n", 
			vid.sz.frames, vid.sz.channels, vid.sz.height, vid.sz.width);
}

int main(int argc, char **argv)
{
	//! Check if there is the right call for the algorithm
	clo_usage("Unit test for Video NL-Bayes methods");

	//! Get command line inputs
	const string i_video_path = clo_option("-i", "", "Patch to input sequence (printf format)");
	const int i_firstFrame = clo_option("-f", 0, "First frame");
	const int i_lastFrame  = clo_option("-l", 0, "Last frame");
	const int i_frameStep  = clo_option("-s", 1, "Frame step");
	const int i_sigma      = clo_option("-sigma", 0, "Add noise of standard deviation sigma");
	const unsigned px      = clo_option("-px", 0, "Current point of analysis");
	const unsigned py      = clo_option("-py", 0, "Current point of analysis");
	const unsigned pt      = clo_option("-pt", 0, "Current point of analysis");

	//! Load video
	Video<float> vid_ori;
	{
		vid_ori.loadVideo(i_video_path.c_str(), i_firstFrame, i_lastFrame, i_frameStep);

		print_video_size("loaded video1", vid_ori);
	}

	/*/! Test computeNlMeansEstimateStep1
	{
		using std::vector;

		//! RGB to YUV
		Video<float> vid(vid_ori);
		VideoUtils::transformColorSpace(vid, true);

		//! Initialize parameter structures
		VideoNLM::nlbParams prms1, prms2;
		VideoNLM::initializeNlbParameters(prms1, 1, i_sigma, vid.sz, 0, 1);
		VideoNLM::initializeNlbParameters(prms2, 2, i_sigma, vid.sz, 0, 1);
		VideoNLM::printNlbParameters(prms1);
		VideoNLM::printNlbParameters(prms2);

		//! Used matrices during Bayes' estimate
		const unsigned patch_dim = prms1.sizePatch * prms1.sizePatch ;
		const unsigned patch_num = prms1.nSimilarPatches ;

		//! Matrices used for Bayes' estimate
		vector<unsigned> patch_index(patch_num);
		vector<vector<float> > patch_stack(vid.sz.channels,
				                             vector<float>(patch_num * patch_dim));

		//! Point
		unsigned pind = vid.sz.index(px, py, pt, 0);

		//! Compute stack of patches (detailed output)
		{
			printf("computing patch distances for point: [% 3d,% 3d,% 2d]\n\n",px,py,pt);

			//! Compute distances
			vector<std::pair<float, unsigned> > patch_dists = 
				VideoNLM::estimateSimilarPatchesStep1_debug(vid, patch_stack, patch_index, pind, prms1);

//			//! Display seleted patches
//			for (int n = 0; n < prms1.nSimilarPatches; n++)
//				compare_patch_groups(patch_stack, patch_dists, prms1, vid.sz, n);

			//! Show distances in new image
			Video<float> dist(vid.sz.width, vid.sz.height, vid.sz.frames, 1, -1);
			Video<float> idxs(vid.sz.width, vid.sz.height, vid.sz.frames, 1, -1);

			for (int n = 0; n < patch_dists.size(); n++)
			{
				unsigned i3 = patch_dists[n].second;
				unsigned i1 = (i3 / vid.sz.whc) * vid.sz.wh + i3 % vid.sz.wh;
				dist(i1) = patch_dists[n].first;
				idxs(i1) = n;
			}

			//! Save distance video
			char name[1024];
			sprintf(name,"/tmp/dist_to_%03d.%03d.%02d",px,py,pt);
			dist.saveVideoAscii(name, i_firstFrame, i_frameStep);
			sprintf(name,"/tmp/idxs_to_%03d.%03d.%02d",px,py,pt);
			idxs.saveVideoAscii(name, i_firstFrame, i_frameStep);

			//! Save YUV components
			vid.saveVideoAscii("/tmp/vid_yuv", i_firstFrame, i_frameStep);
		}

		//! Compute nlmeans estimate
		{
			std::vector<std::vector<float> > patch_stack1 = patch_stack;
			std::vector<std::vector<float> > patch_stack2 = patch_stack;

			unsigned n_failed = 0;
			computeNlMeansEstimateStep1_debug    (patch_stack1, n_failed, prms1);
			computeNlMeansEstimateStep1_debug_ref(patch_stack2, n_failed, prms1);

			for (int n = 0; n < prms1.nSimilarPatches; n++)
			{
				printf("\nOriginal patch stack:\n");
				print_patch_from_group(patch_stack, prms1, n);

				printf("\nOptimized implementation:\n");
				print_patch_from_group(patch_stack1, prms1, n);

				printf("\nReference (non-optimized) implementation:\n");
				print_patch_from_group(patch_stack2, prms1, n);
			}
		}

		//! Aggregate 
		Video<float> basic(vid.sz);
		Video<float> weight(vid.sz, 0.f);
		Video<char > mask(vid.sz.width, vid.sz.height, vid.sz.frames, 1, true);
		VideoNLM::computeAggregationStep1(basic, weight, mask, patch_stack, patch_index, prms1);

		//! Save outputs to disk
		basic .saveVideoAscii("/tmp/basic" , i_firstFrame, i_frameStep);
		weight.saveVideoAscii("/tmp/weight", i_firstFrame, i_frameStep);
		
		Video<float> mask_out(vid.sz.width, vid.sz.height, vid.sz.frames);
		for (unsigned i = 0; i < mask.sz.whcf; i++)
			mask_out(i) = 255.*(float)mask(i);
		mask_out.saveVideo("/tmp/mask_%02d.png", i_firstFrame, i_frameStep);

		return EXIT_SUCCESS;

	}//*/

	//! Test processNlBayes, both steps
	{
		//! Add noise
		Video<float> vid;
		VideoUtils::addNoise(vid_ori, vid, i_sigma, true);

		//! Save noisy video
		vid.saveVideo("/tmp/vid_%02d.png", i_firstFrame, i_frameStep);

		//! Initialize parameter structures
		VideoNLM::nlbParams prms1, prms2;
		VideoNLM::initializeNlbParameters(prms1, 1, i_sigma, vid.sz, 0, 1);
		VideoNLM::initializeNlbParameters(prms2, 2, i_sigma, vid.sz, 0, 1);
		VideoNLM::printNlbParameters(prms1);
		VideoNLM::printNlbParameters(prms2);

		//! Parallelization: number of subvideos in which video is divided
		const unsigned nParts = 2;

		//! Run first step
		Video<float> basic(vid.sz);
		{
			//! RGB to YUV
			VideoUtils::transformColorSpace(vid, true);

			//! Divide in subvideos
			std::vector<Video<float> > parts_vid(nParts);
			VideoUtils::subDivide(vid, parts_vid, prms1.boundary, nParts);

			//! Process all sub-videos
			std::vector<Video<float> > parts_basic(nParts);
			std::vector<Video<float> > parts_final(nParts);
			for (int n = 0; n < (int)nParts; n++)
				VideoNLM::processNlBayes(parts_vid[n], parts_basic[n], parts_final[n], prms1);

			//! Get the basic estimate
			VideoUtils::subBuild(parts_basic, basic, prms1.boundary);

			//! YUV to RGB
			VideoUtils::transformColorSpace(vid  , false);
			VideoUtils::transformColorSpace(basic, false);
		}

		//! Save video
		basic.saveVideo("/tmp/basic_%02d.png", i_firstFrame, i_frameStep);

		//! Run second step
		Video<float> final(vid.sz);
		{
			//! Divide in subvideos
			std::vector<Video<float> > parts_vid  (nParts);
			std::vector<Video<float> > parts_basic(nParts);
			VideoUtils::subDivide(vid  , parts_vid  , prms2.boundary, nParts);
			VideoUtils::subDivide(basic, parts_basic, prms2.boundary, nParts);

			//! Process all sub-videos
			std::vector<Video<float> > parts_final(nParts);
			for (int n = 0; n < (int)nParts; n++)
				VideoNLM::processNlBayes(parts_vid[n], parts_basic[n], parts_final[n], prms2);

			//! Get the basic estimate
			VideoUtils::subBuild(parts_final, final, prms2.boundary);
		}

		//! Save video
		final.saveVideo("/tmp/final_%02d.png", i_firstFrame, i_frameStep);
		//*/

		return EXIT_SUCCESS;
	}//*/

}
