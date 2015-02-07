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
#include "VideoNLBayes.h"
#include "cmd_option.h"

#include <algorithm>

/**
 * @file   testVideoNLBayes.cpp
 * @brief  Executable file to test VideoNLBayes
 *
 *
 *
 * @author PABLO ARIAS  <pariasm@gmail.com>
 **/

namespace VideoNLB
{
	std::vector<std::pair<float, unsigned> > estimateSimilarPatchesStep1_debug(
		Video_f32 const& i_im
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
		i_im.coords(pidx, px, py, pt, pc);
	
		const unsigned rangex[2] = {px - (sWx - 1)/2, px + (sWx - 1)/2};
		const unsigned rangey[2] = {py - (sWx - 1)/2, py + (sWx - 1)/2};
		const unsigned ranget[2] = {std::max((int)pt - (int)sWt_b, 0),
		                            std::min((int)pt + (int)sWt_f, (int)i_im.frames-1)};

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
				std::make_pair(diff, i_im.index(qx, qy, qt, 0));
		}
	
		//! Keep only the N2 best similar patches
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
		                  distance.end(), comparaisonFirst);
	
		//! Register position of patches
		for (unsigned n = 0; n < nSimP; n++) o_index[n] = distance[n].second;
	
		//! Stack selected patches into the 3D group
		const unsigned w  = i_im.width;
		const unsigned wh = i_im.wh;
		for (unsigned c  = 0; c < i_im.channels; c++)
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
}

void print_patch_group(
	std::vector<std::vector<float> > &group3d
,	std::vector<std::pair<float, unsigned> > &patch_dists
,	VideoNLB::nlbParams prms
,	VideoSize sz
,	unsigned idx
){
	assert(idx < patch_dists.size());

	const unsigned sP    = prms.sizePatch;
	const unsigned nSimP = prms.nSimilarPatches;
	const unsigned w  = sz.width;
	const unsigned wh = sz.wh;

	unsigned px, py, pt, pc;
	sz.coords(patch_dists[0].second, px, py, pt, pc);

	printf("patch %d from stack, anchored at: [%d, %d, %d]:\n", 0, px, py, pt);
	for (unsigned hy = 0, k = 0; hy < sP; hy++)
	{
		for (unsigned hx = 0; hx < sP; hx++, k++)
		{
			printf("% 3f,% 3f,% 3f   ", group3d[0][k*nSimP + 0], 
			                            group3d[1][k*nSimP + 0],
			                            group3d[2][k*nSimP + 0]);
		}
		printf("\n");
	}
	printf("\n");

	sz.coords(patch_dists[idx].second, px, py, pt, pc);
	printf("patch %d from stack, anchored at: [%d, %d, %d]:\n", idx, px, py, pt);
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
	printf("\ndistance between them: %f\n\n\n", patch_dists[idx].first);
}



void print_video_size(const std::string& name, const Video_f32& vid)
{
	printf("%s", name.c_str());
	printf("\tFrames: %d - Channels: %d - Height: %d - Width %d\n", 
			vid.frames, vid.channels, vid.height, vid.width);
}

int main(int argc, char **argv)
{
	//! Check if there is the right call for the algorithm
	clo_usage("Unit test for Video NL-Bayes methods");
//	if (argc != 9 && argc != 6)
//	{
//		fprintf(stdout, "Usage: %s path-to-frames first-frame last-frame frame-step sigma [px py pt]\n", argv[0]);
//		return EXIT_FAILURE;
//	}

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
	Video_f32 vid_ori;
	{
		if (vid_ori.loadVideo(i_video_path.c_str(), i_firstFrame, i_lastFrame, i_frameStep)
				== EXIT_FAILURE)
		{
			fprintf(stderr, "Exiting. Failed to load video.\n");
			return EXIT_FAILURE;
		}

		print_video_size("loaded video1", vid_ori);
	}

	/*/! Test estimateSimilarPatchesStep1 and computeAggregationStep1 on a single pixel
	{
		using std::vector;

		//! RGB to YUV
		Video_f32 vid(vid_ori);
		VideoUtils::transformColorSpace(vid, true);

		//! Initialize parameter structures
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, prms2, i_sigma, vid.size(), 0, 0, 1);
		VideoNLB::printNlbParameters(prms1, prms2);

		//! Used matrices during Bayes' estimate
		const unsigned patch_dim = prms1.sizePatch * prms1.sizePatch ;
		const unsigned patch_num = prms1.nSimilarPatches ;

		//! Matrices used for Bayes' estimate
		vector<unsigned> patch_index(patch_num);
		vector<vector<float> > patch_stack(vid.channels,
				                             vector<float>(patch_num * patch_dim));

		//! Compute stack of patches (detailed output)
		{
			printf("computing patch distances for point: [% 3d,% 3d,% 2d]\n\n",px,py,pt);

			//! Point
			unsigned pind = vid.index(px, py, pt, 0);

			//! Compute distances
			vector<std::pair<float, unsigned> > patch_dists = 
				VideoNLB::estimateSimilarPatchesStep1_debug(vid, patch_stack, patch_index, pind, prms1);

			//! Display seleted patches
			for (int n = 0; n < prms1.nSimilarPatches; n++)
				print_patch_group(patch_stack, patch_dists, prms1, vid.size(), n);

			//! Show distances in new image
			Video_f32 dist(vid.width, vid.height, vid.frames, 1, -1);
			Video_f32 idxs(vid.width, vid.height, vid.frames, 1, -1);

			for (int n = 0; n < patch_dists.size(); n++)
			{
				unsigned i3 = patch_dists[n].second;
				unsigned i1 = (i3 / vid.whc) * vid.wh + i3 % vid.wh;
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

		//! Aggregate 
		Video_f32 basic(vid.size());
		Video_f32 weight(vid.size(), 0.f);
		std::vector<bool> mask(vid.whf, true);
		VideoNLB::computeAggregationStep1(basic, weight, mask, patch_stack, patch_index, prms1);

		//! Save outputs to disk
		basic .saveVideoAscii("/tmp/basic" , i_firstFrame, i_frameStep);
		weight.saveVideoAscii("/tmp/weight", i_firstFrame, i_frameStep);
		
		Video_f32 mask_out(vid.width, vid.height, vid.frames);
		for (unsigned i = 0; i < mask.size(); i++)
			mask_out(i) = 255.*(float)mask[i];
		mask_out.saveVideo("/tmp/mask_%02d.png", i_firstFrame, i_frameStep);

	}//*/

	/*/! Test estimateSimilarPatchesStep1 and computeAggregationStep1 on all video
	{
		using std::vector;

		//! RGB to YUV
		Video_f32 vid(vid_ori);
		VideoUtils::transformColorSpace(vid, true);

		//! Initialize parameter structures
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, prms2, i_sigma, vid.size(), 0, 0, 1);
		VideoNLB::printNlbParameters(prms1, prms2);

		//! Used matrices during Bayes' estimate
		const unsigned patch_dim = prms1.sizePatch * prms1.sizePatch ;
		const unsigned patch_num = prms1.nSimilarPatches ;

		//! Matrices used for Bayes' estimate
		vector<unsigned> patch_index(patch_num);
		vector<vector<float> > patch_stack(vid.channels,
				                             vector<float>(patch_num * patch_dim));

		Video_f32 basic(vid.size());
		Video_f32 weight(vid.size(), 0.f);
		std::vector<bool> mask(vid.whf, false);

		//! Only pixels of the center of the image must be processed (not the boundaries)
		unsigned sWx = prms1.sizeSearchWindow;
		for (unsigned f =   0; f < vid.frames      ; f++)
		for (unsigned i = sWx; i < vid.height - sWx; i++)
		for (unsigned j = sWx; j < vid.width  - sWx; j++)
			mask[f * vid.wh + i * vid.width + j] = true;


		//! Aggregate
		for (int pt = 0; pt < vid.frames; pt++)
		for (int py = 0; py < vid.height; py++)
		for (int px = 0; px < vid.width ; px++)
		if  (mask[pt * vid.wh + py * vid.width + px])
		{
			//! Point
			unsigned pind = vid.index(px, py, pt, 0);

			//! Compute stack of patches
			VideoNLB::estimateSimilarPatchesStep1(vid, patch_stack, patch_index, pind, prms1);

			//! Aggregate 
			VideoNLB::computeAggregationStep1(basic, weight, mask, patch_stack, patch_index, prms1);
		}

		//! Normalize
		for (int pt = 0; pt < vid.frames  ; pt++)
		for (int pc = 0; pc < vid.channels; pc++)
		for (int py = 0; py < vid.height  ; py++)
		for (int px = 0; px < vid.width   ; px++)
			if (weight(px,py,pt,pc))
				basic(px,py,pt,pc) /= weight(px,py,pt,pc);

		//! YUV to RGB
		VideoUtils::transformColorSpace(basic, false);

		//! Save outputs to disk
		basic .saveVideo("/tmp/basic_%02d.png" , i_firstFrame, i_frameStep);
		weight.saveVideoAscii("/tmp/weight", i_firstFrame, i_frameStep);
		
		Video_f32 mask_out(vid.width, vid.height, vid.frames);
		for (unsigned i = 0; i < mask.size(); i++)
			mask_out(i) = 255.*(float)mask[i];
		mask_out.saveVideo("/tmp/mask_%02d.png", i_firstFrame, i_frameStep);

	}//*/

	//! Test processNlBayes, both steps
	{
		//! Add noise
		Video_f32 vid;
		VideoUtils::addNoise(vid_ori, vid, i_sigma, true);

		//! Save noisy video
		vid.saveVideo("/tmp/vid_%02d.png", i_firstFrame, i_frameStep);

		//! Initialize parameter structures
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, prms2, i_sigma, vid.size(), 0, 0, 1);
		VideoNLB::printNlbParameters(prms1, prms2);

		//! Parallelization: number of subvideos in which video is divided
		const unsigned nParts = 2;

		//! Run first step
		Video_f32 basic(vid.size());
		{
			//! RGB to YUV
			VideoUtils::transformColorSpace(vid, true);

			//! Divide in subvideos
			std::vector<Video_f32> parts_vid(nParts);
			VideoUtils::subDivide(vid, parts_vid, prms1.boundary, nParts);

			//! Process all sub-videos
			std::vector<Video_f32> parts_basic(nParts);
			std::vector<Video_f32> parts_final(nParts);
			for (int n = 0; n < (int)nParts; n++)
				processNlBayes(parts_vid[n], parts_basic[n], parts_final[n], prms1);

			//! Get the basic estimate
			VideoUtils::subBuild(parts_basic, basic, prms1.boundary);

			//! YUV to RGB
			VideoUtils::transformColorSpace(vid  , false);
			VideoUtils::transformColorSpace(basic, false);
		}

		//! Run second step
		Video_f32 final(vid.size());
		{
			//! Divide in subvideos
			std::vector<Video_f32> parts_vid  (nParts);
			std::vector<Video_f32> parts_basic(nParts);
			VideoUtils::subDivide(vid  , parts_vid  , prms2.boundary, nParts);
			VideoUtils::subDivide(basic, parts_basic, prms2.boundary, nParts);

			//! Process all sub-videos
			std::vector<Video_f32> parts_final(nParts);
			for (int n = 0; n < (int)nParts; n++)
				processNlBayes(parts_vid[n], parts_basic[n], parts_final[n], prms2);

			//! Get the basic estimate
			VideoUtils::subBuild(parts_final, final, prms2.boundary);
		}

		//! Save video
		basic.saveVideo("/tmp/basic_%02d.png", i_firstFrame, i_frameStep);
		final.saveVideo("/tmp/final_%02d.png", i_firstFrame, i_frameStep);

		return EXIT_SUCCESS;
	}//*/

	/*/! Test indexing operations
	{
		VideoSize size3 = vid_ori.size();
		VideoSize size1(size3);
		size1.channels = 1;
		size1.update_fields();

		//! Index for a 3 channel video
		unsigned idx3 = size3.index(px, py, pt, 0);

		//! Index for a 1 channel video
		unsigned idx1 = size1.index(px, py, pt);

		//! Conversion formulas
		unsigned idx1_from_idx3 = (idx3 / size3.whc) * size1.wh  + idx3 % size3.wh ;
		unsigned idx3_from_idx1 = (idx1 / size1.wh ) * size3.whc + idx3 % size1.wh ;

		printf("idx3           = %d - idx1           = %d\n", idx3, idx1);
		printf("idx3 from idx1 = %d - idx1 from idx3 = %d\n", idx3_from_idx1, idx1_from_idx3);
	}//*/

}
