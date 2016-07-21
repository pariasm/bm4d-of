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
#include "nldct.h"
#include "cmd_option.h"

#include <algorithm>

/**
 * @file   patch_group.cpp
 * @brief  Retrieves the group of similar of a given one.
 *
 *
 *
 * @author PABLO ARIAS  <pariasm@gmail.com>
 **/

namespace VideoNLB
{
	std::vector<std::pair<float, unsigned> > estimateSimilarPatchesStep1_debug(
		Video<float> const& i_im
	,	std::vector<std::vector<float> > &o_group3d
	,	std::vector<unsigned> &o_index
	,	const unsigned pidx
	,	const nlbParams &p_params
	){
		//! Initialization
		int sWx   = p_params.sizeSearchWindow;
		int sWy   = p_params.sizeSearchWindow;
		const int sWt_f = p_params.sizeSearchTimeRangeFwd;
		const int sWt_b = p_params.sizeSearchTimeRangeBwd;
		const int sPx   = p_params.sizePatch;
		const int sPt   = p_params.sizePatchTime;

		//! Coordinates of center of search box
		unsigned px, py, pt, pc;
		i_im.sz.coords(pidx, px, py, pt, pc);

		unsigned rangex[2];
		unsigned rangey[2];
		unsigned ranget[2];

#ifdef CENTRED_SEARCH
		rangex[0] = std::max(0, (int)px - (sWx-1)/2);
		rangey[0] = std::max(0, (int)py - (sWy-1)/2);
		ranget[0] = std::max(0, (int)pt -  sWt_b   );

		rangex[1] = std::min((int)i_im.sz.width  - sPx, (int)px + (sWx-1)/2);
		rangey[1] = std::min((int)i_im.sz.height - sPx, (int)py + (sWy-1)/2);
		ranget[1] = std::min((int)i_im.sz.frames - sPt, (int)pt +  sWt_f   );
#else
		int shift_x = std::min(0, (int)px - (sWx-1)/2); 
		int shift_y = std::min(0, (int)py - (sWy-1)/2); 
		int shift_t = std::min(0, (int)pt -  sWt_b   ); 

		shift_x += std::max(0, (int)px + (sWx-1)/2 - (int)i_im.sz.width  + sPx); 
		shift_y += std::max(0, (int)py + (sWy-1)/2 - (int)i_im.sz.height + sPx); 
		shift_t += std::max(0, (int)pt +  sWt_f    - (int)i_im.sz.frames + sPt); 

		rangex[0] = std::max(0, (int)px - (sWx-1)/2 - shift_x);
		rangey[0] = std::max(0, (int)py - (sWy-1)/2 - shift_y);
		ranget[0] = std::max(0, (int)pt -  sWt_b    - shift_t);

		rangex[1] = std::min((int)i_im.sz.width  - sPx, (int)px + (sWx-1)/2 - shift_x);
		rangey[1] = std::min((int)i_im.sz.height - sPx, (int)py + (sWy-1)/2 - shift_y);
		ranget[1] = std::min((int)i_im.sz.frames - sPt, (int)pt +  sWt_f    - shift_t);
#endif

		//! Redefine size of search range
		sWx = rangex[1] - rangex[0] + 1;
		sWy = rangey[1] - rangey[0] + 1;
		int sWt = ranget[1] - ranget[0] + 1;

		std::vector<std::pair<float, unsigned> > distance(sWx * sWy * sWt);

		//! Compute distance between patches in search range
		for (unsigned qt = ranget[0], dt = 0; qt <= ranget[1]; qt++, dt++)
		for (unsigned qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
		for (unsigned qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
		{
			//! Squared L2 distance
			float dist = 0.f, dif;
			for (unsigned ht = 0; ht < sPt; ht++)
			for (unsigned hy = 0; hy < sPx; hy++)
			for (unsigned hx = 0; hx < sPx; hx++)
				dist += (dif = i_im(px + hx, py + hy, pt + ht)
								 - i_im(qx + hx, qy + hy, qt + ht)) * dif;

			//! Save distance and corresponding patch index
			distance[dt * sWx*sWy + dy * sWx + dx] = 
				std::make_pair(dist, i_im.sz.index(qx, qy, qt, 0));
		}

		//! Keep only the N2 best similar patches
		const unsigned nSimP = std::min(p_params.nSimilarPatches, (unsigned)distance.size());
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
								distance.end(), comparaisonFirst);

		distance.resize(nSimP);

		if (nSimP <  p_params.nSimilarPatches)
		{
			printf("SR1 [%d,%d,%d] ~ [%d-%d, %d-%d, %d-%d] - nsim = %d\n", 
					px,py,pt,rangex[0], rangex[1], rangey[0], rangey[1], ranget[0], ranget[1], nSimP);
		}

		//! Register position of patches
		for (unsigned n = 0; n < nSimP; n++) o_index[n] = distance[n].second;

		//! Stack selected patches into the 3D group
		const unsigned w   = i_im.sz.width;
		const unsigned wh  = i_im.sz.wh;
		const unsigned whc = i_im.sz.whc;
		for (unsigned c  = 0; c < i_im.sz.channels; c++)
		for (unsigned ht = 0, k = 0; ht < sPt; ht++)
		for (unsigned hy = 0;        hy < sPx; hy++)
		for (unsigned hx = 0;        hx < sPx; hx++)
		for (unsigned n  = 0; n < nSimP; n++, k++)
			o_group3d[c][k] = i_im(c * wh + o_index[n] + ht * whc + hy * w + hx);

		/* 000  pixels from all patches
		 * 001  pixels from all patches
		 * ...
		 * spt,spx,spx pixels from all patches
		 */

		return distance;
	}

	std::vector<std::pair<float, unsigned> > estimateSimilarPatchesStep2_debug(
		Video<float> const& i_imNoisy
	,	Video<float> const& i_imBasic
	,	std::vector<float> &o_group3dNoisy
	,	std::vector<float> &o_group3dBasic
	,	std::vector<unsigned> &o_index
	,	const unsigned pidx
	,	const nlbParams &p_params
	){
		//! Initialization
		int sWx   = p_params.sizeSearchWindow;
		int sWy   = p_params.sizeSearchWindow;
		const int sWt_f = p_params.sizeSearchTimeRangeFwd;
		const int sWt_b = p_params.sizeSearchTimeRangeBwd;
		const int sPx   = p_params.sizePatch;
		const int sPt   = p_params.sizePatchTime;

		//! Coordinates of center of search box
		unsigned px, py, pt, pc;
		i_imBasic.sz.coords(pidx, px, py, pt, pc);

		unsigned rangex[2];
		unsigned rangey[2];
		unsigned ranget[2];

#ifdef CENTRED_SEARCH
		rangex[0] = std::max(0, (int)px - (sWx-1)/2);
		rangey[0] = std::max(0, (int)py - (sWy-1)/2);
		ranget[0] = std::max(0, (int)pt -  sWt_b   );

		rangex[1] = std::min((int)i_imNoisy.sz.width  - sPx, (int)px + (sWx-1)/2);
		rangey[1] = std::min((int)i_imNoisy.sz.height - sPx, (int)py + (sWy-1)/2);
		ranget[1] = std::min((int)i_imNoisy.sz.frames - sPt, (int)pt +  sWt_f   );
#else
		int shift_x = std::min(0, (int)px - (sWx-1)/2); 
		int shift_y = std::min(0, (int)py - (sWy-1)/2); 
		int shift_t = std::min(0, (int)pt -  sWt_b   ); 

		shift_x += std::max(0, (int)px + (sWx-1)/2 - (int)i_imNoisy.sz.width  + sPx); 
		shift_y += std::max(0, (int)py + (sWy-1)/2 - (int)i_imNoisy.sz.height + sPx); 
		shift_t += std::max(0, (int)pt +  sWt_f    - (int)i_imNoisy.sz.frames + sPt); 

		rangex[0] = std::max(0, (int)px - (sWx-1)/2 - shift_x);
		rangey[0] = std::max(0, (int)py - (sWy-1)/2 - shift_y);
		ranget[0] = std::max(0, (int)pt -  sWt_b    - shift_t);

		rangex[1] = std::min((int)i_imNoisy.sz.width  - sPx, (int)px + (sWx-1)/2 - shift_x);
		rangey[1] = std::min((int)i_imNoisy.sz.height - sPx, (int)py + (sWy-1)/2 - shift_y);
		ranget[1] = std::min((int)i_imNoisy.sz.frames - sPt, (int)pt +  sWt_f    - shift_t);
#endif

		//! Redefine size of search range
		sWx = rangex[1] - rangex[0] + 1;
		sWy = rangey[1] - rangey[0] + 1;
		int sWt = ranget[1] - ranget[0] + 1;

		std::vector<std::pair<float, unsigned> > distance(sWx * sWy * sWt);

		//! Compute distance between patches in search range
		const int chnls = i_imNoisy.sz.channels;
		for (unsigned qt = ranget[0], dt = 0; qt <= ranget[1]; qt++, dt++)
		for (unsigned qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
		for (unsigned qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
		{
			//! Squared L2 distance between color patches of basic image
			float dist = 0.f, dif;
			for (unsigned c = 0; c < chnls; c++)
			for (unsigned ht = 0; ht < sPt; ht++)
			for (unsigned hy = 0; hy < sPx; hy++)
			for (unsigned hx = 0; hx < sPx; hx++)
				dist += (dif = i_imBasic(px + hx, py + hy, pt + ht, c)
								 - i_imBasic(qx + hx, qy + hy, qt + ht, c) ) * dif;

			//! Save distance and corresponding patch index
			distance[dt * sWx*sWy + dy * sWx + dx] = 
				std::make_pair(dist, i_imBasic.sz.index(qx, qy, qt, 0));
		}

		//! Keep only the nSimilarPatches best similar patches
		unsigned nSimP = std::min(p_params.nSimilarPatches, (unsigned)distance.size());
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
								distance.end(), comparaisonFirst);

		distance.resize(nSimP);

		if (nSimP <  p_params.nSimilarPatches)
		{
			printf("SR2 [%d,%d,%d] ~ [%d-%d, %d-%d, %d-%d] - nsim = %d\n", 
					px,py,pt,rangex[0], rangex[1], rangey[0], rangey[1], ranget[0], ranget[1], nSimP);
		}

		//! Save index of similar patches
		const float threshold = (p_params.tau > distance[nSimP - 1].first ?
										 p_params.tau : distance[nSimP - 1].first);
		nSimP = 0;

		//! Register position of similar patches
		for (unsigned n = 0; n < distance.size(); n++)
			if (distance[n].first < threshold)
				o_index[nSimP++] = distance[n].second;

		//! Save similar patches into 3D groups
		const unsigned w   = i_imNoisy.sz.width;
		const unsigned wh  = i_imNoisy.sz.wh;
		const unsigned whc = i_imNoisy.sz.whc;
		for (unsigned c = 0, k = 0; c < chnls; c++)
		for (unsigned ht = 0; ht < sPt; ht++)
		for (unsigned hy = 0; hy < sPx; hy++)
		for (unsigned hx = 0; hx < sPx; hx++)
		for (unsigned n = 0; n < nSimP; n++, k++)
		{
			o_group3dNoisy[k] = i_imNoisy(c * wh + o_index[n] + ht * whc + hy * w + hx);
			o_group3dBasic[k] = i_imBasic(c * wh + o_index[n] + ht * whc + hy * w + hx);
		}

		return distance;
	}
}

void print_group_coordinates(
 	std::vector<std::pair<float, unsigned> > &patch_dists
,	VideoSize sz
,	std::string out
){
	FILE *fid = out.size() ? fopen(out.c_str(),"a") : stdout;
	for (unsigned idx = 0; idx < patch_dists.size(); idx++) 
	{
		unsigned px, py, pt, pc;
		sz.coords(patch_dists[idx].second, px, py, pt, pc);
		fprintf(fid, "[COORD]% 3d % 3d % 3d\n", px, py, pt);
		fprintf(fid, "[DISTA]%f\n", patch_dists[idx].first);
	}
	if (out.size()) fclose(fid);

}

void print_group(
	std::vector<float> &group3d
,	const unsigned patch_dim
,	const unsigned patch_num
,	VideoSize sz
,	std::string prefix
,	std::string out
){
	FILE *fid = out.size() ? fopen(out.c_str(),"a") : stdout;
	for (unsigned idx = 0; idx < patch_num; idx++) 
	{
		fprintf(fid, "%s", prefix.c_str());
		for (unsigned k = 0; k < patch_dim; k++)
			fprintf(fid, "% 3f   ", group3d[k*patch_num + idx]);
		fprintf(fid, "\n");
	}
	if (out.size()) fclose(fid);

}



void print_video_size(const std::string& name, const Video<float>& vid)
{
	printf("%s", name.c_str());
	printf("\tFrames: %d - Channels: %d - Height: %d - Width %d\n", 
			vid.sz.frames, vid.sz.channels, vid.sz.height, vid.sz.width);
}

int main(int argc, char **argv)
{
	clo_usage("Retrieves group of similar patches and outputs it in ASCII files");

	//! Get command line inputs
	const string i_noisy_path = clo_option("-i", "", "Patch to noisy sequence (printf format)");
	const string i_basic_path = clo_option("-b", "", "Patch to basic sequence (printf format)");

	const string output_file  = clo_option("-out", "stdout", "Output file");

	const int i_firstFrame = clo_option("-f", 0, "First frame");
	const int i_lastFrame  = clo_option("-l", 0, "Last frame");
	const int i_frameStep  = clo_option("-s", 1, "Frame step");

	const int sigma        = clo_option("-sigma", 0, "Standard deviation of the noise");
	const unsigned px      = clo_option("-PAx", 0, "Point of analysis");
	const unsigned py      = clo_option("-PAy", 0, "Point of analysis");
	const unsigned pt      = clo_option("-PAt", 0, "Point of analysis");

	//! Video NLB parameters relevant to the computation of distances
	const int time_search  = clo_option("-wt", 2, "> Search window temporal radius, step 1");
	const int space_search = clo_option("-wx",-1, "> Search window spatial radius, step 1");
	const int patch_sizex  = clo_option("-px",-1, "> Spatial patch size, step 1");
	const int patch_sizet  = clo_option("-pt", 1, "> Temporal patch size, step 1");
	const int num_patches  = clo_option("-np",-1, "> Number of similar patches, step 1");

	const int step = (i_basic_path.size()) ? 2 : 1;


	//! Load video
	Video<float> vid_nsy;
	Video<float> vid_bsc;
	{
		vid_nsy.loadVideo(i_noisy_path.c_str(), i_firstFrame, i_lastFrame, i_frameStep);
//		print_video_size("loaded noisy video: ", vid_nsy);
		
		if (step == 2)
		{
			vid_bsc.loadVideo(i_basic_path.c_str(), i_firstFrame, i_lastFrame, i_frameStep);
//			print_video_size("loaded basic video: ", vid_bsc);
		}
	}

	//! Set parameters
	VideoNLB::nlbParams prms;
	{
		//! Initialize parameter structures
		VideoNLB::initializeNlbParameters(prms, step, sigma, vid_nsy.sz, 0, 1, time_search, time_search, patch_sizet);

		//! Override with command line parameters
		if (space_search >= 0) VideoNLB::setSizeSearchWindow(prms, (unsigned)space_search);
		if (patch_sizex  >= 0) VideoNLB::setSizePatch(prms, vid_nsy.sz, (unsigned)patch_sizex);;
		if (num_patches  >= 0) VideoNLB::setNSimilarPatches(prms, (unsigned)num_patches);

		//! Print parameters
//		VideoNLB::printNlbParameters(prms);
	}

	//! Compute group of similar patches on selected pixel single pixel
	{
		using std::vector;

//		printf("computing patch distances for point: [% 3d,% 3d,% 2d]\n\n",px,py,pt);

		//! Point
		unsigned pind = vid_nsy.sz.index(px, py, pt, 0);

		//! Dimensions of patch groups
		unsigned patch_dim;
		unsigned patch_num;
		{
			const unsigned sWx = prms.sizeSearchWindow;
			const unsigned sWt = prms.sizeSearchTimeRangeFwd +
										prms.sizeSearchTimeRangeBwd + 1;
			const unsigned sPx = prms.sizePatch;
			const unsigned sPt = prms.sizePatchTime;

			patch_dim = (step == 1) ? sPx * sPx * sPt : sPx * sPx * sPt * vid_nsy.sz.channels;
			patch_num = (step == 1) ? prms.nSimilarPatches : sWx * sWx * sWt;
		}

		//! Compute distances
		vector<std::pair<float, unsigned> > patch_dists;
		vector<unsigned> patch_index(patch_num);

		if (step == 1)
		{
			//! RGB to YUV
			Video<float> vid(vid_nsy);
			VideoUtils::transformColorSpace(vid, true);

			vector<vector<float> > patch_stack_nsy(vid.sz.channels, 
			                                       vector<float>(patch_num * patch_dim));

			patch_dists = VideoNLB::estimateSimilarPatchesStep1_debug(vid,
					patch_stack_nsy, patch_index, pind, prms);

//			for (unsigned c = 0; c < vid.sz.channels; c++)
//			{
//				char prefix[1024];
//				sprintf(prefix, "[YUV%d]", c);
//				print_group(patch_stack_nsy[c], patch_dim,
//						patch_dists.size(), vid_nsy.sz, prefix, output_file);
//			}
		}
		else
		{
			vector<float> patch_stack_nsy(patch_num * patch_dim);
			vector<float> patch_stack_bsc(patch_num * patch_dim);

			patch_dists = VideoNLB::estimateSimilarPatchesStep2_debug(vid_nsy,
					vid_bsc, patch_stack_nsy, patch_stack_bsc, patch_index, pind,
					prms);

//			char prefix[1024];
//			sprintf(prefix, "[BSIC]");
//			print_group(patch_stack_bsc, patch_dim,
//					patch_dists.size(), vid_nsy.sz, prefix, output_file);
//			sprintf(prefix, "[NISY]");
//			print_group(patch_stack_nsy, patch_dim,
//					patch_dists.size(), vid_nsy.sz, prefix, output_file);
		}

		print_group_coordinates(patch_dists, vid_nsy.sz, output_file);

		/*! Compute stack of patches (detailed output)
		{
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
		}//*/

	}

}
