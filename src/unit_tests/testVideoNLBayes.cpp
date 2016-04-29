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
#include "LibMatrix.h"
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

// Assumes row-major ordering
void print_matrix(
	std::vector<float> &matrix
,	unsigned rows
,	unsigned cols
,	std::string filename
){
	FILE *file = fopen(filename.c_str(),"w");

	// print output
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			fprintf(file, "%.16g ", matrix[i*cols + j]);
		}
		fprintf(file, "\n");
	}

	fclose(file);
}

// Assumes row-major ordering
void print_matrix(
	std::vector<double> &matrix
,	unsigned rows
,	unsigned cols
,	std::string filename
){
	FILE *file = fopen(filename.c_str(),"w");

	// print output
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			fprintf(file, "%.16g ", matrix[i*cols + j]);
		}
		fprintf(file, "\n");
	}

	fclose(file);
}

namespace VideoNLB
{
	std::vector<std::pair<float, unsigned> > estimateSimilarPatchesStep1_debug(
		Video<float> const& i_im
	,	std::vector<std::vector<float> > &o_group
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
			o_group[c][k] = i_im.data[o_index[n] + hy * w + hx + c * wh];

		return distance;
	
		/* 00  pixels from all patches
		 * 01  pixels from all patches
		 * ...
		 * 0sp pixels from all patches
		 */
	}

	/**
	 * @brief Debug version of the version of
	 * VideoNLB::computeBayesEstimateStep2_LRSVD that uses iddist to compute
	 * a truncated SVD.
	 *
	 * @param io_groupNoisy: inputs all similar patches in the noisy image,
	 *                         outputs their denoised estimates.
	 * @param i_groupBasic: contains all similar patches in the basic image.
	 * @param i_mat: contains :
	 *    - groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
	 *    - baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
	 *    - covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
	 *    - covMatTmp: allocated memory. Used to process the Bayes estimate;
	 *    - tmpMat: allocated memory. Used to process the Bayes estimate;
	 * @param io_nInverseFailed: update the number of failed matrix inversion;
	 * @param p_imSize: size of the image;
	 * @param p_params: see processStep2 for more explanations;
	 * @param p_nSimP: number of similar patches.
	 *
	 * @return none.
	 **/
	float computeBayesEstimateStep2_LRSVD_id(
		std::vector<float> &io_groupNoisy
	,	std::vector<float>  &i_groupBasic
	,	matWorkspace &i_mat
	,	unsigned &io_nInverseFailed
	,	const VideoSize &p_imSize
	,	nlbParams const& p_params
	,	const unsigned p_nSimP
	){
		//! Parameters initialization
		const float diagVal = p_params.beta * p_params.sigma * p_params.sigma;
		const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
		                    * p_params.sizePatchTime * p_imSize.channels;
		const unsigned r    = p_params.rank;
	
		//! Center 3D groups around their baricenter
		centerData( i_groupBasic, i_mat.baricenter, p_nSimP, sPC);
		centerData(io_groupNoisy, i_mat.baricenter, p_nSimP, sPC);
	
		//! Compute total variance HOW?
		float total_variance = 0.f;
	
		//! Compute SVD
		{
			// convert data to double
			i_mat.svd_ddata.resize(i_groupBasic.size());
			std::vector<double>::iterator ddata = i_mat.svd_ddata.begin();
			std::vector<float >::iterator fdata =  i_groupBasic.begin();
			for (int i = 0; i < i_groupBasic.size(); ++i)
				*ddata++ = (double)*fdata++;

			print_matrix(i_mat.svd_ddata, sPC, p_nSimP, "/tmp/Xt.asc");
	
			// compute low rand SVD
			int info = matrixLRSVD(i_mat.svd_ddata, p_nSimP, sPC, 4*r,
					i_mat.svd_dS, i_mat.svd_dV, i_mat.svd_dU,
					i_mat.svd_dwork);

			print_matrix(i_mat.svd_dU , sPC, r, "/tmp/Ut.asc");
			print_matrix(i_mat.svd_dS , 1  , r, "/tmp/S.asc");
			print_matrix(i_mat.svd_dV , p_nSimP, r, "/tmp/Vt.asc");
	
			// convert SVD matrices to float
			i_mat.svd_S.resize(i_mat.svd_dS.size());
			i_mat.svd_V.resize(i_mat.svd_dV.size());
			i_mat.svd_U.resize(i_mat.svd_dU.size());
	
			std::vector<float >::iterator to;
			std::vector<double>::iterator from;
	
			to   = i_mat.svd_S .begin();
			from = i_mat.svd_dS.begin();
			for (int i = 0; i < i_mat.svd_dS.size(); ++i) *to++ = *from++;
	
			to   = i_mat.svd_V .begin();
			from = i_mat.svd_dV.begin();
			for (int i = 0; i < i_mat.svd_dV.size(); ++i) *to++ = *from++;
	
			to   = i_mat.svd_U .begin();
			from = i_mat.svd_dU.begin();
			for (int i = 0; i < i_mat.svd_dU.size(); ++i) *to++ = *from++;
		}
	
		//! Compute variance captured by the r leading eigenvectors
		float r_variance = 0.f;
		for (int i = 0; i < r; ++i)
		{
			i_mat.svd_S[i] = (i_mat.svd_S[i] * i_mat.svd_S[i]) / (float)p_nSimP;
			r_variance += i_mat.svd_S[i];
		}
	
		//! Compute eigenvalues-based coefficients of Bayes' filter
		float sigma2 = p_params.sigma * p_params.sigma;
		for (unsigned k = 0; k < r; ++k)
			i_mat.svd_S[k] = 1.f / sqrtf( 1.f + sigma2 / i_mat.svd_S[k] );
	
		//! Multiply sing. vector matrix by singular value matrix
		{
			float *svdU = i_mat.svd_U.data();
			for (unsigned k = 0; k < r  ; ++k)
			for (unsigned i = 0; i < sPC; ++i)
				*svdU++ *= i_mat.svd_S[k];
		}
	
		//! Z' = X'*V
		productMatrix(i_mat.groupTranspose,
		              io_groupNoisy,
		              i_mat.svd_U,
		              p_nSimP, r, sPC,
		              false, false);
	
		//! hX' = Z'*V'
		productMatrix(io_groupNoisy,
		              i_mat.groupTranspose,
		              i_mat.svd_U,
		              p_nSimP, sPC, r,
		              false, true);
	
		//! Add baricenter
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				io_groupNoisy[k] += i_mat.baricenter[j];
	
		// return percentage of captured variance
		return r_variance / total_variance;
	
	}
}

void print_patch_group(
	std::vector<std::vector<float> > &group
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
			printf("% 3f,% 3f,% 3f   ", group[0][k*nSimP + 0], 
			                            group[1][k*nSimP + 0],
			                            group[2][k*nSimP + 0]);
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
			printf("% 3f,% 3f,% 3f   ", group[0][k*nSimP + idx], 
			                            group[1][k*nSimP + idx],
			                            group[2][k*nSimP + idx]);
		}
		printf("\n");
	}
	printf("\ndistance between them: %f\n\n\n", patch_dists[idx].first);
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
	const string i_noisy_path = clo_option("-i", "", "Patch to noisy sequence (printf format)");
	const string i_basic_path = clo_option("-b", "", "Patch to basic sequence (printf format)");

	const string output_file  = clo_option("-out", "stdout", "Output file");

	const int i_firstFrame = clo_option("-f", 0, "First frame");
	const int i_lastFrame  = clo_option("-l", 0, "Last frame");
	const int i_frameStep  = clo_option("-s", 1, "Frame step");

	const int i_sigma      = clo_option("-sigma", 0, "Add noise of standard deviation sigma");
	const unsigned px      = clo_option("-PAx", 100, "Point of analysis");
	const unsigned py      = clo_option("-PAy", 100, "Point of analysis");
	const unsigned pt      = clo_option("-PAt",   2, "Point of analysis");

	//! Video NLB parameters
	const int time_search1  = clo_option("-wt1", 2, "> Search window temporal radius, step 1");
	const int time_search2  = clo_option("-wt2", 2, "> Search window temporal radius, step 2");
	const int space_search1 = clo_option("-wx1",-1, "> Search window spatial radius, step 1");
	const int space_search2 = clo_option("-wx2",-1, "> Search window spatial radius, step 2");
	const int patch_sizex1  = clo_option("-px1",-1, "> Spatial patch size, step 1");
	const int patch_sizex2  = clo_option("-px2",-1, "> Spatial patch size, step 2");
	const int patch_sizet1  = clo_option("-pt1", 1, "> Temporal patch size, step 1");
	const int patch_sizet2  = clo_option("-pt2", 1, "> Temporal patch size, step 2");
	const int num_patches1  = clo_option("-np1",-1, "> Number of similar patches, step 1");
	const int num_patches2  = clo_option("-np2",-1, "> Number of similar patches, step 2");
	const int rank1         = clo_option("-r1" , 4, "> Rank or covariance matrix, step 1");
	const int rank2         = clo_option("-r2" , 4, "> Rank or covariance matrix, step 2");

	const int step = (i_basic_path.size()) ? 2 : 1;

	//! Load video
	Video<float> vid_nsy;
	Video<float> vid_bsc;
	{
		vid_nsy.loadVideo(i_noisy_path.c_str(), i_firstFrame, i_lastFrame, i_frameStep);

		if (step == 2)
		{
			vid_bsc.loadVideo(i_basic_path.c_str(), i_firstFrame, i_lastFrame, i_frameStep);
			print_video_size("loaded basic video: ", vid_bsc);
		}
		else
		{
			Video<float> vid_ori(vid_nsy);
			VideoUtils::addNoise(vid_ori, vid_nsy, i_sigma, 1);
		}
		
	}

	//! Set parameters
	VideoNLB::nlbParams prms1, prms2;
	{
		VideoNLB::initializeNlbParameters(prms1, 1, i_sigma, vid_nsy.sz, 0, 1, time_search1, time_search1, patch_sizet1);
		VideoNLB::initializeNlbParameters(prms2, 2, i_sigma, vid_nsy.sz, 0, 1, time_search2, time_search2, patch_sizet2);

		//! Override with command line parameters
		if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
		if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
		if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, vid_nsy.sz, (unsigned)patch_sizex1);;
		if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, vid_nsy.sz, (unsigned)patch_sizex2);;
		if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
		if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);

		prms1.rank = rank1;
		prms2.rank = rank2;
	}

	VideoNLB::printNlbParameters(prms1);
	VideoNLB::printNlbParameters(prms2);


	/*/! Test estimateSimilarPatchesStep1 and computeAggregationStep1 on a single pixel
	{
		using std::vector;

		//! RGB to YUV
		Video<float> vid(vid_nsy);
		VideoUtils::transformColorSpace(vid, true);

		//! Initialize parameter structures
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, 1, i_sigma, vid.sz, 0, 1);
		VideoNLB::initializeNlbParameters(prms2, 2, i_sigma, vid.sz, 0, 1);
		VideoNLB::printNlbParameters(prms1);
		VideoNLB::printNlbParameters(prms2);

		//! Used matrices during Bayes' estimate
		const unsigned patch_dim = prms1.sizePatch * prms1.sizePatch ;
		const unsigned patch_num = prms1.nSimilarPatches ;

		//! Matrices used for Bayes' estimate
		vector<unsigned> patch_index(patch_num);
		vector<vector<float> > patch_stack(vid.sz.channels,
				                             vector<float>(patch_num * patch_dim));

		//! Compute stack of patches (detailed output)
		{
			printf("computing patch distances for point: [% 3d,% 3d,% 2d]\n\n",px,py,pt);

			//! Point
			unsigned pind = vid.sz.index(px, py, pt, 0);

			//! Compute distances
			vector<std::pair<float, unsigned> > patch_dists = 
				VideoNLB::estimateSimilarPatchesStep1_debug(vid, patch_stack, patch_index, pind, prms1);

			//! Display seleted patches
			for (int n = 0; n < prms1.nSimilarPatches; n++)
				print_patch_group(patch_stack, patch_dists, prms1, vid.sz, n);

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

		//! Aggregate 
		Video<float> basic(vid.sz);
		Video<float> weight(vid.sz, 0.f);
		std::vector<bool> mask(vid.sz.whf, true);
		VideoNLB::computeAggregationStep1(basic, weight, mask, patch_stack, patch_index, prms1);

		//! Save outputs to disk
		basic .saveVideoAscii("/tmp/basic" , i_firstFrame, i_frameStep);
		weight.saveVideoAscii("/tmp/weight", i_firstFrame, i_frameStep);
		
		Video<float> mask_out(vid.sz.width, vid.sz.height, vid.sz.frames);
		for (unsigned i = 0; i < mask.size(); i++)
			mask_out(i) = 255.*(float)mask[i];
		mask_out.saveVideo("/tmp/mask_%02d.png", i_firstFrame, i_frameStep);

	}//*/

	/*/! Test estimateSimilarPatchesStep1 and computeAggregationStep1 on all video
	{
		using std::vector;

		//! RGB to YUV
		Video<float> vid(vid_nsy);
		VideoUtils::transformColorSpace(vid, true);

		//! Initialize parameter structures
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, 1, i_sigma, vid.sz, 0, 1);
		VideoNLB::initializeNlbParameters(prms2, 2, i_sigma, vid.sz, 0, 1);
		VideoNLB::printNlbParameters(prms1);
		VideoNLB::printNlbParameters(prms2);

		//! Used matrices during Bayes' estimate
		const unsigned patch_dim = prms1.sizePatch * prms1.sizePatch ;
		const unsigned patch_num = prms1.nSimilarPatches ;

		//! Matrices used for Bayes' estimate
		vector<unsigned> patch_index(patch_num);
		vector<vector<float> > patch_stack(vid.sz.channels,
				                             vector<float>(patch_num * patch_dim));

		Video<float> basic(vid.sz);
		Video<float> weight(vid.sz, 0.f);
		std::vector<bool> mask(vid.sz.whf, false);

		//! Only pixels of the center of the image must be processed (not the boundaries)
		unsigned sWx = prms1.sizeSearchWindow;
		for (unsigned f =   0; f < vid.sz.frames      ; f++)
		for (unsigned i = sWx; i < vid.sz.height - sWx; i++)
		for (unsigned j = sWx; j < vid.sz.width  - sWx; j++)
			mask[f * vid.sz.wh + i * vid.sz.width + j] = true;


		//! Aggregate
		for (int pt = 0; pt < vid.sz.frames; pt++)
		for (int py = 0; py < vid.sz.height; py++)
		for (int px = 0; px < vid.sz.width ; px++)
		if  (mask[pt * vid.sz.wh + py * vid.sz.width + px])
		{
			//! Point
			unsigned pind = vid.sz.index(px, py, pt, 0);

			//! Compute stack of patches
			VideoNLB::estimateSimilarPatchesStep1(vid, patch_stack, patch_index, pind, prms1);

			//! Aggregate 
			VideoNLB::computeAggregationStep1(basic, weight, mask, patch_stack, patch_index, prms1);
		}

		//! Normalize
		for (int pt = 0; pt < vid.sz.frames  ; pt++)
		for (int pc = 0; pc < vid.sz.channels; pc++)
		for (int py = 0; py < vid.sz.height  ; py++)
		for (int px = 0; px < vid.sz.width   ; px++)
			if (weight(px,py,pt,pc))
				basic(px,py,pt,pc) /= weight(px,py,pt,pc);

		//! YUV to RGB
		VideoUtils::transformColorSpace(basic, false);

		//! Save outputs to disk
		basic .saveVideo("/tmp/basic_%02d.png" , i_firstFrame, i_frameStep);
		weight.saveVideoAscii("/tmp/weight", i_firstFrame, i_frameStep);
		
		Video<float> mask_out(vid.sz.width, vid.sz.height, vid.sz.frames);
		for (unsigned i = 0; i < mask.size(); i++)
			mask_out(i) = 255.*(float)mask[i];
		mask_out.saveVideo("/tmp/mask_%02d.png", i_firstFrame, i_frameStep);

	}//*/

	//! Compute full basic estimate, and test estimateSimilarPatchesStep2 on a single pixel
	{
		using std::vector;

		//! RGB to YUV
		Video<float> vid(vid_nsy);

		//! Run first step
		Video<float> basic(vid.sz);
		{
			//! Number of subvideos in which video is divided
			const unsigned nParts = 2;

			if (prms1.verbose)
			{
				printf("1st Step\n");
				for (int p = 0; p < nParts; ++p) printf("\n");
			}

			//! RGB to YUV
			VideoUtils::transformColorSpace(vid, true);

			//! Divide in subvideos
			std::vector<Video<float> > parts_vid(nParts);
			std::vector<VideoUtils::CropPosition > crops(nParts);
			VideoUtils::subDivideTight(vid, parts_vid, crops, prms1.boundary,
					nParts);

			//! Process all sub-videos
			std::vector<Video<float> > parts_basic(nParts);
			std::vector<Video<float> > parts_final(nParts);
			for (int n = 0; n < (int)nParts; n++)
				processNlBayes(parts_vid[n], parts_basic[n], parts_final[n],
						prms1, crops[n]);

			//! Get the basic estimate
			VideoUtils::subBuildTight(parts_basic, basic, prms1.boundary);

			//! YUV to RGB
			VideoUtils::transformColorSpace(vid  , false);
			VideoUtils::transformColorSpace(basic, false);

			//! Save video
			basic.saveVideo("/tmp/basic_%02d.png", i_firstFrame, i_frameStep);
		}

		//! Compute second stage Bayes estimate for the point of analysis 
		{

			printf("2n step for point p = [%d,%d,%d]\n",px, py, pt);

			//! Allocate
			Video<float> final(vid.sz);

			const unsigned sWx = prms2.sizeSearchWindow;
			const unsigned sWt = prms2.sizeSearchTimeRangeFwd +
										prms2.sizeSearchTimeRangeBwd + 1;
			const unsigned sPx = prms2.sizePatch;
			const unsigned sPt = prms2.sizePatchTime;
			const VideoSize sz = vid.sz;

			//! Matrices used during Bayes' estimate
			const unsigned patch_dim = sPx * sPx * sPt * sz.channels;
			const unsigned patch_num = sWx * sWx * sWt;

			//! Matrices used for Bayes' estimate
			VideoNLB::matWorkspace mat;
			mat.groupTranspose.resize(patch_num * patch_dim);
			mat.tmpMat        .resize(patch_dim * patch_dim);
			mat.covMat        .resize(patch_dim * patch_dim);
			mat.covMatTmp     .resize(patch_dim * patch_dim);
			mat.baricenter    .resize(patch_dim);
	
			//! Matrices used for Bayes' estimate
			vector<unsigned> index(patch_num);
			vector<float> group_nsy(patch_num * patch_dim);
			vector<float> group_bsc(patch_num * patch_dim);
			vector<float> agg_weights(patch_num);
	
			//! Indices of the point of analysis in a color and scalar video
			const unsigned ij  = sz.index(px,py,pt);
			const unsigned ij3 = sz.index(px,py,pt, 0);

			//! Search for similar patches around the reference one
			unsigned nSimP = estimateSimilarPatchesStep2(vid_nsy, basic,
					group_nsy, group_bsc, index, ij3, prms2);

			printf("\tFound %d similar patches at p = [%d,%d,%d]\n", nSimP, px, py, pt);

			//! Bayes' estimate
			unsigned failed_inv = 0;
			const float var_id = computeBayesEstimateStep2_LRSVD_id(group_nsy, group_bsc, mat, failed_inv, sz, prms2, nSimP);
			printf("\tid done\n");
			const float var_fl = computeBayesEstimateStep2_FR      (group_nsy, group_bsc, mat, failed_inv, sz, prms2, nSimP, agg_weights);
			printf("\tfull rank done\n");

			printf("\tVariances of group at p = [%d,%d,%d]\n",px, py, pt);
			printf("\t\tw/id      var = %f\n",var_id);
			printf("\t\tw/full    var = %f\n",var_fl);
		}

	}//*/

	/*/! Test processNlBayes, both steps
	{
		//! Add noise
		Video<float> vid;
		VideoUtils::addNoise(vid_nsy, vid, i_sigma, true);

		//! Save noisy video
		vid.saveVideo("/tmp/vid_%02d.png", i_firstFrame, i_frameStep);

		//! Initialize parameter structures
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, 1, i_sigma, vid.sz, 0, 1);
		VideoNLB::initializeNlbParameters(prms2, 2, i_sigma, vid.sz, 0, 1);
		VideoNLB::printNlbParameters(prms1);
		VideoNLB::printNlbParameters(prms2);

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
				processNlBayes(parts_vid[n], parts_basic[n], parts_final[n], prms1);

			//! Get the basic estimate
			VideoUtils::subBuild(parts_basic, basic, prms1.boundary);

			//! YUV to RGB
			VideoUtils::transformColorSpace(vid  , false);
			VideoUtils::transformColorSpace(basic, false);
		}

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
		VideoSize size3 = vid_nsy.sz;
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
