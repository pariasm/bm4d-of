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
#include <cstdio>
#include <cstdlib>
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

void print_video_size(const std::string& name, const Video<float>& vid)
{
	printf("%s", name.c_str());
	printf("\tFrames: %d - Channels: %d - Height: %d - Width %d\n", 
			vid.sz.frames, vid.sz.channels, vid.sz.height, vid.sz.width);
}

/**
 * @brief For a point in the image, computes the patch group and the
 * corresponding covariance matrix. Then some leading eigenvectors and
 * eigenvalues of this matrix are computed.
 *
 * @return Exit status.
 **/
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
	const unsigned px      = clo_option("-PAx", 0, "Point of analysis");
	const unsigned py      = clo_option("-PAy", 0, "Point of analysis");
	const unsigned pt      = clo_option("-PAt", 0, "Point of analysis");

	//! Video NLB parameters relevant to the computation of distances
	const int time_search  = clo_option("-wt", 2, "> Search window temporal radius");
	const int space_search = clo_option("-wx",-1, "> Search window spatial  radius");
	const int patch_sizex  = clo_option("-px",-1, "> Spatial  patch size");
	const int patch_sizet  = clo_option("-pt", 1, "> Temporal patch size");
	const int num_patches  = clo_option("-np",-1, "> Number of similar patches");
	const int rank         = clo_option("-r" ,-1, "> Covariance matrix rank");

	//! Load video
	Video<float> vid_ori;
	{
		vid_ori.loadVideo(i_video_path.c_str(), i_firstFrame, i_lastFrame, i_frameStep);
		print_video_size("loaded video1", vid_ori);
	}

	//! Test estimateSimilarPatchesStep1 and computeAggregationStep1 on a single pixel
	{
		using std::vector;

		//! RGB to YUV
		Video<float> vid;
		VideoUtils::addNoise(vid_ori, vid, i_sigma, 1);
		VideoUtils::transformColorSpace(vid, true);

		//! Initialize parameter structures
		VideoNLB::nlbParams prms1;
		{
			//! Initialize parameter structures
			VideoNLB::initializeNlbParameters(prms1, 1, i_sigma, vid.sz, 0, 1, time_search, time_search, patch_sizet);

			//! Override with command line parameters
			if (space_search >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search);
			if (patch_sizex  >= 0) VideoNLB::setSizePatch(prms1, vid.sz, (unsigned)patch_sizex);;
			if (num_patches  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches);

			VideoNLB::printNlbParameters(prms1);
		}

		//! Used matrices during Bayes' estimate
		const unsigned patch_dim = prms1.sizePatch * prms1.sizePatch * prms1.sizePatchTime ;
		const unsigned patch_num = prms1.nSimilarPatches ;

		//! Matrices used for Bayes' estimate
		vector<unsigned> index(patch_num);
		VideoNLB::matWorkspace mat;
		{
			mat.groupTranspose.resize(patch_num * patch_dim);
			mat.tmpMat        .resize(patch_dim * patch_dim);
			mat.covMat        .resize(patch_dim * patch_dim);
			mat.covMatTmp     .resize(patch_dim * patch_dim);
			mat.baricenter    .resize(patch_dim);
		}

		//! Matrices used for Bayes' estimate
		vector<unsigned> patch_index(patch_num);
		vector<vector<float> > patch_stack(vid.sz.channels,
		                                   vector<float>(patch_num * patch_dim));

		//! Compute stack of patches
		printf("computing patch group for point: [% 3d,% 3d,% 2d]\n\n",px,py,pt);

		//! Point
		unsigned pind = vid.sz.index(px, py, pt, 0);

		//! Compute distances
		unsigned n_similar = 
			VideoNLB::estimateSimilarPatchesStep1(vid, patch_stack, patch_index, pind, prms1);

		//! Compute covariance matrix (for luminance channel only)
		unsigned channel = 0;
		{
			//! Center data around the baricenter
			centerData(patch_stack[channel], mat.baricenter, n_similar, patch_dim);

			//! Compute the covariance matrix of the set of similar patches
			covarianceMatrix(patch_stack[channel], mat.covMat, n_similar, patch_dim);
		}

		//! Print patch group and covariance matrix to a file
		print_matrix(mat.covMat, patch_dim, patch_dim, "/tmp/covariance_matrix.asc");
		print_matrix(patch_stack[channel], patch_dim, n_similar, "/tmp/data_matrix.asc");

		//! Compute first 4 eigenvectors and eigenvalues
		const int r = rank;
		int info = matrixEigs(mat.covMat, patch_dim, r, mat.covEigVals, mat.covEigVecs);

		printf("matrixEigs exited with status: %d\n", info);

		print_matrix(mat.covEigVals, 1, r,         "/tmp/eigenvals.asc");
		print_matrix(mat.covEigVecs, r, patch_dim, "/tmp/eigenvecs.asc");

		//! Compute SVD with LAPACK
		int info_svd = matrixSVD(patch_stack[channel], patch_dim, n_similar,
		                         mat.svd_S, mat.svd_VT, mat.svd_U,
		                         mat.svd_work, mat.svd_iwork);

		printf("matrixSVD exited with status: %d\n", info);

		int min_dim = std::min(patch_dim, n_similar);
		print_matrix(mat.svd_S , min_dim, 1        , "/tmp/svdS.asc");
		print_matrix(mat.svd_U , min_dim, patch_dim, "/tmp/svdU.asc");
		print_matrix(mat.svd_VT, min_dim, n_similar, "/tmp/svdVT.asc");

//		//! Compute LR SVD with ID
//		{
//			mat.svd_ddata.resize(patch_stack[channel].size());
//			std::vector<double>::iterator ddata = mat.svd_ddata.begin();
//			std::vector<float >::iterator fdata = patch_stack[channel].begin();
//			for (int i = 0; i < patch_stack[channel].size(); ++i)
//				*ddata++ = (double)*fdata++;
//		}
//		const int l = 10;
//		printf("size ddata = %d\n",(int)mat.svd_ddata.size());
//		int info_lrsvd = matrixLRSVD(mat.svd_ddata, patch_dim, n_similar, r + l,
//		                             mat.svd_dS, mat.svd_dV, mat.svd_dU,
//		                             mat.svd_dwork);
//
//		printf("matrixLRSVD exited with status: %d\n", info);
//
//		print_matrix(mat.svd_dS, r + l, 1        , "/tmp/svddS.asc");
//		print_matrix(mat.svd_dU, r + l, patch_dim, "/tmp/svddU.asc");
//		print_matrix(mat.svd_dV, r + l, n_similar, "/tmp/svddV.asc");


	}
}
