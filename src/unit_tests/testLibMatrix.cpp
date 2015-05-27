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
	const unsigned px      = clo_option("-px", 0, "Current point of analysis");
	const unsigned py      = clo_option("-py", 0, "Current point of analysis");
	const unsigned pt      = clo_option("-pt", 0, "Current point of analysis");

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
		Video<float> vid(vid_ori);
		VideoUtils::transformColorSpace(vid, true);

		//! Initialize parameter structures
		VideoNLB::nlbParams prms1, prms2;
		{
			VideoNLB::initializeNlbParameters(prms1, 1, i_sigma, vid.sz, 0, 1);
			VideoNLB::initializeNlbParameters(prms2, 2, i_sigma, vid.sz, 0, 1);
			VideoNLB::printNlbParameters(prms1);
			VideoNLB::printNlbParameters(prms2);
		}

		//! Used matrices during Bayes' estimate
		const unsigned patch_dim = prms1.sizePatch * prms1.sizePatch * prms1.sizePatchTime ;
		const unsigned patch_num = prms1.nSimilarPatches ;

		//! Matrices used for Bayes' estimate
		vector<unsigned> index(patch_num);
		VideoNLB::matWorkspace mat;
		{
			mat.group3dTranspose.resize(patch_num * patch_dim);
			mat.tmpMat          .resize(patch_dim * patch_dim);
			mat.covMat          .resize(patch_dim * patch_dim);
			mat.covMatTmp       .resize(patch_dim * patch_dim);
			mat.baricenter      .resize(patch_dim);
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

		//! Print covariance matrix to a file
		print_matrix(mat.covMat, patch_dim, patch_dim, "covariance_matrix.asc");

		//! Compute first 4 eigenvectors and eigenvalues
		int info = matrixEigs(mat.covMat, patch_dim, 4, mat.covMatTmp, mat.tmpMat);

		print_matrix(mat.covMatTmp, 1, 4,      "eigenvals.asc");
		print_matrix(mat.tmpMat, 4, patch_dim, "eigenvecs.asc");

	}
}
