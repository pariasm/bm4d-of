/*
 * Original work Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Modified work Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
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
 * @file VideoNLBayes.cpp
 * @brief Video NL-Bayes denoising functions
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 * @author Pablo Arias <pariasm@gmail.com>
 **/

#include <stdexcept>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>

#include "VideoNLBayes.h"
#include "LibMatrix.h"
#include "../Utilities/Utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

//#define DEBUG_SHOW_PATCH_GROUPS
#define DEBUG_SHOW_WEIGHT
#define CENTRED_SEARCH

namespace VideoNLB
{

/**
 * @brief Initialize Parameters of the NL-Bayes algorithm.
 *
 * @param o_params   : will contain the nlbParams for the first step of the algorithm;
 * @param p_step     : select first or second step;
 * @param p_sigma    : standard deviation of the noise;
 * @param p_size     : size of the video;
 * @param p_flatArea : if true, use the homogeneous area trick for the first step;
 * @param p_verbose  : if true, print some informations.
 * @param p_timeSearchRagneFwd : temporal search range forwards.
 * @param p_timeSearchRagneBwd : temporal search range backwards.
 *
 * @return none.
 **/
void initializeNlbParameters(
	nlbParams &o_params
,	const unsigned p_step
,	const float p_sigma
,	const VideoSize &p_size
,	const bool p_flatArea
,	const bool p_verbose
,	const unsigned timeSearchRangeFwd
,	const unsigned timeSearchRangeBwd
,	const unsigned sizePatchTime
){
	const bool s1 = (p_step == 1);

	//! Standard deviation of the noise
	o_params.sigma = p_sigma;

	//! Patch size, spatial dimension
	if (p_size.channels == 1)
	{
		if(s1) o_params.sizePatch = (p_sigma < 30.f ? 5 : 7);
		else   o_params.sizePatch = 5;
	}
	else
	{
		if(s1) o_params.sizePatch = (p_sigma < 20.f ? 3 : (p_sigma < 50.f ? 5 : 7));
		else   o_params.sizePatch = (p_sigma < 50.f ? 3 : (p_sigma < 70.f ? 5 : 7));
	}

	//! Patch size, temporal dimension
	o_params.sizePatchTime = sizePatchTime;

	//! Number of similar patches
	if (p_size.channels == 1)
	{
		if(s1) o_params.nSimilarPatches = (p_sigma < 10.f ?  35 : 
		                                  (p_sigma < 30.f ?  45 : 
		                                  (p_sigma < 80.f ?  90 : 
		                                                    100))); // ASK MARC differs from manuscript
		else   o_params.nSimilarPatches = (p_sigma < 20.f ?  15 : 
		                                  (p_sigma < 40.f ?  25 :
		                                  (p_sigma < 80.f ?  30 :
		                                                     45))); // ASK MARC differs from manuscript
	}
	else
		o_params.nSimilarPatches = o_params.sizePatch * o_params.sizePatch * 3;

	//! Offset: step between two similar patches
	o_params.offSet = o_params.sizePatch / 2;

	//! Use the homogeneous area detection trick
	o_params.useHomogeneousArea = p_flatArea;

	//! Size of the search window around the reference patch (must be odd)
	o_params.sizeSearchWindow = o_params.nSimilarPatches / 2; // ASK MARC and 7*sizePatch1 in IPOL manuscript
	if (o_params.sizeSearchWindow % 2 == 0) o_params.sizeSearchWindow++;

	//! Search window, temporal search radii
	o_params.sizeSearchTimeRangeFwd = timeSearchRangeFwd;
	o_params.sizeSearchTimeRangeBwd = timeSearchRangeBwd;
	o_params.nSimilarPatches *= timeSearchRangeFwd + timeSearchRangeBwd + 1; // FIXME: this is just a test

	//! Size of boundaries used during the sub division
	o_params.boundary = 2*(o_params.sizeSearchWindow/2) + (o_params.sizePatch - 1);

	//! Parameter used to determine if an area is homogeneous
	o_params.gamma = 1.05f;

	//! Parameter used to estimate the covariance matrix
	if (p_size.channels == 1)
	{
		if(s1) o_params.beta = (p_sigma < 15.f ? 1.1f : (p_sigma < 70.f ? 1.f : 0.9f));
		else   o_params.beta = (p_sigma < 15.f ? 1.1f : (p_sigma < 35.f ? 1.f : 0.9f));
	}
	else
	{
		if(s1) o_params.beta = 1.f;
		else   o_params.beta = (p_sigma < 50.f ? 1.2f : 1.f);
	}

	//! Parameter used to determine similar patches
	//  Differs from manuscript (tau = 4) because (1) this threshold is to be used 
	//  with the squared distance and (2) the distance is not normalized
	if(s1) o_params.tau = 0; // not used
	else   o_params.tau = 16.f * o_params.sizePatch * o_params.sizePatch * p_size.channels;

	//! Print information?
	o_params.verbose = p_verbose;

	//! Is first step?
	o_params.isFirstStep = (p_step == 1);

	//! Boost the paste trick
	o_params.doPasteBoost = true;
}

//	depend on sigma:
//		  sizePatch
//		  nSimilarPatches
//		  beta
//
// depend on sizePatch:
//		  nSimilarPatches (if channels == 3)
//		  sizeSearchWindow
//		  offSet
//		  tau
//
// depend on nSimilarPatches
//		  sizeSearchWindow
//
// depend on sizeSearchWindow
//		  boundary
//		  nSimilarPatches (cannot be more than total number of searchable patches)
//
// depend on sizeSearchTimeRangeFwd/Bwd
//		  nSimilarPatches

/**
 * @brief Sets size of spatial search window. It sets the border width accordingly,
 * and also ensures that the number of similar patches is not larger that the 
 * total number of available patches.
 *
 * @param prms             : nlbParams for first or second step of the algorithm;
 * @param sizeSearchWindow : size of search window;
 *
 * @return none.
 **/
void setSizeSearchWindow(nlbParams& prms, unsigned sizeSearchWindow)
{
	prms.sizeSearchWindow = sizeSearchWindow;
	prms.boundary = 2*(sizeSearchWindow/2) + (prms.sizePatch - 1);
	prms.nSimilarPatches = std::min(prms.nSimilarPatches, sizeSearchWindow *
	                                                      sizeSearchWindow *
	                                                     (prms.sizeSearchTimeRangeFwd +
	                                                      prms.sizeSearchTimeRangeBwd + 1));
}

/**
 * @brief Sets size of the patch. It sets the pixel offset as half the patch
 * size (this is BM3D speed-up).
 *
 * @param prms      : nlbParams for first or second step of the algorithm;
 * @param sizePatch : size of the patch;
 *
 * @return none.
 **/
void setSizePatch(nlbParams& prms, const VideoSize &size, unsigned sizePatch)
{
	prms.sizePatch = sizePatch;
	prms.boundary = 2*(prms.sizeSearchWindow/2) + (prms.sizePatch - 1);
	prms.offSet = sizePatch/2;

	//! Update number of similar patches, only if it is less than recommended value
	if (size.channels == 3)
//		prms.nSimilarPatches = std::max(prms.sizePatch * prms.sizePatch * 3, prms.nSimilarPatches);
		prms.nSimilarPatches = prms.sizePatch * prms.sizePatch * 3 * 
		                      (prms.sizeSearchTimeRangeFwd + 
		                       prms.sizeSearchTimeRangeBwd + 1);
	
}

/**
 * @brief Sets number of similar patches, ensuring that the number of similar
 * patches is not larger that the total number of available patches.
 *
 * @param prms            : nlbParams for first or second step of the algorithm;
 * @param nSimilarPatches : number of similar patches;
 *
 * @return none.
 **/
void setNSimilarPatches(nlbParams& prms, unsigned nSimilarPatches)
{
	prms.nSimilarPatches = nSimilarPatches;
	prms.nSimilarPatches = std::min(nSimilarPatches, prms.sizeSearchWindow *
	                                                 prms.sizeSearchWindow *
	                                                (prms.sizeSearchTimeRangeFwd +
	                                                 prms.sizeSearchTimeRangeBwd + 1));
}

/**
 * @brief Display parameters of the NL-Bayes algorithm.
 *
 * @param i_params : nlbParams for first or second step of the algorithm;
 *
 * @return none.
 **/
void printNlbParameters(
	const nlbParams &i_prms
){
//	printf("\nVideo NLBayes parameters\n");
//	printf("------------------------\n\n");
//	printf("Noise sigma = %g\n", i_prms1.sigma);

	printf("Parameters for step %d:\n", i_prms.isFirstStep ? 1 : 2);
	printf("\tPatch search:\n");
	printf("\t\tPatch size                  = %d\n"       , i_prms.sizePatch);
	printf("\t\tPatch size temporal         = %d\n"       , i_prms.sizePatchTime);
	printf("\t\tNumber of patches           = %d\n"       , i_prms.nSimilarPatches);
	if (!i_prms.isFirstStep) printf("\t\tDistance threshold (tau)    = %g\n"       , i_prms.tau);
	else                     printf("\t\tDistance threshold (tau)    = N/A\n"      );
	printf("\t\tSpatial search window       = %dx%d\n"    , i_prms.sizeSearchWindow, i_prms.sizeSearchWindow);
	printf("\t\tTemporal search range       = [-%d,%d]\n" , i_prms.sizeSearchTimeRangeBwd, i_prms.sizeSearchTimeRangeBwd);
	printf("\t\tSpatial border added        = %d\n"       , i_prms.boundary);
	printf("\tGroup filtering:\n");
	printf("\t\tBeta                        = %g\n"       , i_prms.beta);
	if (i_prms.useHomogeneousArea)
		printf("\t\tFlat area trick with gamma  = %g\n"       , i_prms.gamma);
	else
		printf("\t\tFlat area trick             = inactive\n");
	printf("\tSpeed-ups:\n");
	printf("\t\tOffset                      = %d\n"       , i_prms.offSet);
	printf("\t\tPasteBoost                  = %s\n\n"     , i_prms.doPasteBoost ? "active" : "inactive");
}

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * This function splits the image in several subimages when using
 * OpenMP. Each subimage is then assinged to a thread, which runs
 * processNlBayes on the subimage. The subimages are then assembled
 * into the final image.
 *
 * @param i_noisy: contains the noisy image;
 * @param o_basic: will contain the basic estimate image after the first step;
 * @param o_final: will contain the final denoised image after the second step;
 * @param p_useArea1 : if true, use the flat area trick for the first step;
 * @param p_useArea2 : if true, use the flat area trick for the second step;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if something wrong happens during the whole process.
 **/
int runNlBayes(
	Video<float> const& i_noisy
,	Video<float> &o_basic
,	Video<float> &o_final
,	const bool p_useArea1
,	const bool p_useArea2
,	const float p_sigma
,	const bool p_verbose
){
	//! Video size
	VideoSize size = i_noisy.sz;

	//! Parameters Initialization
	nlbParams p_prms1, p_prms2;
	initializeNlbParameters(p_prms1, 1, p_sigma, size, p_useArea1, p_verbose);
	initializeNlbParameters(p_prms2, 2, p_sigma, size, p_useArea2, p_verbose);

	//! NL-Bayes
	return runNlBayes(i_noisy, o_basic, o_final, p_prms1, p_prms2);
}

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * This function splits the image in several subimages when using
 * OpenMP. Each subimage is then assinged to a thread, which runs
 * processNlBayes on the subimage. The subimages are then assembled
 * into the final image.
 *
 * @param i_imNoisy: contains the noisy image;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_prms1  : parameters for first  step;
 * @param p_prms2  : parameters for second step;
 *
 * @return EXIT_FAILURE if something wrong happens during the whole process.
 **/
int runNlBayes(
	Video<float> const& i_imNoisy
,	Video<float> &o_imBasic
,	Video<float> &o_imFinal
,	const nlbParams p_prms1
,	const nlbParams p_prms2
){
	//! Only 1, 3 or 4-channels images can be processed.
	const unsigned chnls = i_imNoisy.sz.channels;
	if (! (chnls == 1 || chnls == 3 || chnls == 4))
		throw std::runtime_error("VideoNLB::runNlBayes: Wrong number of "
				"channels. Must be 1, 3 or 4.");

	//! Number of available cores
	unsigned nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
	if (p_prms1.verbose) printf("OpenMP is using %d threads\n", nThreads);
#endif
	const unsigned nParts = 2 * nThreads;

	//! Video size
	VideoSize imSize = i_imNoisy.sz;

	//! Initialization
	o_imBasic.resize(imSize);
	o_imFinal.resize(imSize);

	//! Print parameters
	if (p_prms1.verbose) printNlbParameters(p_prms1);
	if (p_prms2.verbose) printNlbParameters(p_prms2);

	//! Step 1
	{
		if (p_prms1.verbose) printf("1st Step\n");

		//! RGB to YUV
		Video<float> imNoisy = i_imNoisy;
		VideoUtils::transformColorSpace(imNoisy, true);

		//! Divide the noisy image into sub-images in order to easier parallelize the process
		std::vector<Video<float> > imNoisySub(nParts);
		std::vector<VideoUtils::CropPosition > imCrops(nParts);
		VideoUtils::subDivideTight(imNoisy, imNoisySub, imCrops, p_prms1.boundary, nParts);

		//! Process all sub-images
		std::vector<Video<float> > imBasicSub(nParts);
		std::vector<Video<float> > imFinalSub(nParts);
#ifdef _OPENMP
		// we make a copy of prms structure because, since it is constant,
		// it causes a compilation error with OpenMP (only on IPOL server)
		nlbParams prms1(p_prms1);
#pragma omp parallel for schedule(dynamic, nParts/nThreads) \
		shared(imNoisySub, imBasicSub, imFinalSub) \
		firstprivate (prms1)
#endif
		for (int n = 0; n < (int)nParts; n++)
			processNlBayes(imNoisySub[n], imBasicSub[n], imFinalSub[n], p_prms1, imCrops[n]);

		//! Get the basic estimate
		VideoUtils::subBuildTight(imBasicSub, o_imBasic, p_prms1.boundary);

		//! YUV to RGB
		VideoUtils::transformColorSpace(o_imBasic, false);

#ifdef DEBUG_SHOW_WEIGHT
		{
			std::vector<Video<float> > subWeights(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].origin_x;
				int part_y = imCrops[n].origin_y;
				int part_t = imCrops[n].origin_t;
				char name[1024];
				sprintf(name, "weight_step1_partx%dy%dt%d_%%03d.png", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subWeights[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> weight(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subWeights, weight, p_prms1.boundary);

			// Write to disk
			weight.saveVideo("weight_step1_%03d.png", 1);
		}
#endif

	}

	//! Step 2
	{
		if (p_prms2.verbose) printf("2nd Step\n");

		//! Divide the noisy and basic images into sub-images in order to easier parallelize the process
		std::vector<Video<float> > imNoisySub(nParts);
		std::vector<Video<float> > imBasicSub(nParts);
		std::vector<VideoUtils::CropPosition > imCrops(nParts);

		VideoUtils::subDivideTight(i_imNoisy, imNoisySub, imCrops, p_prms2.boundary, nParts);
		VideoUtils::subDivideTight(o_imBasic, imBasicSub, imCrops, p_prms2.boundary, nParts);

		//! Process all sub-images
		std::vector<Video<float> > imFinalSub(nParts);
#ifdef _OPENMP
		// we make a copy of prms structure because, since it is constant,
		// it causes a compilation error with OpenMP (only on IPOL server)
		nlbParams prms2(p_prms2);
#pragma omp parallel for schedule(dynamic, nParts/nThreads) \
		shared(imNoisySub, imBasicSub, imFinalSub) \
		firstprivate (prms2)
#endif
		for (int n = 0; n < (int) nParts; n++)
			processNlBayes(imNoisySub[n], imBasicSub[n], imFinalSub[n], p_prms2, imCrops[n]);

		//! Get the final result
		VideoUtils::subBuildTight(imFinalSub, o_imFinal, p_prms2.boundary);

#ifdef DEBUG_SHOW_WEIGHT
		{
			std::vector<Video<float> > subWeights(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].origin_x;
				int part_y = imCrops[n].origin_y;
				int part_t = imCrops[n].origin_t;
				char name[1024];
				sprintf(name, "weight_step2_partx%dy%dt%d_%%03d.png", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subWeights[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> weight(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subWeights, weight, p_prms1.boundary);

			// Write to disk
			weight.saveVideo("weight_step2_%03d.png", 1);
		}
#endif
	}

	return EXIT_SUCCESS;
}

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy image;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_params: see nlbParams.
 *
 * @return none.
 **/
void processNlBayes(
	Video<float> const& i_imNoisy
,	Video<float> &io_imBasic
,	Video<float> &o_imFinal
,	nlbParams const& p_params
,	VideoUtils::CropPosition p_crop
){
	using std::vector;

	//! Parameters initialization
	const bool step1 = p_params.isFirstStep;
	const unsigned sWx = p_params.sizeSearchWindow;
	const unsigned sWt = p_params.sizeSearchTimeRangeFwd +
	                     p_params.sizeSearchTimeRangeBwd + 1;// VIDEO
	const unsigned sPx = p_params.sizePatch;
	const unsigned sPt = p_params.sizePatchTime;

	unsigned nInverseFailed = 0;
	const float threshold = p_params.sigma * p_params.sigma * p_params.gamma *
	                       (p_params.isFirstStep ? i_imNoisy.sz.channels : 1.f);

	//! Weight sum per pixel // FIXME: only one channel
	Video<float> weight(i_imNoisy.sz.width, i_imNoisy.sz.height, i_imNoisy.sz.frames, 1, 0.f);

	//! Mask: true for pixels that still need to be processed
	Video<char> mask(i_imNoisy.sz.width, i_imNoisy.sz.height, i_imNoisy.sz.frames, 1, false);

	//! Only pixels of the center of the image must be processed (not the boundaries)
	int ori_x =                                  (p_crop.origin_x > 0                      ) ? sPx-1 + sWx/2 : 0 ; 
	int ori_y =                                  (p_crop.origin_y > 0                      ) ? sPx-1 + sWx/2 : 0 ; 
	int ori_f =                                  (p_crop.origin_t > 0                      ) ? sPt-1 + sWt/2 : 0 ; 
	int end_x = (int)i_imNoisy.sz.width  - (int)((p_crop.ending_x < p_crop.source_sz.width ) ? sPx-1 + sWx/2 : sPx-1); 
	int end_y = (int)i_imNoisy.sz.height - (int)((p_crop.ending_y < p_crop.source_sz.height) ? sPx-1 + sWx/2 : sPx-1); 
	int end_f = (int)i_imNoisy.sz.frames - (int)((p_crop.ending_t < p_crop.source_sz.frames) ? sPt-1 + sWt/2 : sPt-1); 
	for (int f = ori_f; f < end_f; f++)
	for (int y = ori_y; y < end_y; y++)
	for (int x = ori_x; x < end_x; x++)
		mask(x,y,f) = true;

	//! Used matrices during Bayes' estimate
	const unsigned patch_dim = step1 ? sPx * sPx * sPt : sPx * sPx * sPt * i_imNoisy.sz.channels;
	const unsigned patch_num = step1 ? p_params.nSimilarPatches : sWx * sWx * sWt;

	//! Matrices used for Bayes' estimate
	vector<unsigned> index(patch_num);
	matParams mat;
	mat.group3dTranspose.resize(patch_num * patch_dim);
	mat.tmpMat          .resize(patch_dim * patch_dim);
	mat.covMat          .resize(patch_dim * patch_dim);
	mat.covMatTmp       .resize(patch_dim * patch_dim);
	mat.baricenter      .resize(patch_dim);

	if (step1)
	{
		//! Allocate Sizes
		io_imBasic.resize(i_imNoisy.sz);

		//! Matrices used for Bayes' estimate
		vector<vector<float> > group3d(i_imNoisy.sz.channels,
				                         vector<float>(patch_num * patch_dim));

		for (unsigned ij = 0; ij < i_imNoisy.sz.whf; ij += p_params.offSet)
			if (mask(ij)) //< Only non-seen patches are processed
			{
				if (p_params.verbose && (ij % 10000 == 0))
				{
					printf("\rprocessing step1 %05.1f", (float)ij/(float)(i_imNoisy.sz.whf)*100.f);
					std::cout << std::flush;
				}

				const unsigned ij3 = (ij / i_imNoisy.sz.wh) * i_imNoisy.sz.whc + ij % i_imNoisy.sz.wh;

				//! Search for similar patches around the reference one
				unsigned nSimP = estimateSimilarPatchesStep1(i_imNoisy, group3d, index, ij3, p_params);

				//! If we use the homogeneous area trick
				bool doBayesEstimate = true;
				if (p_params.useHomogeneousArea)
					doBayesEstimate = !computeHomogeneousAreaStep1(group3d, sPx, patch_num,
							threshold, i_imNoisy.sz);

				//! Else, use Bayes' estimate
				if (doBayesEstimate)
					computeBayesEstimateStep1(group3d, mat, nInverseFailed, p_params, nSimP);

				//! Aggregation
				computeAggregationStep1(io_imBasic, weight, mask, group3d, index, p_params, nSimP);
//				//! TEST: 'temporal' aggregation
//				computeTemporalAggregationStep1(io_imBasic, weight, mask, group3d, index, p_params, nSimP);
			}

		if (p_params.verbose) printf("\n");

		//! Weighted aggregation
		computeWeightedAggregation(i_imNoisy, io_imBasic, weight);
	}
	else
	{
		//! Allocate Sizes
		o_imFinal.resize(i_imNoisy.sz);

		//! Matrices used for Bayes' estimate
		vector<float> group3dNoisy(patch_num * patch_dim);
		vector<float> group3dBasic(patch_num * patch_dim);

		for (unsigned ij = 0; ij < i_imNoisy.sz.whf; ij += p_params.offSet)
			if (mask(ij)) //< Only non-seen patches are processed
			{
				if (p_params.verbose && (ij % 10000 == 0))
				{
					printf("\rprocessing step2 %05.1f", (float)ij/(float)(i_imNoisy.sz.whf)*100.f);
					std::cout << std::flush;
				}

				const unsigned ij3 = (ij / i_imNoisy.sz.wh) * i_imNoisy.sz.whc + ij % i_imNoisy.sz.wh;

				//! Search for similar patches around the reference one
				unsigned nSimP = estimateSimilarPatchesStep2(i_imNoisy, io_imBasic,
						group3dNoisy, group3dBasic, index, ij3, p_params);

				//! If we use the homogeneous area trick
				bool doBayesEstimate = true;
				if (p_params.useHomogeneousArea)
					doBayesEstimate = !computeHomogeneousAreaStep2(group3dNoisy,
							group3dBasic, sPx, nSimP, threshold, i_imNoisy.sz);

				//! Else, use Bayes' estimate
				if (doBayesEstimate)
					computeBayesEstimateStep2(group3dNoisy, group3dBasic, mat,
							nInverseFailed, i_imNoisy.sz, p_params, nSimP);

				//! Aggregation
				computeAggregationStep2(o_imFinal, weight, mask, group3dNoisy,
						index, p_params, nSimP);
				//! TEST: 'temporal' aggregation
//				computeTemporalAggregationStep2(o_imFinal, weight, mask, group3dNoisy,
//						index, p_params, nSimP);
			}

		if (p_params.verbose) printf("\n");

		//! Weighted aggregation
		computeWeightedAggregation(i_imNoisy, o_imFinal, weight);
	}

	if (nInverseFailed > 0 && p_params.verbose)
		std::cout << "nInverseFailed = " << nInverseFailed << std::endl;

#ifdef DEBUG_SHOW_WEIGHT
	{
		int part_x = p_crop.origin_x;
		int part_y = p_crop.origin_y;
		int part_t = p_crop.origin_t;
		char name[1024];
		sprintf(name, "weight_step%d_partx%dy%dt%d_%%03d.png", step1 ? 1 : 2, part_x, part_y, part_t);
		weight.saveVideo(name, 1, 1);
	}
#endif
}

/**
 * @brief Estimate the best similar patches to a reference one.
 *
 * @param i_im: contains the noisy image on which distances are processed;
 * @param o_group3d: will contain values of similar patches;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
unsigned estimateSimilarPatchesStep1(
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

	rangex[0] = (int)px - (sWx-1)/2 - shift_x;
	rangey[0] = (int)py - (sWy-1)/2 - shift_y;
	ranget[0] = (int)pt -  sWt_b    - shift_t;

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

//	printf("distance.size() = %d", distance.size());
//	for (int i = 0; i < distance.size(); i++)
//	{
//		unsigned cx, cy, ct, cc;
//		i_im.sz.coords(distance[i].second, cx, cy, ct, cc);
//		printf("d[%03d] = %g - p = [%02d,%02d,%02d,%2d]\n", i, distance[i].first,
//				cx, cy, ct, cc);
//	}

	//! Keep only the N2 best similar patches
	const unsigned nSimP = std::min(p_params.nSimilarPatches, (unsigned)distance.size());
	std::partial_sort(distance.begin(), distance.begin() + nSimP,
	                  distance.end(), comparaisonFirst);

//	for (int i = 0; i < nSimP; i++)
//	{
//		unsigned cx, cy, ct, cc;
//		i_im.sz.coords(distance[i].second, cx, cy, ct, cc);
//		printf("d[%03d] = %g - p = [%02d,%02d,%02d,%2d]\n", i, distance[i].first,
//				cx, cy, ct, cc);
//	}


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

	return nSimP;
}

/**
 * @brief Keep from all near patches the similar ones to the reference patch for the second step.
 *
 * @param i_imNoisy: contains the original noisy image;
 * @param i_imBasic: contains the basic estimation;
 * @param o_group3dNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_group3dBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_imSize: size of images;
 * @param p_params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep2(
	Video<float> const& i_imNoisy
,	Video<float> const& i_imBasic
,	std::vector<float> &o_group3dNoisy
,	std::vector<float> &o_group3dBasic
,	std::vector<unsigned> &o_index
,	const unsigned pidx
,	const nlbParams &p_params
){
	//! Initialization
	const int chnls = i_imNoisy.sz.channels;
	const int sPx   = p_params.sizePatch;
	const int sPt   = p_params.sizePatchTime;
	int sWx   = p_params.sizeSearchWindow;
	int sWy   = p_params.sizeSearchWindow;
	const int sWt_f = p_params.sizeSearchTimeRangeFwd;
	const int sWt_b = p_params.sizeSearchTimeRangeBwd;

	//! Coordinates of center of search box
	unsigned px, py, pt, pc;
	i_imBasic.sz.coords(pidx, px, py, pt, pc);

	unsigned rangex[2];
	unsigned rangey[2];
	unsigned ranget[2];

#ifdef CENTRED_SEARCH
	rangex[0] = std::max(0, (int)px - (sWx-1)/2);
	rangey[0] = std::max(0, (int)py - (sWx-1)/2);
	ranget[0] = std::max(0, (int)pt -  sWt_b   );

	rangex[1] = std::min((int)i_imNoisy.sz.width  - sPx, (int)px + (sWx-1)/2);
	rangey[1] = std::min((int)i_imNoisy.sz.height - sPx, (int)py + (sWx-1)/2);
	ranget[1] = std::min((int)i_imNoisy.sz.frames - sPt, (int)pt +  sWt_f   );
#else
	int shift_x = std::min(0, (int)px - (sWx-1)/2); 
	int shift_y = std::min(0, (int)py - (sWx-1)/2); 
	int shift_t = std::min(0, (int)pt -  sWt_b   ); 

	rangex[0] = (int)px - (sWx-1)/2 - shift_x;
	rangey[0] = (int)py - (sWx-1)/2 - shift_y;
	ranget[0] = (int)pt -  sWt_b    - shift_t;

	rangex[1] = std::min((int)i_imNoisy.sz.width  - sPx, (int)px + (sWx-1)/2 - shift_x);
	rangey[1] = std::min((int)i_imNoisy.sz.height - sPx, (int)py + (sWx-1)/2 - shift_y);
	ranget[1] = std::min((int)i_imNoisy.sz.frames - sPt, (int)pt +  sWt_f    - shift_t);
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

	//! Save index of similar patches
	const float threshold = (p_params.tau > distance[p_params.nSimilarPatches - 1].first ?
	                         p_params.tau : distance[p_params.nSimilarPatches - 1].first);
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

	return nSimP;
}

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group3d: contains for each channels values of similar patches. If an homogeneous area
 *			is detected, will contain the average of all pixels in similar patches;
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_doLinearRegression: if true, apply a linear regression to average value of pixels;
 * @param p_imSize: size of the image.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep1(
	std::vector<std::vector<float> > &io_group3d
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
){
	//! Initialization
	const unsigned N = p_sP * p_sP * p_nSimP;

	//! Compute the standard deviation of the set of patches
	float stdDev = 0.f;
	for (unsigned c = 0; c < p_imSize.channels; c++)
		stdDev += computeStdDeviation(io_group3d[c], p_sP * p_sP, p_nSimP, 1);

	//! If we are in an homogeneous area
	if (stdDev < p_threshold)
	{
		for (unsigned c = 0; c < p_imSize.channels; c++)
		{
			float mean = 0.f;

			for (unsigned k = 0; k < N; k++) mean += io_group3d[c][k];
			mean /= (float) N;

			for (unsigned k = 0; k < N; k++) io_group3d[c][k] = mean;
		}
		return 1;
	}
	else return 0;
}

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group3dNoisy: inputs values of similar patches for the noisy video;
 *                         if the area is classified as homogeneous, outputs the
 *                         average of all pixels in all patches.
 * @param i_group3dBasic: contains values of similar patches for the basic video.
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_imSize: size of the video.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep2(
	std::vector<float> &io_group3dNoisy
,	std::vector<float> const &i_group3dBasic
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
){
	//! Parameters
	const unsigned sP2 = p_sP * p_sP;
	const unsigned sPC = sP2 * p_imSize.channels;

	//! Compute the standard deviation of the set of patches
	const float stdDev = computeStdDeviation(io_group3dNoisy, sP2, p_nSimP, p_imSize.channels);

	//! If we are in an homogeneous area
	if (stdDev < p_threshold)
	{
		for (unsigned c = 0; c < p_imSize.channels; c++)
		{
				float mean = 0.f;
				for (unsigned n = 0; n < p_nSimP; n++)
				for (unsigned k = 0; k < sP2; k++)
					mean += i_group3dBasic[n * sPC + c * sP2 + k];

				mean /= float(sP2 * p_nSimP);

				for (unsigned n = 0; n < p_nSimP; n++)
				for (unsigned k = 0; k < sP2; k++)
					io_group3dNoisy[n * sPC + c * sP2 + k] = mean;
		}
		return 1;
	}
	else return 0;
}

/**
 * @brief Compute the Bayes estimation.
 *
 * @param io_group3d: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *    - group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *    - baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 *    - covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *    - covMatTmp: allocated memory. Used to process the Bayes estimate;
 *    - tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeBayesEstimateStep1(
	std::vector<std::vector<float> > &io_group3d
,	matParams &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters
	const unsigned chnls = io_group3d.size();
	const unsigned sP2   = p_params.sizePatch * p_params.sizePatch * p_params.sizePatchTime;
	const float valDiag  = p_params.beta * p_params.sigma * p_params.sigma;

#ifndef DEBUG_SHOW_PATCH_GROUPS
	//! Bayes estimate
	if (p_nSimP >= sP2)
	{
		for (unsigned c = 0; c < chnls; c++)
		{
			//! Center data around the baricenter
			centerData(io_group3d[c], i_mat.baricenter, p_nSimP, sP2);

			//! Compute the covariance matrix of the set of similar patches
			covarianceMatrix(io_group3d[c], i_mat.covMat, p_nSimP, sP2);

			//! Bayes' Filtering
			if (inverseMatrix(i_mat.covMat, sP2) == EXIT_SUCCESS)
			{
				productMatrix(i_mat.group3dTranspose, i_mat.covMat, io_group3d[c],
								  sP2, sP2, p_nSimP);
				for (unsigned k = 0; k < sP2 * p_nSimP; k++)
					io_group3d[c][k] -= valDiag * i_mat.group3dTranspose[k];
			}
			else io_nInverseFailed++;

			//! Add baricenter
			for (unsigned j = 0, k = 0; j < sP2; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					io_group3d[c][k] += i_mat.baricenter[j];
		}
	}
	else io_nInverseFailed++;
#else
	//! Show color-code representing patch status
	if (p_nSimP >= sP2)
	{
		float color[3];
		// green
		color[0] =  147.2243f;
		color[1] =  0;
		color[2] = -208.2066f;

		for (unsigned c = 0; c < chnls; c++)
		{
			//! Center data around the baricenter
			centerData(io_group3d[c], i_mat.baricenter, p_nSimP, sP2);

			//! Compute the covariance matrix of the set of similar patches
			covarianceMatrix(io_group3d[c], i_mat.covMat, p_nSimP, sP2);

			//! Bayes' Filtering
			if (inverseMatrix(i_mat.covMat, sP2) == EXIT_FAILURE)
			{
				io_nInverseFailed++;
				// blue
				color[0] =  147.2243f;
				color[1] = -180.3122f;
				color[2] =  104.1033f;
			}
		}

		for (unsigned k = 0; k < sP2 * p_nSimP; k++)
		{
			io_group3d[0][k] =  color[0];
			io_group3d[1][k] =  color[1];
			io_group3d[2][k] =  color[2];
		}
	}
	else
	{
		io_nInverseFailed++;
		for (unsigned k = 0; k < sP2 * p_nSimP; k++)
		{
			// red
			io_group3d[0][k] = 147.2243f;
			io_group3d[1][k] = 180.3122f;
			io_group3d[2][k] = 104.1033f;
		}
	}
#endif
}

/**
 * @brief Compute the Bayes estimation.
 *
 * @param io_group3dNoisy: inputs all similar patches in the noisy image,
 *                         outputs their denoised estimates.
 * @param i_group3dBasic: contains all similar patches in the basic image.
 * @param i_mat: contains :
 *    - group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *    - baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
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
void computeBayesEstimateStep2(
	std::vector<float> &io_group3dNoisy
,	std::vector<float>  &i_group3dBasic
,	matParams &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float diagVal = p_params.beta * p_params.sigma * p_params.sigma;
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime * p_imSize.channels;

	//! Center 3D groups around their baricenter
	centerData( i_group3dBasic, i_mat.baricenter, p_nSimP, sPC);
	centerData(io_group3dNoisy, i_mat.baricenter, p_nSimP, sPC);

	//! Compute the covariance matrix of the set of similar patches
	covarianceMatrix(i_group3dBasic, i_mat.covMat, p_nSimP, sPC);

	//! Bayes' Filtering
	for (unsigned k = 0; k < sPC; k++)
		i_mat.covMat[k * sPC + k] += diagVal;

	//! Compute the estimate
	if (inverseMatrix(i_mat.covMat, sPC) == EXIT_SUCCESS)
	{
		productMatrix(i_group3dBasic, i_mat.covMat, io_group3dNoisy, sPC, sPC, p_nSimP);
		for (unsigned k = 0; k < sPC * p_nSimP; k++)
			io_group3dNoisy[k] -= diagVal * i_group3dBasic[k];
	}
	else 
		io_nInverseFailed++;

	//! Add baricenter
	for (unsigned j = 0, k = 0; j < sPC; j++)
		for (unsigned i = 0; i < p_nSimP; i++, k++)
			io_group3dNoisy[k] += i_mat.baricenter[j];
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeAggregationStep1(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<std::vector<float> > const& i_group3d
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls = io_im.sz.channels;
	const unsigned sPx   = p_params.sizePatch;
	const unsigned sPt   = p_params.sizePatchTime;

	const unsigned w   = io_im.sz.width;
	const unsigned h   = io_im.sz.height;
	const unsigned wh  = io_im.sz.wh;
	const unsigned whc = io_im.sz.whc;

#ifndef DEBUG_SHOW_PATCH_GROUPS
	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / whc) * wh + ind % wh;
		for (unsigned pt = 0; pt < sPt; pt++)
		for (unsigned py = 0; py < sPx; py++)
		for (unsigned px = 0; px < sPx; px++)
		{
			for (unsigned c = 0; c < chnls; c++)
			{
				const unsigned ij = ind  + c * wh;
				io_im(ij + pt * whc + py * w + px) += 
					i_group3d[c][(pt * sPx*sPx + py * sPx + px) * p_nSimP + n];
			}
			io_weight(ind1 + pt * wh + py * w + px)++;
		}

		//! Use Paste Trick
		io_mask(ind1) = false;

		unsigned px, py, pt;
		io_mask.sz.coords(ind1, px, py, pt);

		if (p_params.doPasteBoost)
		{
			if (py >     2*sPx) io_mask(ind1 - w) = false;
			if (py < w - 2*sPx) io_mask(ind1 + w) = false;
			if (px >     2*sPx) io_mask(ind1 - 1) = false;
			if (px > h - 2*sPx) io_mask(ind1 + 1) = false;
		}
	}
#else
	for (unsigned n = 0; n < 1; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / whc) * wh + ind % wh;
		for (unsigned pt = 0; pt < sPt; pt++)
		for (unsigned py = 0; py < sPx; py++)
		for (unsigned px = 0; px < sPx; px++)
		{
			for (unsigned c = 0; c < chnls; c++)
			{
				const unsigned ij = ind  + c * wh;
				io_im(ij + pt * whc + py * w + px) += 
					i_group3d[c][(pt * sPx*sPx + py * sPx + px) * p_nSimP + n];
			}
			io_weight(ind1 + pt * wh + py * w + px)++;
		}
	}

	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / whc) * wh + ind % wh;

		//! Use Paste Trick
		io_mask(ind1) = false;

		unsigned px, py, pt;
		io_mask.sz.coords(ind1, px, py, pt);

		if (p_params.doPasteBoost)
		{
			if (py >     2*sPx) io_mask(ind1 - w) = false;
			if (py < w - 2*sPx) io_mask(ind1 + w) = false;
			if (px >     2*sPx) io_mask(ind1 - 1) = false;
			if (px > h - 2*sPx) io_mask(ind1 + 1) = false;
		}
	}
#endif
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeAggregationStep2(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<float> const& i_group3d
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls = io_im.sz.channels;
	const unsigned sPx   = p_params.sizePatch;
	const unsigned sPt   = p_params.sizePatchTime;

	const unsigned w   = io_im.sz.width;
	const unsigned h   = io_im.sz.height;
	const unsigned wh  = io_im.sz.wh;
	const unsigned whc = io_im.sz.whc;

	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			const unsigned ij = ind + c * wh;
			for (unsigned pt = 0; pt < sPt; pt++)
			for (unsigned py = 0; py < sPx; py++)
			for (unsigned px = 0; px < sPx; px++, k++)
				io_im(ij + pt * whc + py * w + px) +=
					i_group3d[k * p_nSimP + n];
		}

		const unsigned ind1 = (ind / whc) * wh + ind % wh;
		for (unsigned pt = 0; pt < sPt; pt++)
		for (unsigned py = 0; py < sPx; py++)
		for (unsigned px = 0; px < sPx; px++)
			io_weight(ind1 + pt * wh + py * w + px)++;

		//! Apply Paste Trick
		io_mask(ind1) = false;

		unsigned px, py, pt;
		io_mask.sz.coords(ind1, px, py, pt);

		if (p_params.doPasteBoost)
		{
			if (py >   0) io_mask(ind1 - w) = false;
			if (py < w-1) io_mask(ind1 + w) = false;
			if (px >   0) io_mask(ind1 - 1) = false;
			if (px > h-1) io_mask(ind1 + 1) = false;
		}
	}
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D
 * group. This version is for a test: in the original version, all patches
 * in the group are marked as processed, and cannot be origins of a patch
 * group. In this version we only mark as processed the patches of the 
 * group which are nearby frames to the group origin.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeTemporalAggregationStep1(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<std::vector<float> > const& i_group3d
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls  = io_im.sz.channels;
	const unsigned width  = io_im.sz.width;
	const unsigned height = io_im.sz.height;
	const unsigned sP     = p_params.sizePatch;

	//! Compute coordinates of group origin
	int FRAME_STEP = 2;
	unsigned ox, oy, ot, oc;
	io_im.sz.coords(i_index[0], ox, oy, ot, oc);
	
#ifndef DEBUG_SHOW_PATCH_GROUPS
	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / io_im.sz.whc) * io_im.sz.wh + ind % io_im.sz.wh;
		for (unsigned p = 0; p < sP; p++)
		for (unsigned q = 0; q < sP; q++)
		{
			for (unsigned c = 0; c < chnls; c++)
			{
				const unsigned ij = ind  + c * width * height;
				io_im(ij + p * width + q) += i_group3d[c][(p * sP + q) * p_nSimP + n];
			}
			io_weight(ind1 + p * width + q)++;
		}

		// compute coordinates of current pixel
		unsigned px, py, pt, pc;
		io_im.sz.coords(ind, px, py, pt, pc);

		if (abs((int)pt - (int)ot) < FRAME_STEP)
		{
			//! Use Paste Trick
			io_mask(ind1) = false;

			if (p_params.doPasteBoost)
			{
				io_mask(ind1 - width) = false;
				io_mask(ind1 + width) = false;
				io_mask(ind1 - 1    ) = false;
				io_mask(ind1 + 1    ) = false;
			}
		}
	}
#else
	for (unsigned n = 0; n < 1; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / io_im.sz.whc) * io_im.sz.wh + ind % io_im.sz.wh;
		for (unsigned p = 0; p < 1; p++)
		for (unsigned q = 0; q < 1; q++)
		{
			for (unsigned c = 0; c < chnls; c++)
			{
				const unsigned ij = ind  + c * width * height;
				io_im(ij + p * width + q) += i_group3d[c][(p * sP + q) * p_nSimP + n];
			}
			io_weight(ind1 + p * width + q)++;
		}
	}

	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / io_im.sz.whc) * io_im.sz.wh + ind % io_im.sz.wh;

		// compute coordinates of current pixel
		unsigned px, py, pt, pc;
		io_im.sz.coords(ind, px, py, pt, pc);

		if (abs((int)pt - (int)ot) < FRAME_STEP)
		{
			//! Use Paste Trick
			io_mask(ind1) = false;

			if (p_params.doPasteBoost)
			{
				io_mask(ind1 - width) = false;
				io_mask(ind1 + width) = false;
				io_mask(ind1 - 1    ) = false;
				io_mask(ind1 + 1    ) = false;
			}
		}
	}
#endif
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D
 * group. This version is for a test: in the original version, all patches
 * in the group are marked as processed, and cannot be origins of a patch
 * group. In this version we only mark as processed the patches of the 
 * group which are nearby frames to the group origin.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeTemporalAggregationStep2(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<float> const& i_group3d
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
){
	//! Parameters initializations
	const unsigned chnls = io_im.sz.channels;
	const unsigned width = io_im.sz.width;
	const unsigned wh    = width * io_im.sz.height;
	const unsigned sP    = p_params.sizePatch;

	//! Compute coordinates of group origin
	int FRAME_STEP = 2;
	unsigned ox, oy, ot, oc;
	io_im.sz.coords(i_index[0], ox, oy, ot, oc);
	
	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind  = i_index[n];
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			const unsigned ij = ind + c * wh;
			for (unsigned p = 0; p < sP; p++)
			for (unsigned q = 0; q < sP; q++, k++)
				io_im(ij + p * width + q) += i_group3d[k * p_nSimP + n];
		}

		const unsigned ind1 = (ind / io_im.sz.whc) * io_im.sz.wh + ind % io_im.sz.wh;
		for (unsigned p = 0; p < sP; p++)
		for (unsigned q = 0; q < sP; q++)
			io_weight(ind1 + p * width + q)++;

		// compute coordinates of current pixel
		unsigned px, py, pt, pc;
		io_im.sz.coords(ind, px, py, pt, pc);

		if (abs((int)pt - (int)ot) < FRAME_STEP)
		{
			//! Apply Paste Trick
			io_mask(ind1) = false;

			if (p_params.doPasteBoost)
			{
				io_mask(ind1 - width) = false;
				io_mask(ind1 + width) = false;
				io_mask(ind1 - 1    ) = false;
				io_mask(ind1 + 1    ) = false;
			}
		}
	}
}

/**
 * @brief Compute the final weighted aggregation.
 *
 * i_im: image of reference, when the weight if null;
 * io_im: will contain the final image;
 * i_weight: associated weight for each estimate of pixels.
 *
 * @return : none.
 **/
void computeWeightedAggregation(
	Video<float> const& i_im
,	Video<float> &io_im
,	Video<float> const& i_weight
){
	for (unsigned f = 0; f < i_im.sz.frames  ; f++)
	for (unsigned c = 0; c < i_im.sz.channels; c++)
	for (unsigned y = 0; y < i_im.sz.height  ; y++)
	for (unsigned x = 0; x < i_im.sz.width   ; x++)
		if (i_weight(x,y,f) > 0.f) io_im(x,y,f,c) /= i_weight(x,y,f);
		else                       io_im(x,y,f,c) /= i_im(x,y,f,c);

}

}
