/*
 * Modified work Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
 * Original work Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
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
 * @file nldct.cpp
 * @brief NL-DCT image & video denoising functions
 *
 * @author Pablo Arias <pariasm@gmail.com>
 *
 * Based on NL-Bayes code from
 * Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include <stdexcept>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <math.h>
#include <float.h>

#include <unordered_set>

#include "nldct.h"
#include "matrix_funs.h"
#include "utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define DCT_BASIS
//#define DCT_DONT_CENTER1
//#define DCT_DONT_CENTER2

//#define USE_FFTW

/* Shrinks the mean to 0 with an empirical hyper-Bayesian approach. */
//#define MEAN_HYPERPRIOR1
//#define MEAN_HYPERPRIOR2
#define MEAN_HYPERPRIOR_BM3D1
#define MEAN_HYPERPRIOR_BM3D2

/* Avoid negative weights in the empirical Wiener filter. When the estimated
 * variance of a certain component is lower than the noise variance, the filter
 * coefficient is set to zero. This applies whenever the Gaussian model is
 * estimated from the noisy patches: in the first step, or in the second step
 * if NOISY_COVARIANCE2 option is defined. */
#define THRESHOLD_WEIGHTS1
//#define THRESHOLD_WEIGHTS2

/* Uses an adaptation of Li,Zhand,Dai fixed point iteration to estimate the
 * signal power in the empirical Wiener filter. It applies whenever the Gaussian
 * model is learnt from the noisy patches, but currently it is implemented only 
 * in the second step (it should always be used together with the NOISY_COVARIANCE2
 * option defined). */
//#define LI_ZHANG_DAI1
//#define LI_ZHANG_DAI2

/* Use nonlinear coefficient thresholding instead of empirical Wiener filter in
 * Step 1. Parameter beta is used to control the threshold beta*sigma². */
#define THRESHOLDING1
//#define SOFT_THRESHOLD1
//#define SOFT_THRESHOLD1_BAYES 

/* Use nonlinear coefficient thresholding instead of empirical Wiener filter in
 * Step 2. Parameter beta is used to control the threshold beta*sigma². */
//#define THRESHOLDING2
//#define SOFT_THRESHOLD2
//#define SOFT_THRESHOLD2_BAYES 

/* Use a linear coefficient thresholding instead of empirical Wiener filter in
 * Step 1. The thresholding is applied to a certain coefficient based on its 
 * sample variance over the group of similar patches. The threshold is linear
 * since it is actually a 0-1 filter. Each patch is processed with the same
 * threshold. Parameter beta is used to control the threshold beta*sigma².
 * This option is disabled if THRESHOLDING1 is defined. */
//#define LINEAR_THRESHOLDING1
//#define LINEAR_HARD_THRESHOLDING2
//#define LINEAR_SOFT_THRESHOLDING2

/* Corrects the 'centering bug' discovered by Nicola. In the second step, basic
 * and noisy patches are centered using the basic baricenter. If left undefined,
 * each set of patches (noisy and basic) are centered using their own baricenter. */
//#define BARICENTER_BASIC

/* The parameter beta is used as a noise correction factor. It provides a way 
 * to control the thresholding/filtering strength of the algorithm, by modifying
 * sigma as beta * sigma. By default, this modifications should only be applied 
 * to the sigma in the threshold/filter operator. If the following flag is defined,
 * the modified sigma is also used to estimate the a priori variances in the first
 * step. This option causes a 0.05dB increase in the results for a particular value
 * of beta. However, the resultng method seems much more sensitive to beta. 
 */
//#define USE_BETA_FOR_VARIANCE
//TODO The Bayes estimation with non-fixed DCT basis uses by default the modified
//TODO beta in the variance. We should test if applying it only to the filter is
//TODO better.

/* Compute the 2nd step covariance matrix from the noisy patches. In this way,
 * the basic estimate is used only in the computation of the patch distances.*/
//#define NOISY_COVARIANCE2

/* Use Gaussian (of the group variance) decay for the group aggregation weights.
 * The default is a logistic threshold. */
//#define GAUSSIAN_GROUP_AGGREGATION_WEIGHTS


//#define CENTRED_SEARCH

// colors
#define ANSI_BLK  "\x1b[30m"
#define ANSI_RED  "\x1b[31m"
#define ANSI_GRN  "\x1b[32m"
#define ANSI_YLW  "\x1b[33m"
#define ANSI_BLU  "\x1b[34m"
#define ANSI_MAG  "\x1b[35m"
#define ANSI_CYN  "\x1b[36m"
#define ANSI_WHT  "\x1b[37m"
#define ANSI_BBLK "\x1b[30;01m"
#define ANSI_BRED "\x1b[31;01m"
#define ANSI_BGRN "\x1b[32;01m"
#define ANSI_BYLW "\x1b[33;01m"
#define ANSI_BBLU "\x1b[34;01m"
#define ANSI_BMAG "\x1b[35;01m"
#define ANSI_BCYN "\x1b[36;01m"
#define ANSI_BWHT "\x1b[37;01m"
#define ANSI_RST "\x1b[0m"

namespace VideoNLB
{

#ifdef USE_FFTW
// define global dct thread handler
DCTThreadsHandler globalDCT;
#endif

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
,	const unsigned rank
){
	const bool s1 = (p_step == 1);

	o_params.transform = (p_step == 1) ? bior1_5 : dct;
//	o_params.transform = dct;

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
	o_params.offSet = (p_step == 2) ? std::max((unsigned)1, 3 * o_params.sizePatch / 4)
	                                : std::max((unsigned)1, 2 * o_params.sizePatch / 4);

//	o_params.offSetTime = std::max((unsigned)1, o_params.sizePatchTime / 2);
	o_params.offSetTime = 1;

	//! Use the homogeneous area detection trick
	o_params.useHomogeneousArea = p_flatArea;

	//! Size of the search window around the reference patch (must be odd)
	o_params.sizeSearchWindow = o_params.nSimilarPatches / 2; // ASK MARC and 7*sizePatch1 in IPOL manuscript
	if (o_params.sizeSearchWindow % 2 == 0) o_params.sizeSearchWindow++;

	//! Search window, temporal search radii
	o_params.sizeSearchTimeRangeFwd = timeSearchRangeFwd;
	o_params.sizeSearchTimeRangeBwd = timeSearchRangeBwd;
	o_params.nSimilarPatches *= timeSearchRangeFwd + timeSearchRangeBwd + 1;

	//! Size of boundaries used during the sub division
	o_params.boundary = 2*(o_params.sizeSearchWindow/2) + (o_params.sizePatch - 1);

	//! Parameter used to determine if an area is homogeneous
	o_params.gamma = 1.05f;

	//! Parameter used to estimate the covariance matrix
	if (p_size.channels == 1)
	{
		o_params.beta = 1.f;
//		if(s1) o_params.beta = (p_sigma < 15.f ? 1.1f : (p_sigma < 70.f ? 1.f : 0.9f));
//		else   o_params.beta = (p_sigma < 15.f ? 1.1f : (p_sigma < 35.f ? 1.f : 0.9f));
	}
	else
	{
		o_params.beta = 1.f;
//		if(s1) o_params.beta = 1.f;
//		else   o_params.beta = (p_sigma < 50.f ? 1.2f : 1.f);
	}

	//! Noise correction factor in case of filtering the group baricenter
	o_params.betaMean = 1.f;

	//! Maximum rank of covariance matrix
	o_params.rank = rank;

	//! Decay of per-patch aggregation weights
	o_params.aggreGammaPatch = 0.f;

	//! Decay of per-group aggregation weights
	o_params.aggreGammaGroup = 0.f;

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

	//! Color space to be used
	o_params.colorSpace = YUV;
}

/*	depend on sigma:
 *		  sizePatch
 *		  nSimilarPatches
 *		  beta
 *
 * depend on sizePatch:
 *		  nSimilarPatches (if channels == 3)
 *		  sizeSearchWindow
 *		  offSet
 *		  tau
 *
 * depend on nSimilarPatches
 *		  sizeSearchWindow
 *
 * depend on sizeSearchWindow
 *		  boundary
 *		  nSimilarPatches (cannot be more than total number of searchable patches)
 *
 * depend on sizeSearchTimeRangeFwd/Bwd
 *		  nSimilarPatches
 */

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
	int sizePatch_old = prms.sizePatch;
	prms.sizePatch = sizePatch;
	prms.boundary = 2*(prms.sizeSearchWindow/2) + (prms.sizePatch - 1);
	prms.offSet = sizePatch/2;
	prms.offSet = (prms.isFirstStep) ? std::max((unsigned)1, 3 * sizePatch / 4)
	                                 : std::max((unsigned)1, 2 * sizePatch / 4);

	//! Update number of similar patches, only if it is less than recommended value
	if (size.channels == 3)
//		prms.nSimilarPatches = std::max(prms.sizePatch * prms.sizePatch * 3, prms.nSimilarPatches);
		prms.nSimilarPatches = prms.sizePatch * prms.sizePatch * 3 * 
		                      (prms.sizeSearchTimeRangeFwd + 
		                       prms.sizeSearchTimeRangeBwd + 1);

	if (!prms.isFirstStep)
		prms.tau *= (float)(prms.sizePatch * prms.sizePatch) 
		          / (float)( sizePatch_old *  sizePatch_old);
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
 * @brief Sets the distance threshold relative to the patch size.
 *
 * @param prms : nlbParams for first or second step of the algorithm;
 * @param tau  : distance threshold;
 *
 * @return none.
 **/
void setTau(nlbParams& prms, const VideoSize &size, float tau)
{
	prms.tau = tau * tau * prms.sizePatch * prms.sizePatch * 
	           prms.sizePatchTime * (prms.isFirstStep ? 1 : size.channels);
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

	printf("\x1b[37;01m" "Parameters for step %d:" ANSI_RST "\n" , i_prms.isFirstStep ? 1 : 2);
	printf("\tPatch search:\n");
	printf("\t\tPatch size                  = %d\n"       , i_prms.sizePatch);
	printf("\t\tPatch size temporal         = %d\n"       , i_prms.sizePatchTime);
	printf("\t\tNumber of patches           = %d\n"       , i_prms.nSimilarPatches);
	printf("\t\tDistance threshold (tau)    = %g\n"       , i_prms.tau);
	printf("\t\tSpatial search window       = %dx%d\n"    , i_prms.sizeSearchWindow, i_prms.sizeSearchWindow);
	printf("\t\tTemporal search range       = [-%d,%d]\n" , i_prms.sizeSearchTimeRangeBwd, i_prms.sizeSearchTimeRangeBwd);
	printf("\t\tSpatial border added        = %d\n"       , i_prms.boundary);
	printf("\tGroup filtering:\n");
	printf("\t\tBeta                        = %g\n"       , i_prms.beta);
	printf("\t\tRank                        = %d\n"       , i_prms.rank);
	printf("\t\tPatch aggre. weights decay  = %g\n"       , i_prms.aggreGammaPatch);
	printf("\t\tGroup aggre. weights decay  = %g\n"       , i_prms.aggreGammaGroup);
#if defined(MEAN_HYPERPRIOR1) || defined(MEAN_HYPERPRIOR2)
	printf("\t\tBeta (mean)                 = %g\n"       , i_prms.betaMean);
#endif
	if (i_prms.useHomogeneousArea)
		printf("\t\tFlat area trick with gamma  = %g\n"       , i_prms.gamma);
	else
		printf("\t\tFlat area trick             = inactive\n");
	printf("\tSpeed-ups:\n");
	printf("\t\tOffset                      = %d\n"       , i_prms.offSet);
	printf("\t\tOffsetTime                  = %d\n"       , i_prms.offSetTime);
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
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> const& i_noisy
,	Video<float> &o_basic
,	Video<float> &o_final
,	const bool p_useArea1
,	const bool p_useArea2
,	const float p_sigma
,	const bool p_verbose
#ifdef DEBUG_COMPUTE_GROUP_ERROR
,	Video<float> & i_imClean
#endif
){
	//! Video size
	VideoSize size = i_noisy.sz;

	//! Parameters Initialization
	nlbParams p_prms1, p_prms2;
	initializeNlbParameters(p_prms1, 1, p_sigma, size, p_useArea1, p_verbose);
	initializeNlbParameters(p_prms2, 2, p_sigma, size, p_useArea2, p_verbose);

	//! NL-Bayes
#ifndef DEBUG_COMPUTE_GROUP_ERROR
	return runNlBayes(i_noisy, o_basic, o_final, p_prms1, p_prms2);
#else
	return runNlBayes(i_noisy, o_basic, o_final, p_prms1, p_prms2, i_imClean);
#endif
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
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> const &i_imNoisy
,	Video<float> &o_imBasic
,	Video<float> &o_imFinal
,	const nlbParams p_prms1
,	const nlbParams p_prms2
#ifdef DEBUG_COMPUTE_GROUP_ERROR
,  Video<float> &i_imClean
#endif
){
#ifndef DEBUG_COMPUTE_GROUP_ERROR
	return runNlBayes(i_imNoisy, Video<float>(), Video<float>(), 
			o_imBasic, o_imFinal, p_prms1, p_prms2);
#else
	return runNlBayes(i_imNoisy, Video<float>(), Video<float>(), 
			o_imBasic, o_imFinal, p_prms1, p_prms2, i_imClean);
#endif
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
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> const& i_imNoisy
,  Video<float> const& i_fflow
,  Video<float> const& i_bflow
,	Video<float> &o_imBasic
,	Video<float> &o_imFinal
,	const nlbParams p_prms1
,	const nlbParams p_prms2
#ifdef DEBUG_COMPUTE_GROUP_ERROR
,	Video<float> & i_imClean
#endif
){
//	printf("Noisy video: %dx%dx%d channels %d\n", i_imNoisy.sz.width,
//	                                              i_imNoisy.sz.height,
//	                                              i_imNoisy.sz.frames,
//	                                              i_imNoisy.sz.channels);
	//! Only 1, 3 or 4-channels images can be processed.
	const unsigned chnls = i_imNoisy.sz.channels;
	if (! (chnls == 1 || chnls == 3 || chnls == 4))
		throw std::runtime_error("VideoNLB::runNlBayes: Wrong number of "
				"channels. Must be 1, 3 or 4.");

	//! Check if optical flow is valid
	const bool use_oflow = (i_fflow.sz.width > 0);
	if (use_oflow && 
		(i_fflow.sz.channels != 2 || 
		 i_fflow.sz.width  != i_imNoisy.sz.width  ||
		 i_fflow.sz.height != i_imNoisy.sz.height ||
		 i_fflow.sz.frames != i_imNoisy.sz.frames))
		throw std::runtime_error("VideoNLB::runNlBayes: Wrong size of fwd optical flow.");

	if (use_oflow && 
		(i_bflow.sz.channels != 2 || 
		 i_bflow.sz.width  != i_imNoisy.sz.width  ||
		 i_bflow.sz.height != i_imNoisy.sz.height ||
		 i_bflow.sz.frames != i_imNoisy.sz.frames))
		throw std::runtime_error("VideoNLB::runNlBayes: Wrong size of bwd optical flow.");

	//! Print compiler options
	if (p_prms1.verbose)
	{
#ifdef DCT_BASIS
		printf(ANSI_BCYN "DCT_BASIS > Using DCT basis\n" ANSI_RST);
#endif
#ifdef USE_FFTW 
		printf(ANSI_BCYN "USE_FFTW > \n" ANSI_RST);
#endif
#ifdef DCT_DONT_CENTER1
		printf(ANSI_BCYN "DCT_DONT_CENTER1 > Centering DCT step 1\n" ANSI_RST);
#endif
#ifdef DCT_DONT_CENTER2
		printf(ANSI_BCYN "DCT_DONT_CENTER2 > Centering DCT step 2\n" ANSI_RST);
#endif
#ifdef NOISY_COVARIANCE2
		printf(ANSI_BCYN "NOISY_COVARIANCE2 > Computing 2nd step cov. matrix from noisy patches.\n" ANSI_RST);
#endif
#ifdef LI_ZHANG_DAI1
		printf(ANSI_BCYN "LI_ZHANG_DAI1 > Using Li-Zhang-Dai's empirical Wiener, step 1.\n" ANSI_RST);
#endif
#if defined(LI_ZHANG_DAI2) && defined(NOISY_COVARIANCE2)
		printf(ANSI_BCYN "LI_ZHANG_DAI2 > Using Li-Zhang-Dai's empirical Wiener, step 2.\n" ANSI_RST);
#endif
#ifdef BARICENTER_BASIC
		printf(ANSI_BCYN "BARICENTER_BASIC > Centering noisy patches with basic baricenter.\n" ANSI_RST);
#endif
#if defined(THRESHOLD_WEIGHTS1) && !(defined(LINEAR_THRESHOLDING1) || defined(THRESHOLDING1))
		printf(ANSI_BCYN "THRESHOLD_WEIGHTS1 > Thresholding step 1 negative Wiener weights\n" ANSI_RST);
#endif
#ifdef THRESHOLDING1
		printf(ANSI_BCYN "THRESHOLDING1 > Coefficient thresholding step 1 instead of Wiener weights\n" ANSI_RST);
#endif
#ifdef USE_BETA_FOR_VARIANCE
		printf(ANSI_BCYN "USE_BETA_FOR_VARIANCE > Noise correction in step 1 applied both MAP and variances.\n" ANSI_RST);
#endif
#if defined(THRESHOLDING1) && defined(SOFT_THRESHOLD1)
		printf(ANSI_BCYN "SOFT_THRESHOLD1 > Coefficient thresholding step 1 is soft\n" ANSI_RST);
#endif
#if defined(THRESHOLDING1) && defined(SOFT_THRESHOLD1_BAYES)
		printf(ANSI_BCYN "SOFT_THRESHOLD1_BAYES > Coefficient thresholding step 1 is soft, Bayesian\n" ANSI_RST);
#endif
#if defined(LINEAR_THRESHOLDING1) && !defined(THRESHOLDING1) 
		printf(ANSI_BCYN "LINEAR_THRESHOLDING1 > Variances thresholding step 1 instead of Wiener weights\n" ANSI_RST);
#endif
#ifdef THRESHOLDING2
		printf(ANSI_BCYN "THRESHOLDING2 > Coefficient thresholding step 2 instead of Wiener weights\n" ANSI_RST);
#endif
#if defined(THRESHOLDING2) && defined(SOFT_THRESHOLD2)
		printf(ANSI_BCYN "SOFT_THRESHOLD2 > Coefficient thresholding step 2 is soft\n" ANSI_RST);
#endif
#if defined(THRESHOLDING2) && defined(SOFT_THRESHOLD2_BAYES)
		printf(ANSI_BCYN "SOFT_THRESHOLD2_BAYES > Coefficient thresholding step 2 is soft, Bayesian\n" ANSI_RST);
#endif
#if defined(LINEAR_SOFT_THRESHOLDING2) && !defined(THRESHOLDING2) 
		printf(ANSI_BCYN "LINEAR_SOFT_THRESHOLDING2 > Variances soft thresholding step 2 instead of Wiener weights\n" ANSI_RST);
#endif
#if defined(LINEAR_HARD_THRESHOLDING2) && !defined(THRESHOLDING2) 
		printf(ANSI_BCYN "LINEAR_HARD_THRESHOLDING2 > Variances hard thresholding step 2 instead of Wiener weights\n" ANSI_RST);
#endif
#if defined(THRESHOLD_WEIGHTS2) && defined(NOISY_COVARIANCE2)
		printf(ANSI_BCYN "THRESHOLD_WEIGHTS2 > Thresholding step 2 negative Wiener weights\n" ANSI_RST);
#endif
#if defined(MEAN_HYPERPRIOR1) && !defined(MEAN_HYPERPRIOR_BM3D1)
		printf(ANSI_BCYN "MEAN_HYPERPRIOR1 > Assuming a Bayesian prior over the sample mean, step 1\n" ANSI_RST);
#endif
#if defined(MEAN_HYPERPRIOR1) && defined(MEAN_HYPERPRIOR_BM3D1)
		printf(ANSI_BCYN "MEAN_HYPERPRIOR_BM3D1 > Assuming a BM3D prior over the sample mean, step 1\n" ANSI_RST);
#endif
#if defined(MEAN_HYPERPRIOR2) && !defined(MEAN_HYPERPRIOR_BM3D2)
		printf(ANSI_BCYN "MEAN_HYPERPRIOR2 > Assuming a Bayesian prior over the sample mean, step 2\n" ANSI_RST);
#endif
#if defined(MEAN_HYPERPRIOR2) && defined(MEAN_HYPERPRIOR_BM3D2)
		printf(ANSI_BCYN "MEAN_HYPERPRIOR_BM3D2 > Assuming a BM3D prior over the sample mean, step 2\n" ANSI_RST);
#endif
	}


	//! Number of available cores
	unsigned nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
	if (p_prms1.verbose) printf(ANSI_CYN "OpenMP is using %d threads\n" ANSI_RST, nThreads);
#endif
	const unsigned nParts = 2 * nThreads;

	//! Video size
	VideoSize imSize = i_imNoisy.sz;

	//! Determine steps
	int step = p_prms1.sizePatch ?
	          (p_prms2.sizePatch ? 0 : 1) :
	          (p_prms2.sizePatch ? 2 :-1) ;

	if (step == -1)
		throw std::runtime_error("VideoNLB::runNlBayes: Both patch sizes are zero.");

	//! Initialization
	if (step != 2) o_imBasic.resize(imSize);
	if (step != 1) o_imFinal.resize(imSize);

	if (step == 2 && imSize != o_imBasic.sz) 
		throw std::runtime_error("VideoNLB::runNlBayes: sizes of noisy and "
				"basic videos don't match");

	//! Print parameters
	if (p_prms1.verbose) if (step != 2) printNlbParameters(p_prms1);
	if (p_prms2.verbose) if (step != 1) printNlbParameters(p_prms2);

	//! Percentage of processed groups over total number of pixels
	std::vector<float> groupsRatio(2,0.f);

	//! Step 1
	if (step != 2)
	{
		if (p_prms1.verbose)
		{
			printf("1st Step\n");
			for (int p = 0; p < nParts; ++p) printf("\n");
		}

		//! RGB to YUV
		Video<float> imNoisy = i_imNoisy;
		if (p_prms1.colorSpace == YUV)
			VideoUtils::transformColorSpace(imNoisy, true);

		//! Divide the noisy image into sub-images in order to easier parallelize the process
		std::vector<Video<float> > imNoisySub(nParts);
		std::vector<VideoUtils::CropPosition > imCrops(nParts);
		VideoUtils::subDivideTight(imNoisy, imNoisySub, imCrops, p_prms1.boundary, nParts);

		//! Divide optical flow 
		std::vector<Video<float> > fflowSub(nParts), bflowSub(nParts);
		if (use_oflow)
		{
			std::vector<VideoUtils::CropPosition > oflowCrops(nParts);
			VideoUtils::subDivideTight(i_fflow, fflowSub, oflowCrops, p_prms1.boundary, nParts);
			VideoUtils::subDivideTight(i_bflow, bflowSub, oflowCrops, p_prms1.boundary, nParts);
		}

#ifdef DEBUG_COMPUTE_GROUP_ERROR
		if (p_prms1.colorSpace == YUV)
			VideoUtils::transformColorSpace(i_imClean, true);
		std::vector<Video<float> > imCleanSub(nParts);
		VideoUtils::subDivideTight(i_imClean, imCleanSub, imCrops, p_prms1.boundary, nParts);
#endif

#ifdef USE_FFTW
		//! Initialize DCT algorithms
		globalDCT.init(p_prms1.sizePatch,
		               p_prms1.sizePatch,
		               p_prms1.sizePatchTime,
		               p_prms1.nSimilarPatches,
		               nThreads);
#endif

		//! Process all sub-images
		std::vector<Video<float> > imBasicSub(nParts);
		std::vector<Video<float> > imFinalSub(nParts);
		std::vector<unsigned> groupsProcessedSub(nParts);
#ifdef _OPENMP
		// we make a copy of prms structure because, since it is constant,
		// it causes a compilation error with OpenMP (only on IPOL server)
		nlbParams prms1(p_prms1);
#pragma omp parallel for schedule(dynamic, nParts/nThreads) \
		shared(imNoisySub, imBasicSub, imFinalSub) \
		firstprivate (prms1)
#endif
		for (int n = 0; n < (int)nParts; n++)
			groupsProcessedSub[n] = 
#ifndef DEBUG_COMPUTE_GROUP_ERROR
				processNlBayes(imNoisySub[n], fflowSub[n], bflowSub[n],
				               imBasicSub[n], imFinalSub[n], p_prms1, imCrops[n]);
#else
				processNlBayes(imNoisySub[n], fflowSub[n], bflowSub[n],
				               imBasicSub[n], imFinalSub[n], imCleanSub[n],
				               p_prms1, imCrops[n]);
#endif

#ifdef USE_FFTW
		//! Destroy DCT algorithms
		globalDCT.destroy();
#endif

		//! Get the basic estimate
		VideoUtils::subBuildTight(imBasicSub, o_imBasic, p_prms1.boundary);

		//! YUV to RGB
		if (p_prms1.colorSpace == YUV)
			VideoUtils::transformColorSpace(o_imBasic, false);

#ifdef DEBUG_COMPUTE_GROUP_ERROR
		if (p_prms1.colorSpace == YUV)
			VideoUtils::transformColorSpace(i_imClean, false);
#endif

		for (int n = 0; n < (int)nParts; n++)
			groupsRatio[0] += 100.f * (float)groupsProcessedSub[n]/(float)imSize.whf;

#ifdef DEBUG_SHOW_WEIGHT
		{
			std::vector<Video<float> > subWeights(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].tile_x;
				int part_y = imCrops[n].tile_y;
				int part_t = imCrops[n].tile_t;
				char name[1024];
				sprintf(name, "dump/weight_step1_%d.%d.%d_%%03d.tif", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subWeights[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> weight(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subWeights, weight, p_prms1.boundary);

			// Write to disk
			weight.saveVideo("wei1_%03d.tif", 1);
		}
		{
			std::vector<Video<float> > subVars(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].tile_x;
				int part_y = imCrops[n].tile_y;
				int part_t = imCrops[n].tile_t;
				char name[1024];
				sprintf(name, "dump/var_step1_%d.%d.%d_%%03d.tif", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subVars[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> variance(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subVars, variance, p_prms1.boundary);

			// Write to disk
			variance.saveVideo("var1_%03d.tif", 1);
		}
 #ifdef DEBUG_COMPUTE_GROUP_ERROR
		{
			std::vector<Video<float> > subVars(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].tile_x;
				int part_y = imCrops[n].tile_y;
				int part_t = imCrops[n].tile_t;
				char name[1024];
				sprintf(name, "dump/pge_step1_%d.%d.%d_%%03d.tif", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subVars[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> variance(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subVars, variance, p_prms1.boundary);

			// Write to disk
			variance.saveVideo("pge1_%03d.tif", 1);
		}
 #endif
#endif

	}

	//! Step 2
	if (step != 1)
	{
		if (p_prms2.verbose)
		{
			if (step == 2) for (int p = 0; p <= nParts; ++p) printf("\n");
			printf("\x1b[%dF2nd Step\n",nParts+1);
			for (int p = 0; p < nParts; ++p) printf("\x1b[2K\n");
		}

		//! RGB to YUV
		Video<float> imNoisy = i_imNoisy;
		if (p_prms2.colorSpace == YUV)
		{
			VideoUtils::transformColorSpace(  imNoisy, true);
			VideoUtils::transformColorSpace(o_imBasic, true);
		}

		//! Divide the noisy and basic images into sub-images in order to easier parallelize the process
		std::vector<Video<float> > imNoisySub(nParts);
		std::vector<Video<float> > imBasicSub(nParts);
		std::vector<VideoUtils::CropPosition > imCrops(nParts);

		VideoUtils::subDivideTight(  imNoisy, imNoisySub, imCrops, p_prms2.boundary, nParts);
		VideoUtils::subDivideTight(o_imBasic, imBasicSub, imCrops, p_prms2.boundary, nParts);

		//! Divide optical flow 
		std::vector<Video<float> > fflowSub(nParts), bflowSub(nParts);
		if (use_oflow)
		{
			std::vector<VideoUtils::CropPosition > oflowCrops(nParts);
			VideoUtils::subDivideTight(i_fflow, fflowSub, oflowCrops, p_prms2.boundary, nParts);
			VideoUtils::subDivideTight(i_bflow, bflowSub, oflowCrops, p_prms2.boundary, nParts);
		}

#ifdef DEBUG_COMPUTE_GROUP_ERROR
		if (p_prms2.colorSpace == YUV)
			VideoUtils::transformColorSpace(i_imClean, true);
		std::vector<Video<float> > imCleanSub(nParts);
		VideoUtils::subDivideTight(i_imClean, imCleanSub, imCrops, p_prms2.boundary, nParts);
#endif

#ifdef USE_FFTW
		//! Initialize DCT algorithms
		globalDCT.init(p_prms2.sizePatch,
		               p_prms2.sizePatch,
		               p_prms2.sizePatchTime,
		               p_prms2.nSimilarPatches,
		               nThreads);
#endif

		//! Process all sub-images
		std::vector<Video<float> > imFinalSub(nParts);
		std::vector<float> groupsProcessedSub(nParts);
#ifdef _OPENMP
		// we make a copy of prms structure because, since it is constant,
		// it causes a compilation error with OpenMP (only on IPOL server)
		nlbParams prms2(p_prms2);
#pragma omp parallel for schedule(dynamic, nParts/nThreads) \
		shared(imNoisySub, imBasicSub, imFinalSub) \
		firstprivate (prms2)
#endif
		for (int n = 0; n < (int) nParts; n++)
		{
			groupsProcessedSub[n] = 
#ifndef DEBUG_COMPUTE_GROUP_ERROR
				processNlBayes(imNoisySub[n], fflowSub[n], bflowSub[n],
				               imBasicSub[n], imFinalSub[n], p_prms2, imCrops[n]);
#else
				processNlBayes(imNoisySub[n], fflowSub[n], bflowSub[n],
				               imBasicSub[n], imFinalSub[n], imCleanSub[n],
				               p_prms2, imCrops[n]);
#endif
		}

#ifdef USE_FFTW
		//! Destroy DCT algorithms
		globalDCT.destroy();
#endif

		//! Get the final result
		VideoUtils::subBuildTight(imFinalSub, o_imFinal, p_prms2.boundary);

		//! Undo color transform
		if (p_prms2.colorSpace == YUV)
		{
			VideoUtils::transformColorSpace(o_imBasic, false);
			VideoUtils::transformColorSpace(o_imFinal, false);
		}

#ifdef DEBUG_COMPUTE_GROUP_ERROR
		if (p_prms2.colorSpace == YUV)
			VideoUtils::transformColorSpace(i_imClean, false);
#endif

		for (int n = 0; n < (int)nParts; n++)
			groupsRatio[1] += 100.f * (float)groupsProcessedSub[n]/(float)imSize.whf;

#ifdef DEBUG_SHOW_WEIGHT
		{
			std::vector<Video<float> > subWeights(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].tile_x;
				int part_y = imCrops[n].tile_y;
				int part_t = imCrops[n].tile_t;
				char name[1024];
				sprintf(name, "dump/weight_step2_%d.%d.%d_%%03d.tif", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subWeights[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> weight(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subWeights, weight, p_prms2.boundary);

			// Write to disk
			weight.saveVideo("wei2_%03d.tif", 1);
		}
		{
			std::vector<Video<float> > subVars(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].tile_x;
				int part_y = imCrops[n].tile_y;
				int part_t = imCrops[n].tile_t;
				char name[1024];
				sprintf(name, "dump/var_step2_%d.%d.%d_%%03d.tif", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subVars[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> variance(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subVars, variance, p_prms2.boundary);

			// Write to disk
			variance.saveVideo("var2_%03d.tif", 1);
		}
 #ifdef DEBUG_COMPUTE_GROUP_ERROR
		{
			std::vector<Video<float> > subVars(nParts);

			// First load all weight sequences
			for (int n = 0; n < (int)nParts; n++)
			{
				// Build file name
				int part_x = imCrops[n].tile_x;
				int part_y = imCrops[n].tile_y;
				int part_t = imCrops[n].tile_t;
				char name[1024];
				sprintf(name, "dump/pge_step2_%d.%d.%d_%%03d.tif", part_x, part_y, part_t);

				int part_frames = imCrops[n].ending_t - part_t;
				subVars[n].loadVideo(name, 1, part_frames, 1);
			}

			// Call set build
			Video<float> variance(imSize.width, imSize.height, imSize.frames);
			VideoUtils::subBuildTight(subVars, variance, p_prms2.boundary);

			// Write to disk
			variance.saveVideo("pge2_%03d.tif", 1);
		}
 #endif
#endif
		if (p_prms2.verbose) printf("\n");
	}

	return groupsRatio;
}

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy image;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_params: see nlbParams.
 *
 * @return Number of processed groups of similar patches.
 **/
unsigned processNlBayes(
	Video<float> const& i_imNoisy
,	Video<float> const& i_fflow
,	Video<float> const& i_bflow
,	Video<float> &io_imBasic
,	Video<float> &o_imFinal
#ifdef DEBUG_COMPUTE_GROUP_ERROR
,	Video<float> const& i_imClean
#endif
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
	const VideoSize sz = i_imNoisy.sz;

	unsigned nInverseFailed = 0;
	const float threshold = p_params.sigma * p_params.sigma * p_params.gamma *
	                       (p_params.isFirstStep ? i_imNoisy.sz.channels : 1.f);

	//! Weight sum per pixel
	Video<float> weight(sz.width, sz.height, sz.frames, 1, 0.f);

	//! Mask: true for pixels that still need to be processed
	Video<char> mask(sz.width, sz.height, sz.frames, 1, false);

	//! There's a border added only if the crop
	//  doesn't touch the source image boundary
	bool border_x0 = p_crop.origin_x > 0;
	bool border_y0 = p_crop.origin_y > 0;
	bool border_t0 = p_crop.origin_t > 0;
	bool border_x1 = p_crop.ending_x < p_crop.source_sz.width;
	bool border_y1 = p_crop.ending_y < p_crop.source_sz.height;
	bool border_t1 = p_crop.ending_t < p_crop.source_sz.frames;

	//! Only pixels of the center of the image must be processed (not the boundaries)
	int n_groups = 0;
	unsigned stepx = p_params.offSet;
	unsigned stepy = p_params.offSet;
	unsigned stepf = p_params.offSetTime;
	int ori_x =                        border_x0 ? sPx-1 + sWx/2 : 0 ;
	int ori_y =                        border_y0 ? sPx-1 + sWx/2 : 0 ;
	int ori_f =                        border_t0 ? sPt-1 + sWt/2 : 0 ;
	int end_x = (int)sz.width  - (int)(border_x1 ? sPx-1 + sWx/2 : sPx-1);
	int end_y = (int)sz.height - (int)(border_y1 ? sPx-1 + sWx/2 : sPx-1);
	int end_f = (int)sz.frames - (int)(border_t1 ? sPt-1 + sWt/2 : sPt-1);
	for (int f = ori_f, df = 0; f < end_f; f++, df++)
	for (int y = ori_y, dy = 0; y < end_y; y++, dy++)
	for (int x = ori_x, dx = 0; x < end_x; x++, dx++)
	{
		if ( (df % stepf == 0) || (!border_t1 && f == end_f - 1))
		{
			int phasey = (!border_t1 && f == end_f - 1) ? 0 : f/stepf;

			if ( (dy % stepy == phasey % stepy) ||
			     (!border_y1 && y == end_y - 1) ||
				  (!border_y0 && y == ori_y    ) )
			{
				int phasex = (!border_y1 && y == end_y - 1) ? 0 : (phasey + y/stepy);

				if ( (dx % stepx == phasex % stepx) ||
				     (!border_x1 && x == end_x - 1) ||
				     (!border_x0 && x == ori_x    ) )
				{
					mask(x,y,f) = true;
					n_groups++;
				}
			}
		}
	}

//	printf("Processing at most %d groups of similar patches\n", n_groups);


#ifdef DEBUG_SHOW_WEIGHT
	{
		int part_x = p_crop.tile_x;
		int part_y = p_crop.tile_y;
		int part_t = p_crop.tile_t;
		char name[1024];
		Video<float> mask_f(mask.sz);
		for (int i = 0; i < mask.sz.whcf; ++i) mask_f(i) = 255*(float)mask(i);
		sprintf(name, "dump/msk_step%d_%d.%d.%d_%%03d.png", step1 ? 1 : 2, part_x, part_y, part_t);
		mask_f.saveVideo(name, 1, 1);
	}
#endif

	//! Used matrices during Bayes' estimate
	const unsigned patch_dim = step1 ? sPx * sPx * sPt : sPx * sPx * sPt * sz.channels;
//	const unsigned patch_num = step1 ? p_params.nSimilarPatches : sWx * sWx * sWt;
	const unsigned patch_num = sWx * sWx * sWt;

	//! Matrices used for Bayes' estimate
	vector<unsigned> index(patch_num);
	matWorkspace mat;
	mat.groupTranspose     .resize(patch_num * patch_dim);
	mat.tmpMat             .resize(patch_dim * patch_dim);
	mat.covMat             .resize(patch_dim * patch_dim);
	mat.covMatTmp          .resize(patch_dim * patch_dim);
	mat.baricenter         .resize(patch_dim);
	mat.groupTransposeNoisy.resize(patch_num * patch_dim); //FIXME
	mat.baricenterNoisy    .resize(patch_dim); //FIXME

	//! Store DCT basis
	const unsigned sPc = 1; // number of channels TODO: obsolete!
	mat.patch_basis    .resize(patch_dim * patch_dim);
	mat.patch_basis_inv.resize(patch_dim * patch_dim);
	mat.patch_basis_x.resize(sPx * sPc * sPx * sPc);
	mat.patch_basis_y.resize(sPx * sPc * sPx * sPc);
	mat.patch_basis_t.resize(sPt * sPc * sPt * sPc);
	if (p_params.transform == dct)
	{
		// 1D DCT basis for signals of length sPx
		std::vector<float> cosx(sPx * sPx);
		for (unsigned kx = 0, i = 0; kx < sPx; kx++)
		for (unsigned nx = 0       ; nx < sPx; nx++, i++)
			cosx[i] = cos( (0.5 + (float)nx)*  (float)kx * M_PI / (float)sPx )
			        * sqrt( (kx == 0 ? 1. : 2.) / (float) sPx );

		// 1D DCT basis for signals of length sPt
		std::vector<float> cost(sPt * sPt);
		for (unsigned kt = 0, i = 0; kt < sPt; kt++)
		for (unsigned nt = 0       ; nt < sPt; nt++, i++)
//			cost[i] = nt == kt ? 1 : 0;
			cost[i] = cos( (0.5 + (float)nt)*  (float)kt * M_PI / (float)sPt )
			        * sqrt( (kt == 0 ? 1. : 2.) / (float) sPt );

		//! 3D DCT basis
		for (unsigned kc = 0, i = 0; kc < sPc; kc++)
		for (unsigned kt = 0       ; kt < sPt; kt++) // triplets (kt,ky,kx) are
		for (unsigned ky = 0       ; ky < sPx; ky++) // frequencies and index
		for (unsigned kx = 0       ; kx < sPx; kx++) // basis vectors
			for (unsigned nc = 0    ; nc < sPc; nc++)
			for (unsigned nt = 0    ; nt < sPt; nt++)
			for (unsigned ny = 0    ; ny < sPx; ny++)      // triplets (nt,ny,nx) are
			for (unsigned nx = 0    ; nx < sPx; nx++, i++) // patch positions
				mat.patch_basis[i] = (nc == kc) ? cosx[kx * sPx + nx]
				                                * cosx[ky * sPx + ny]
				                                * cost[kt * sPt + nt] : 0;

		//! inverse transform
		mat.patch_basis_inv = mat.patch_basis;

		//! 1D DCT bases
		for (unsigned kc = 0, i = 0; kc < sPc; kc++)
		for (unsigned kx = 0       ; kx < sPx; kx++) // basis vectors
			for (unsigned nc = 0    ; nc < sPc; nc++)
			for (unsigned nx = 0    ; nx < sPx; nx++, i++) // patch positions
				mat.patch_basis_x[i] = (nc == kc) ? cosx[kx * sPx + nx] : 0;

		for (unsigned kc = 0, i = 0; kc < sPc; kc++)
		for (unsigned ky = 0       ; ky < sPx; ky++) // basis vectors
			for (unsigned nc = 0    ; nc < sPc; nc++)
			for (unsigned ny = 0    ; ny < sPx; ny++, i++) // patch positions
				mat.patch_basis_y[i] = (nc == kc) ? cosx[ky * sPx + ny] : 0;

		for (unsigned kc = 0, i = 0; kc < sPc; kc++)
		for (unsigned kt = 0       ; kt < sPt; kt++) // basis vectors
			for (unsigned nc = 0    ; nc < sPc; nc++)
			for (unsigned nt = 0    ; nt < sPt; nt++, i++) // patch positions
				mat.patch_basis_t[i] = (nc == kc) ? cost[kt * sPt + nt] : 0;
	
		/*! Verify
		for (unsigned kx = 0, i = 0; kx < sPx; kx++)
		{
			float norm = 0.;
			for (unsigned nx = 0; nx < sPx; nx++, i++) norm += cosx[i]*cosx[i];
			printf("|cosx[%d]| = %f\n",kx,norm);
		}

		for (unsigned kt = 0, i = 0; kt < sPt; kt++)
		{
			float norm = 0.;
			for (unsigned nt = 0; nt < sPt; nt++, i++) norm += cost[i]*cost[i];
			printf("|cost[%d]| = %f\n",kt,norm);
		}

		for (unsigned kc = 0, i = 0; kc < sPc; kc++)
		for (unsigned kt = 0       ; kt < sPt; kt++) // triplets (kt,ky,kx) are
		for (unsigned ky = 0       ; ky < sPx; ky++) // frequencies and index
		for (unsigned kx = 0       ; kx < sPx; kx++) // basis vectors
		{
			float norm = 0.;
			for (unsigned nc = 0; nc < sPc; nc++)
			for (unsigned nt = 0; nt < sPt; nt++)
			for (unsigned ny = 0; ny < sPx; ny++)      // triplets (nt,ny,nx) are
			for (unsigned nx = 0; nx < sPx; nx++, i++) // patch positions
				norm += mat.covEigVecs[i]*mat.covEigVecs[i];

			printf("|v[%d,%d,%d,%d]| = %f\n",kc,kt,ky,kx,norm);
		}*/
	}
	else
	{
		if (sPx != 8)
			throw std::runtime_error("VideoNLB::processNlBayes: bior1_5 transform is only "
					"implemented for patches of size 8 x 8 x pt.");

		// 1D bior1.5 basis for signals of length sPx=8
		std::vector<float> bx(sPx * sPx);
		int i = 0;
		bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274;
		bx[i++]= 0.219417649252501; bx[i++]= 0.449283757993216; bx[i++]= 0.449283757993216; bx[i++]= 0.219417649252501; bx[i++]=-0.219417649252501; bx[i++]=-0.449283757993216; bx[i++]=-0.449283757993216; bx[i++]=-0.219417649252501;
		bx[i++]= 0.569359398342846; bx[i++]= 0.402347308162278; bx[i++]=-0.402347308162278; bx[i++]=-0.569359398342846; bx[i++]=-0.083506045090284; bx[i++]= 0.083506045090284; bx[i++]=-0.083506045090284; bx[i++]= 0.083506045090284;
		bx[i++]=-0.083506045090284; bx[i++]= 0.083506045090284; bx[i++]=-0.083506045090284; bx[i++]= 0.083506045090284; bx[i++]= 0.569359398342846; bx[i++]= 0.402347308162278; bx[i++]=-0.402347308162278; bx[i++]=-0.569359398342846;
		bx[i++]= 0.707106781186547; bx[i++]=-0.707106781186547; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0;
		bx[i++]=                 0; bx[i++]=                 0; bx[i++]= 0.707106781186547; bx[i++]=-0.707106781186547; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0;
		bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]= 0.707106781186547; bx[i++]=-0.707106781186547; bx[i++]=                 0; bx[i++]=                 0;
		bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]=                 0; bx[i++]= 0.707106781186547; bx[i++]=-0.707106781186547;

		// 1D basis for signals of length sPt (identity matrix)
		std::vector<float> bt(sPt * sPt);
		for (unsigned kt = 0, i = 0; kt < sPt; kt++)
		for (unsigned nt = 0       ; nt < sPt; nt++, i++)
			bt[i] = nt == kt ? 1 : 0;

		//! 3D DCT basis
		for (unsigned kc = 0, i = 0; kc < sPc; kc++)
		for (unsigned kt = 0       ; kt < sPt; kt++) // triplets (kt,ky,kx) are
		for (unsigned ky = 0       ; ky < sPx; ky++) // frequencies and index
		for (unsigned kx = 0       ; kx < sPx; kx++) // basis vectors
			for (unsigned nc = 0    ; nc < sPc; nc++)
			for (unsigned nt = 0    ; nt < sPt; nt++)
			for (unsigned ny = 0    ; ny < sPx; ny++)      // triplets (nt,ny,nx) are
			for (unsigned nx = 0    ; nx < sPx; nx++, i++) // patch positions
				mat.patch_basis[i] = (nc == kc) ? bx[kx * sPx + nx]
				                                * bx[ky * sPx + ny]
				                                * bt[kt * sPt + nt] : 0;

		// same for the inverse transform
		i = 0;
		bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]= 0.5; bx[i++]= 0.0; bx[i++]= 0.707106781186547; bx[i++]=-0.121533978016438; bx[i++]=-0.000000000000000; bx[i++]= 0.121533978016438;
		bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]= 0.5; bx[i++]=-0.0; bx[i++]=-0.707106781186547; bx[i++]=-0.121533978016438; bx[i++]= 0.000000000000000; bx[i++]= 0.121533978016438;
		bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]=-0.5; bx[i++]=-0.0; bx[i++]= 0.121533978016438; bx[i++]= 0.707106781186548; bx[i++]=-0.121533978016438; bx[i++]= 0.000000000000000;
		bx[i++]= 0.353553390593274; bx[i++]= 0.353553390593274; bx[i++]=-0.5; bx[i++]= 0.0; bx[i++]= 0.121533978016438; bx[i++]=-0.707106781186548; bx[i++]=-0.121533978016438; bx[i++]= 0.000000000000000;
		bx[i++]= 0.353553390593274; bx[i++]=-0.353553390593274; bx[i++]= 0.0; bx[i++]= 0.5; bx[i++]= 0.000000000000000; bx[i++]= 0.121533978016438; bx[i++]= 0.707106781186547; bx[i++]=-0.121533978016438;
		bx[i++]= 0.353553390593274; bx[i++]=-0.353553390593274; bx[i++]= 0.0; bx[i++]= 0.5; bx[i++]= 0.000000000000000; bx[i++]= 0.121533978016438; bx[i++]=-0.707106781186547; bx[i++]=-0.121533978016438;
		bx[i++]= 0.353553390593274; bx[i++]=-0.353553390593274; bx[i++]= 0.0; bx[i++]=-0.5; bx[i++]=-0.121533978016438; bx[i++]=-0.000000000000000; bx[i++]= 0.121533978016438; bx[i++]= 0.707106781186547;
		bx[i++]= 0.353553390593274; bx[i++]=-0.353553390593274; bx[i++]= 0.0; bx[i++]=-0.5; bx[i++]=-0.121533978016438; bx[i++]=-0.000000000000000; bx[i++]= 0.121533978016438; bx[i++]=-0.707106781186547;

		// 1D basis for signals of length sPt (identity matrix)
		for (unsigned kt = 0, i = 0; kt < sPt; kt++)
		for (unsigned nt = 0       ; nt < sPt; nt++, i++)
			bt[i] = nt == kt ? 1 : 0;

		//! 3D DCT basis (stored transposed, so we can use the dct code)
		for (unsigned nc = 0, i = 0; nc < sPc; nc++)
		for (unsigned nt = 0       ; nt < sPt; nt++) // triplets (kt,ky,kx) are
		for (unsigned ny = 0       ; ny < sPx; ny++) // frequencies and index
		for (unsigned nx = 0       ; nx < sPx; nx++) // basis vectors
			for (unsigned kc = 0    ; kc < sPc; kc++)
			for (unsigned kt = 0    ; kt < sPt; kt++)
			for (unsigned ky = 0    ; ky < sPx; ky++)      // triplets (nt,ny,nx) are
			for (unsigned kx = 0    ; kx < sPx; kx++, i++) // patch positions
				mat.patch_basis_inv[i] = (nc == kc) ? bx[kx * sPx + nx]
				                                    * bx[ky * sPx + ny]
				                                    * bt[kt * sPt + nt] : 0;
	}

	mat.agg_window.resize(patch_dim);
	{
		if (sPx == 8) // use BM3D's Kaiser window
		{
			float *z = mat.agg_window.data();

			for (unsigned t = 0; t < sPt; t++)
			{
				*z++ = 0.1924; *z++ = 0.2989; *z++ = 0.3846; *z++ = 0.4325; *z++ = 0.4325; *z++ = 0.3846; *z++ = 0.2989; *z++ = 0.1924;
				*z++ = 0.2989; *z++ = 0.4642; *z++ = 0.5974; *z++ = 0.6717; *z++ = 0.6717; *z++ = 0.5974; *z++ = 0.4642; *z++ = 0.2989;
				*z++ = 0.3846; *z++ = 0.5974; *z++ = 0.7688; *z++ = 0.8644; *z++ = 0.8644; *z++ = 0.7688; *z++ = 0.5974; *z++ = 0.3846;
				*z++ = 0.4325; *z++ = 0.6717; *z++ = 0.8644; *z++ = 0.9718; *z++ = 0.9718; *z++ = 0.8644; *z++ = 0.6717; *z++ = 0.4325;
				*z++ = 0.4325; *z++ = 0.6717; *z++ = 0.8644; *z++ = 0.9718; *z++ = 0.9718; *z++ = 0.8644; *z++ = 0.6717; *z++ = 0.4325;
				*z++ = 0.3846; *z++ = 0.5974; *z++ = 0.7688; *z++ = 0.8644; *z++ = 0.8644; *z++ = 0.7688; *z++ = 0.5974; *z++ = 0.3846;
				*z++ = 0.2989; *z++ = 0.4642; *z++ = 0.5974; *z++ = 0.6717; *z++ = 0.6717; *z++ = 0.5974; *z++ = 0.4642; *z++ = 0.2989;
				*z++ = 0.1924; *z++ = 0.2989; *z++ = 0.3846; *z++ = 0.4325; *z++ = 0.4325; *z++ = 0.3846; *z++ = 0.2989; *z++ = 0.1924;
			}
		}
		else
			for (unsigned i = 0; i < sPt*sPx*sPx; i++)
				mat.agg_window[i] = 1.;
	}


	//! Variance captured by the principal components
	Video<float> variance(mask.sz);
#ifdef DEBUG_COMPUTE_GROUP_ERROR
	Video<float> group_error(mask.sz);
#endif


	//! Total number of groups of similar patches processed
	unsigned group_counter = 0;

	if (step1)
	{
		//! Allocate Sizes
		io_imBasic.resize(sz);

		//! Matrices used for Bayes' estimate
		vector<vector<float> > group(sz.channels, vector<float>(patch_num * patch_dim));
		vector<vector<float> > aggreWeights(sz.channels, vector<float>(patch_num, 1.f));
#ifdef DEBUG_COMPUTE_GROUP_ERROR
		vector<vector<float> > groupClean(sz.channels, vector<float>(patch_num * patch_dim));
#endif

		int remaining_groups = n_groups;
		for (unsigned pt = 0; pt < sz.frames; pt++)
		for (unsigned py = 0; py < sz.height; py++)
		for (unsigned px = 0; px < sz.width ; px++)
			if (mask(px,py,pt)) //< Only non-seen patches are processed
			{
				group_counter++;

				const unsigned ij  = sz.index(px,py,pt);
				const unsigned ij3 = sz.index(px,py,pt, 0);
				//const unsigned ij3 = (ij / sz.wh) * sz.whc + ij % sz.wh;

				if (p_params.verbose && (group_counter % 100 == 0))
				{
					int ntiles = p_crop.ntiles_t * p_crop.ntiles_x * p_crop.ntiles_y;
					int part_idx = p_crop.tile_t * p_crop.ntiles_x * p_crop.ntiles_y +
					               p_crop.tile_y * p_crop.ntiles_x + 
					               p_crop.tile_x;

					printf("\x1b[%dF[%d,%d,%d] %05.1f\x1b[%dE", ntiles - part_idx,
							p_crop.tile_x, p_crop.tile_y, p_crop.tile_t,
							100.f - (float)remaining_groups/(float)(n_groups)*100.f,
							ntiles - part_idx);

					std::cout << std::flush;
				}

				//! Search for similar patches around the reference one
#ifndef DEBUG_COMPUTE_GROUP_ERROR
				unsigned nSimP = estimateSimilarPatchesStep1(i_imNoisy, i_fflow, i_bflow,
						group, index, ij3, p_params);
#else
				unsigned nSimP = estimateSimilarPatchesStep1(i_imNoisy, i_fflow, i_bflow,
						group, index, ij3, p_params, i_imClean, groupClean);
#endif

				//! If we use the homogeneous area trick
				bool doBayesEstimate = true;
				if (p_params.useHomogeneousArea)
					doBayesEstimate = !computeHomogeneousAreaStep1(group, sPx,
							patch_num, threshold, sz);

				//! Else, use Bayes' estimate
				if (doBayesEstimate)
					variance(ij) = computeBayesEstimateStep1(group, mat,
							nInverseFailed, p_params, nSimP, aggreWeights);

#ifdef DEBUG_COMPUTE_GROUP_ERROR
				{
					float groupError = 0, tmp;
					for (int c = 0; c < sz.channels; ++c)
					for (int i = 0; i < nSimP * patch_dim; ++i)
						groupError += (tmp = group[c][i] - groupClean[c][i])*tmp;

					groupError /= (float)(sz.channels*nSimP*patch_dim);
					group_error(ij) = sqrtf(groupError);
				}
#endif

				//! Aggregation
				remaining_groups -=
					computeAggregationStep1(io_imBasic, weight, mask, group, aggreWeights, 
							mat.agg_window, index, p_params, nSimP);
			}

		//! Weighted aggregation
		computeWeightedAggregation(i_imNoisy, io_imBasic, weight);

		if (p_params.verbose)
		{
			int ntiles = p_crop.ntiles_t * p_crop.ntiles_x * p_crop.ntiles_y;
			int part_idx = p_crop.tile_t * p_crop.ntiles_x * p_crop.ntiles_y +
			               p_crop.tile_y * p_crop.ntiles_x +
			               p_crop.tile_x;

			printf(ANSI_GRN "\x1b[%dF[%d,%d,%d] %05.1f\x1b[%dE" ANSI_RST,
					ntiles - part_idx, p_crop.tile_x, p_crop.tile_y,
					p_crop.tile_t, 100.f, ntiles - part_idx);

			std::cout << std::flush;

			//printf("Processed %d groups\n", group_counter);
		}
	}
	else
	{
		//! Allocate Sizes
		o_imFinal.resize(sz);

		//! Matrices used for Bayes' estimate
		vector<float> groupNoisy(patch_num * patch_dim);
		vector<float> groupBasic(patch_num * patch_dim);
		vector<float> aggreWeights(patch_num, 1.f);
#ifdef DEBUG_COMPUTE_GROUP_ERROR
		vector<float> groupClean(patch_num * patch_dim);
#endif

		int remaining_groups = n_groups;
		for (unsigned pt = 0; pt < sz.frames; pt++)
		for (unsigned py = 0; py < sz.height; py++)
		for (unsigned px = 0; px < sz.width ; px++)
			if (mask(px,py,pt)) //< Only non-seen patches are processed
			{
				group_counter++;

				const unsigned ij  = sz.index(px,py,pt);
				const unsigned ij3 = sz.index(px,py,pt, 0);
				//const unsigned ij3 = (ij / sz.wh) * sz.whc + ij % sz.wh;

				if (p_params.verbose && (group_counter % 100 == 0))
				{
					int ntiles = p_crop.ntiles_t * p_crop.ntiles_x * p_crop.ntiles_y;
					int part_idx = p_crop.tile_t * p_crop.ntiles_x * p_crop.ntiles_y +
					               p_crop.tile_y * p_crop.ntiles_x +
					               p_crop.tile_x;

					printf("\x1b[%dF[%d,%d,%d] %05.1f\x1b[%dE", ntiles - part_idx,
							p_crop.tile_x, p_crop.tile_y, p_crop.tile_t,
							100.f - (float)remaining_groups/(float)(n_groups)*100.f,
							ntiles - part_idx);

					std::cout << std::flush;
				}

				//! Search for similar patches around the reference one
#ifndef DEBUG_COMPUTE_GROUP_ERROR
				unsigned nSimP = estimateSimilarPatchesStep2(i_imNoisy, io_imBasic,
						i_fflow, i_bflow, groupNoisy, groupBasic, index, ij3, p_params);
#else
				unsigned nSimP = estimateSimilarPatchesStep2(i_imNoisy, io_imBasic,
						i_fflow, i_bflow, groupNoisy, groupBasic, index, ij3, p_params,
						i_imClean, groupClean);
#endif

				//! If we use the homogeneous area trick
				bool doBayesEstimate = true;
				if (p_params.useHomogeneousArea)
					doBayesEstimate = !computeHomogeneousAreaStep2(groupNoisy,
							groupBasic, sPx, nSimP, threshold, sz);

				//! Else, use Bayes' estimate
				if (doBayesEstimate)
					variance(ij) = computeBayesEstimateStep2(groupNoisy, groupBasic,
							mat, nInverseFailed, sz, p_params, nSimP, aggreWeights);

#ifdef DEBUG_COMPUTE_GROUP_ERROR
				{
					float groupError = 0, tmp;
					for (int i = 0; i < nSimP * patch_dim; ++i)
						groupError += (tmp = groupNoisy[i] - groupClean[i])*tmp;

					groupError /= (float)(nSimP*patch_dim);
					group_error(ij) = sqrtf(groupError);
				}
#endif

				//! Aggregation
				remaining_groups -=
					computeAggregationStep2(o_imFinal, weight, mask, groupNoisy,
						aggreWeights, mat.agg_window, variance, index, p_params, nSimP);
			}

		//! Weighted aggregation
		computeWeightedAggregation(i_imNoisy, o_imFinal, weight);

		if (p_params.verbose)
		{
			int ntiles = p_crop.ntiles_t * p_crop.ntiles_x * p_crop.ntiles_y;
			int part_idx = p_crop.tile_t * p_crop.ntiles_x * p_crop.ntiles_y +
			               p_crop.tile_y * p_crop.ntiles_x +
			               p_crop.tile_x;

			printf(ANSI_GRN "\x1b[%dF[%d,%d,%d] %05.1f\x1b[%dE" ANSI_RST,
					ntiles - part_idx, p_crop.tile_x, p_crop.tile_y,
					p_crop.tile_t, 100.f, ntiles - part_idx);

			std::cout << std::flush;
		}
	}


	if (nInverseFailed > 0 && p_params.verbose)
		std::cout << "nInverseFailed = " << nInverseFailed << std::endl;

#ifdef DEBUG_SHOW_WEIGHT
	{
		int part_x = p_crop.tile_x;
		int part_y = p_crop.tile_y;
		int part_t = p_crop.tile_t;
		char name[1024];
		sprintf(name, "dump/weight_step%d_%d.%d.%d_%%03d.tif",
				step1 ? 1 : 2, part_x, part_y, part_t);
		weight.saveVideo(name, 1, 1);
	}
	{
		int part_x = p_crop.tile_x;
		int part_y = p_crop.tile_y;
		int part_t = p_crop.tile_t;
		char name[1024];
		sprintf(name, "dump/var_step%d_%d.%d.%d_%%03d.tif",
				step1 ? 1 : 2, part_x, part_y, part_t);
		variance.saveVideo(name, 1, 1);
	}
 #ifdef DEBUG_COMPUTE_GROUP_ERROR
	{
		int part_x = p_crop.tile_x;
		int part_y = p_crop.tile_y;
		int part_t = p_crop.tile_t;
		char name[1024];
		sprintf(name, "dump/pge_step%d_%d.%d.%d_%%03d.tif",
				step1 ? 1 : 2, part_x, part_y, part_t);
		group_error.saveVideo(name, 1, 1);
	}
 #endif
#endif

	return group_counter;
}

//#define VBM3D_SEARCH
#ifdef VBM3D_SEARCH
/**
 * @brief Estimate the best similar patches to a reference one.
 *
 * @param im: contains the noisy image on which distances are processed;
 * @param group: will contain values of similar patches;
 * @param index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
unsigned estimateSimilarPatchesStep1(
	Video<float> const& im
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<std::vector<float> > &group
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
#ifdef DEBUG_COMPUTE_GROUP_ERROR
,	Video<float> const& imClean
,	std::vector<std::vector<float> > & groupClean
#endif
){
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime;
	unsigned nSimP = 1;

	if (params.nSimilarPatches > 1)
	{
		unsigned nSimPFrame = 2; // number of similar patches per frame

		const VideoSize sz = im.sz;
		const bool use_flow = (fflow.sz.width > 0);

		//! Determine search range
		const int sWt_f = params.sizeSearchTimeRangeFwd;
		const int sWt_b = params.sizeSearchTimeRangeBwd;

		//! Coordinates of center of search box
		unsigned px, py, pt, pc;
		im.sz.coords(pidx, px, py, pt, pc);

		//! Temporal search range
		int ranget[2];
#ifdef CENTRED_SEARCH
		ranget[0] = std::max(0, (int)pt -  sWt_b);
		ranget[1] = std::min((int)sz.frames - sPt, (int)pt + sWt_f);
#else
		int shift_t = std::min(0, (int)pt -  sWt_b)
		            + std::max(0, (int)pt +  sWt_f - (int)sz.frames + sPt);

		ranget[0] = std::max(0, (int)pt - sWt_b - shift_t);
		ranget[1] = std::min((int)sz.frames - sPt, (int)pt +  sWt_f - shift_t);
#endif

		//! Redefine size of temporal search range
		int sWt = ranget[1] - ranget[0] + 1;

		//! Number of patches in search region
		int nsrch = 0;

		//! Schedule frames in temporal range
		std::vector<int> srch_ranget;
		srch_ranget.push_back(pt);
		for (int qt = pt+1; qt <= ranget[1]; ++qt) srch_ranget.push_back(qt);
		for (int qt = pt-1; qt >= ranget[0]; --qt) srch_ranget.push_back(qt);

		//! Trajectory of search center
		std::vector<std::vector<int> > traj_cx(nSimPFrame, std::vector<int>(sWt,0) ),
		                               traj_cy(nSimPFrame, std::vector<int>(sWt,0) ),
		                               traj_ct(nSimPFrame, std::vector<int>(sWt,0) );

		std::vector<std::pair<float, unsigned>> distance;

		//! Search
//		printf("p=%d %d %d rt=[%d, %d]\n", px, py, pt,ranget[0], ranget[1]);
		for (int ii = 0; ii < sWt; ++ii)
		{
			int qt = srch_ranget[ii]; // video frame number
			int dt = qt - ranget[0]; // search region frame number
			int dir = std::max(-1, std::min(1, qt - (int)pt)); // direction (forward or backwards from pt)

			//! Determine the centers of the square search regions
			int sWx, sWy;
			std::vector<int> cx0, cy0, ct0; // coordinates of the centers
			if (dir != 0)
			{
				// use a smaller search window
				sWx = round((float)params.sizeSearchWindow / 7. * 5.); //FIXME
				sWy = round((float)params.sizeSearchWindow / 7. * 5.); //FIXME
				for (int i = 0; i < nSimPFrame; i++)
				{
					int cx = traj_cx[i][dt - dir];
					int cy = traj_cy[i][dt - dir];
					int ct = traj_ct[i][dt - dir];
					float cxf = cx + (use_flow ? (dir > 0 ? fflow(cx,cy,ct,0) : bflow(cx,cy,ct,0)) : 0.f);
					float cyf = cy + (use_flow ? (dir > 0 ? fflow(cx,cy,ct,1) : bflow(cx,cy,ct,1)) : 0.f);
					cx0.push_back(std::max(0., std::min((double)sz.width  - 1, round(cxf))));
					cy0.push_back(std::max(0., std::min((double)sz.height - 1, round(cyf))));
				}
			}
			else
			{
				// only a single square search region
				cx0.push_back(px);
				cy0.push_back(py);
				sWx = params.sizeSearchWindow;
				sWy = params.sizeSearchWindow;
			}

//			printf("qt=%d - dt=%d - dir=%d - sWx=%d\n", qt, dt, dir, sWx);

			//! Search the square sub-regions
			std::vector<std::pair<float, unsigned>> distance_t;
			std::unordered_set<int> search_region;
			for (int i = 0; i < cx0.size(); i++)
			{
				int rangex[2];
				int rangey[2];
#ifdef CENTRED_SEARCH
				rangex[0] = std::max(0, cx0[i] - (sWx-1)/2);
				rangey[0] = std::max(0, cy0[i] - (sWy-1)/2);

				rangex[1] = std::min((int)sz.width  - sPx, cx0[i] + (sWx-1)/2);
				rangey[1] = std::min((int)sz.height - sPx, cy0[i] + (sWy-1)/2);
#else
				int shift_x = std::min(0, cx0[i] - (sWx-1)/2);
				int shift_y = std::min(0, cy0[i] - (sWy-1)/2);

				shift_x += std::max(0, cx0[i] + (sWx-1)/2 - (int)sz.width  + sPx);
				shift_y += std::max(0, cy0[i] + (sWy-1)/2 - (int)sz.height + sPx);

				rangex[0] = std::max(0, cx0[i] - (sWx-1)/2 - shift_x);
				rangey[0] = std::max(0, cy0[i] - (sWy-1)/2 - shift_y);

				rangex[1] = std::min((int)sz.width  - sPx, cx0[i] + (sWx-1)/2 - shift_x);
				rangey[1] = std::min((int)sz.height - sPx, cy0[i] + (sWy-1)/2 - shift_y);
#endif

//				printf("c0=[%d, %d] rx=[%d, %d] - ry=[%d, %d]\n",
//						cx0[i], cy0[i], rangex[0], rangex[1], rangey[0], rangey[1]);

				// add elements of the square search region to the union
				for (int qy = rangey[0]; qy <= rangey[1]; qy++)
				for (int qx = rangex[0]; qx <= rangex[1]; qx++)
				if (search_region.insert(sz.index(qx, qy, qt, 0)).second)
				{
					//! Squared L2 distance
					float dist = 0.f, dif;
					for (int ht = 0; ht < sPt; ht++)
					for (int hy = 0; hy < sPx; hy++)
					for (int hx = 0; hx < sPx; hx++)
						dist += (dif = im(px + hx, py + hy, pt + ht)
										 - im(qx + hx, qy + hy, qt + ht)) * dif;

					//! Save distance and corresponding patch index
					distance_t.push_back(std::make_pair(dist, sz.index(qx, qy, qt, 0)));
				}
			}

			//! Keep only the N2 best similar patches
			std::partial_sort(distance_t.begin(), distance_t.begin() + nSimPFrame,
									distance_t.end(), comparaisonFirst);

			for (int i = 0; i < nSimPFrame; i++)
			{
				//! Store best matches in distance vector
				distance.push_back(distance_t[i]);

				//! Store coordinates of best matches for the trajetories
				unsigned cx, cy, ct, cc;
				im.sz.coords(distance_t[i].second, cx, cy, ct, cc);
				traj_cx[i][dt] = cx;
				traj_cy[i][dt] = cy;
				traj_ct[i][dt] = ct;
//
//				printf("\tc=[%d, %d, %d] ", cx, cy, ct);
//				printf("-> dist = %f\n", distance_t[i].first);
			}
		}

		//! Keep only the N2 best similar patches
		nSimP = std::min(params.nSimilarPatches, (unsigned)distance.size());
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
		                  distance.end(), comparaisonFirst);

//		for (int i = 0; i < nSimP; i++)
//		{
//				//! Store coordinates of best matches for the trajetories
//				unsigned cx, cy, ct, cc;
//				im.sz.coords(distance[i].second, cx, cy, ct, cc);
//				printf("final c=[%d, %d, %d] ", cx, cy, ct);
//				printf("-> dist = %f\n", distance[i].first);
//		}
//
//		while (1);

		//! Store indices of most similar patches
		for (unsigned n = 0; n < nSimP; n++)
			index[n] = distance[n].second;
	}
	else // nSimilarPatches == 1
		index[0] = pidx;


	//! Stack selected patches into the 3D group
	const unsigned w   = im.sz.width;
	const unsigned wh  = im.sz.wh;
	const unsigned whc = im.sz.whc;
#ifdef USE_FFTW
	for (unsigned c  = 0; c < im.sz.channels; c++)
	for (unsigned n  = 0, k = 0; n < nSimP; n++)
	for (unsigned ht = 0;        ht < sPt; ht++)
	for (unsigned hy = 0;        hy < sPx; hy++)
	for (unsigned hx = 0;        hx < sPx; hx++, k++)
		group[c][k] = im(c * wh + index[n] + ht * whc + hy * w + hx);
#else
	for (unsigned c  = 0; c < im.sz.channels; c++)
	for (unsigned ht = 0, k = 0; ht < sPt; ht++)
	for (unsigned hy = 0;        hy < sPx; hy++)
	for (unsigned hx = 0;        hx < sPx; hx++)
	for (unsigned n  = 0; n < nSimP; n++, k++)
		group[c][k] = im(c * wh + index[n] + ht * whc + hy * w + hx);
#endif

	/* 000  pixels from all patches
	 * 001  pixels from all patches
	 * ...
	 * spt,spx,spx pixels from all patches
	 */

#ifdef DEBUG_COMPUTE_GROUP_ERROR
#ifdef USE_FFTW
	for (unsigned c  = 0; c < im.sz.channels; c++)
	for (unsigned n  = 0, k = 0; n < nSimP; n++)
	for (unsigned ht = 0;        ht < sPt; ht++)
	for (unsigned hy = 0;        hy < sPx; hy++)
	for (unsigned hx = 0;        hx < sPx; hx++, k++)
		groupClean[c][k] = imClean(c * wh + index[n] + ht * whc + hy * w + hx);
#else
	for (unsigned c  = 0; c < im.sz.channels; c++)
	for (unsigned ht = 0, k = 0; ht < sPt; ht++)
	for (unsigned hy = 0;        hy < sPx; hy++)
	for (unsigned hx = 0;        hx < sPx; hx++)
	for (unsigned n  = 0; n < nSimP; n++, k++)
		groupClean[c][k] = imClean(c * wh + index[n] + ht * whc + hy * w + hx);
#endif
#endif

	return nSimP;
}
#else
/**
 * @brief Estimate the best similar patches to a reference one.
 *
 * @param im: contains the noisy image on which distances are processed;
 * @param group: will contain values of similar patches;
 * @param index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
unsigned estimateSimilarPatchesStep1(
	Video<float> const& im
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<std::vector<float> > &group
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
#ifdef DEBUG_COMPUTE_GROUP_ERROR
,	Video<float> const& imClean
,	std::vector<std::vector<float> > & groupClean
#endif
){
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime;
	unsigned nSimP = 1;

	if (params.nSimilarPatches > 1)
	{
		const VideoSize sz = im.sz;
		const bool use_flow = (fflow.sz.width > 0);

		//! Determine search range
		int sWx   = params.sizeSearchWindow;
		int sWy   = params.sizeSearchWindow;
		const int sWt_f = params.sizeSearchTimeRangeFwd;
		const int sWt_b = params.sizeSearchTimeRangeBwd;

		//! Coordinates of center of search box
		unsigned px, py, pt, pc;
		im.sz.coords(pidx, px, py, pt, pc);

		//! Temporal search range
		int ranget[2];
#ifdef CENTRED_SEARCH
		ranget[0] = std::max(0, (int)pt -  sWt_b);
		ranget[1] = std::min((int)sz.frames - sPt, (int)pt + sWt_f);
#else
		int shift_t = std::min(0, (int)pt -  sWt_b)
		            + std::max(0, (int)pt +  sWt_f - (int)sz.frames + sPt);

		ranget[0] = std::max(0, (int)pt - sWt_b - shift_t);
		ranget[1] = std::min((int)sz.frames - sPt, (int)pt +  sWt_f - shift_t);
#endif

		//! Redefine size of temporal search range
		int sWt = ranget[1] - ranget[0] + 1;

		//! Allocate vector of patch distances
		std::vector<std::pair<float, unsigned> > distance(sWx * sWy * sWt);

		//! Number of patches in search region
		int nsrch = 0;

		//! Schedule frames in temporal range
		std::vector<int> srch_ranget;
		srch_ranget.push_back(pt);
		for (int qt = pt+1; qt <= ranget[1]; ++qt) srch_ranget.push_back(qt);
		for (int qt = pt-1; qt >= ranget[0]; --qt) srch_ranget.push_back(qt);

		//! Trajectory of search center
		std::vector<int> cx(sWt,0), cy(sWt,0), ct(sWt,0);

		//! Search
		for (int ii = 0; ii < sWt; ++ii)
		{
			int qt = srch_ranget[ii]; // video frame number
			int dt = qt - ranget[0]; // search region frame number
			int dir = std::max(-1, std::min(1, qt - (int)pt)); // direction (forward or backwards from pt)

			//! Integrate optical flow to new center
			if (dir != 0)
			{
				int cx0 = cx[dt - dir];
				int cy0 = cy[dt - dir];
				int ct0 = ct[dt - dir];
				float cx_f = cx0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,0) : bflow(cx0,cy0,ct0,0)) : 0.f);
				float cy_f = cy0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,1) : bflow(cx0,cy0,ct0,1)) : 0.f);
				cx[dt] = std::max(0., std::min((double)sz.width  - 1, round(cx_f)));
				cy[dt] = std::max(0., std::min((double)sz.height - 1, round(cy_f)));
				ct[dt] = qt;
			}
			else
			{
				cx[dt] = px;
				cy[dt] = py;
				ct[dt] = pt;
			}

			//! Spatial search range
			int rangex[2];
			int rangey[2];
#ifdef CENTRED_SEARCH
			rangex[0] = std::max(0, cx[dt] - (sWx-1)/2);
			rangey[0] = std::max(0, cy[dt] - (sWy-1)/2);

			rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2);
			rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2);
#else
			int shift_x = std::min(0, cx[dt] - (sWx-1)/2);
			int shift_y = std::min(0, cy[dt] - (sWy-1)/2);

			shift_x += std::max(0, cx[dt] + (sWx-1)/2 - (int)sz.width  + sPx);
			shift_y += std::max(0, cy[dt] + (sWy-1)/2 - (int)sz.height + sPx);

			rangex[0] = std::max(0, cx[dt] - (sWx-1)/2 - shift_x);
			rangey[0] = std::max(0, cy[dt] - (sWy-1)/2 - shift_y);

			rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2 - shift_x);
			rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2 - shift_y);
#endif

			//! Compute distance between patches in search range
			for (int qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
			for (int qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
			{
				//! Squared L2 distance
				float dist = 0.f, dif;
				for (int ht = 0; ht < sPt; ht++)
				for (int hy = 0; hy < sPx; hy++)
				for (int hx = 0; hx < sPx; hx++)
					dist += (dif = im(px + hx, py + hy, pt + ht)
					             - im(qx + hx, qy + hy, qt + ht)) * dif;

				//! Save distance and corresponding patch index
				distance[nsrch++] = std::make_pair(dist, sz.index(qx, qy, qt, 0));
			}
		}

		distance.resize(nsrch);

		//! Keep only the N2 best similar patches
		nSimP = std::min(params.nSimilarPatches, (unsigned)distance.size());
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
		                  distance.end(), comparaisonFirst);

		//! Add more patches if their distance is below a threshold
		const float threshold = std::max(params.tau, distance[nSimP - 1].first);
		nSimP = 0;
		for (unsigned n = 0; n < distance.size(); n++)
			if (distance[n].first <= threshold)
				index[nSimP++] = distance[n].second;

//		if (nSimP < params.nSimilarPatches)
//		{
//			printf("SR1 [%d,%d,%d] ~ [%d-%d, %d-%d, %d-%d] - nsim = %d\n", 
//					px,py,pt,rangex[0], rangex[1], rangey[0], rangey[1], ranget[0], ranget[1], nSimP);
//		}

//		if (nSimP > params.nSimilarPatches)
//			printf("SR1 [%d,%d,%d] ~ nsim = %d\n", px,py,pt, nSimP);

//		for (int i = 0; i < nSimP; i++)
//		{
//			unsigned cx, cy, ct, cc;
//			im.sz.coords(distance[i].second, cx, cy, ct, cc);
//			printf("d[%03d] = %g - p = [%02d,%02d,%02d,%2d]\n", i, distance[i].first,
//					cx, cy, ct, cc);
//		}
	}
	else // nSimilarPatches == 1
		index[0] = pidx;

	//! Stack selected patches into the 3D group
	const unsigned w   = im.sz.width;
	const unsigned wh  = im.sz.wh;
	const unsigned whc = im.sz.whc;
#ifdef USE_FFTW
	for (unsigned c  = 0; c < im.sz.channels; c++)
	for (unsigned n  = 0, k = 0; n < nSimP; n++)
	for (unsigned ht = 0;        ht < sPt; ht++)
	for (unsigned hy = 0;        hy < sPx; hy++)
	for (unsigned hx = 0;        hx < sPx; hx++, k++)
		group[c][k] = im(c * wh + index[n] + ht * whc + hy * w + hx);
#else
	for (unsigned c  = 0; c < im.sz.channels; c++)
	for (unsigned ht = 0, k = 0; ht < sPt; ht++)
	for (unsigned hy = 0;        hy < sPx; hy++)
	for (unsigned hx = 0;        hx < sPx; hx++)
	for (unsigned n  = 0; n < nSimP; n++, k++)
		group[c][k] = im(c * wh + index[n] + ht * whc + hy * w + hx);
#endif

	/* 000  pixels from all patches
	 * 001  pixels from all patches
	 * ...
	 * spt,spx,spx pixels from all patches
	 */

#ifdef DEBUG_COMPUTE_GROUP_ERROR
#ifdef USE_FFTW
	for (unsigned c  = 0; c < im.sz.channels; c++)
	for (unsigned n  = 0, k = 0; n < nSimP; n++)
	for (unsigned ht = 0;        ht < sPt; ht++)
	for (unsigned hy = 0;        hy < sPx; hy++)
	for (unsigned hx = 0;        hx < sPx; hx++, k++)
		groupClean[c][k] = imClean(c * wh + index[n] + ht * whc + hy * w + hx);
#else
	for (unsigned c  = 0; c < im.sz.channels; c++)
	for (unsigned ht = 0, k = 0; ht < sPt; ht++)
	for (unsigned hy = 0;        hy < sPx; hy++)
	for (unsigned hx = 0;        hx < sPx; hx++)
	for (unsigned n  = 0; n < nSimP; n++, k++)
		groupClean[c][k] = imClean(c * wh + index[n] + ht * whc + hy * w + hx);
#endif
#endif

	return nSimP;
}
#endif

#ifdef VBM3D_SEARCH
/**
 * @brief Keep from all near patches the similar ones to the reference patch
 * for the second step.
 *
 * @param imNoisy: contains the original noisy image;
 * @param imBasic: contains the basic estimation;
 * @param fflow: forward optical flow (optional)
 * @param bflow: backward optical flow (optional)
 * @param groupNoisy: will contain similar patches for all channels of imNoisy;
 * @param groupBasic: will contain similar patches for all channels of i_imBasic;
 * @param index: will contain index of similar patches;
 * @param pidx: index of the reference patch;
 * @param params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep2(
	Video<float> const& imNoisy
,	Video<float> const& imBasic
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<float> &groupNoisy
,	std::vector<float> &groupBasic
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
#ifdef DEBUG_COMPUTE_GROUP_ERROR
,	Video<float> const& imClean
,	std::vector<float> & groupClean
#endif
){
	const VideoSize sz = imNoisy.sz;
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime;
	const int chnls = sz.channels;
	unsigned nSimP = 1;

	if (params.nSimilarPatches > 1)
	{
		unsigned nSimPFrame = 2; // number of similar patches per frame

		const bool use_flow = (fflow.sz.width > 0);

		//! Determine search range
		const int sWt_f = params.sizeSearchTimeRangeFwd;
		const int sWt_b = params.sizeSearchTimeRangeBwd;

		//! Coordinates of center of search box
		unsigned px, py, pt, pc;
		sz.coords(pidx, px, py, pt, pc);

		//! Temporal search range
		int ranget[2];
#ifdef CENTRED_SEARCH
		ranget[0] = std::max(0, (int)pt -  sWt_b);
		ranget[1] = std::min((int)sz.frames - sPt, (int)pt + sWt_f);
#else
		int shift_t = std::min(0, (int)pt -  sWt_b)
		            + std::max(0, (int)pt +  sWt_f - (int)sz.frames + sPt);

		ranget[0] = std::max(0, (int)pt - sWt_b - shift_t);
		ranget[1] = std::min((int)sz.frames - sPt, (int)pt +  sWt_f - shift_t);
#endif

		//! Redefine size of temporal search range
		int sWt = ranget[1] - ranget[0] + 1;

		//! Number of patches in search region
		int nsrch = 0;

		//! Schedule frames in temporal range
		std::vector<int> srch_ranget;
		srch_ranget.push_back(pt);
		for (int qt = pt+1; qt <= ranget[1]; ++qt) srch_ranget.push_back(qt);
		for (int qt = pt-1; qt >= ranget[0]; --qt) srch_ranget.push_back(qt);

		//! Trajectories of search centers
		std::vector<std::vector<int> > traj_cx(nSimPFrame, std::vector<int>(sWt,0) ),
		                               traj_cy(nSimPFrame, std::vector<int>(sWt,0) ),
		                               traj_ct(nSimPFrame, std::vector<int>(sWt,0) );

		std::vector<std::pair<float, unsigned>> distance;

		//! Search
//		printf("p=%d %d %d rt=[%d, %d]\n", px, py, pt,ranget[0], ranget[1]);
		for (int ii = 0; ii < sWt; ++ii)
		{
			int qt = srch_ranget[ii]; // video frame number
			int dt = qt - ranget[0]; // search region frame number
			int dir = std::max(-1, std::min(1, qt - (int)pt)); // direction (forward or backwards from pt)

			//! Determine the centers of the square search regions
			int sWx, sWy;
			std::vector<int> cx0, cy0, ct0; // coordinates of the centers
			if (dir != 0)
			{
				// use a smaller search window
				sWx = round((float)params.sizeSearchWindow / 7. * 5.); //FIXME
				sWy = round((float)params.sizeSearchWindow / 7. * 5.); //FIXME
				for (int i = 0; i < nSimPFrame; i++)
				{
					int cx = traj_cx[i][dt - dir];
					int cy = traj_cy[i][dt - dir];
					int ct = traj_ct[i][dt - dir];
					float cxf = cx + (use_flow ? (dir > 0 ? fflow(cx,cy,ct,0) : bflow(cx,cy,ct,0)) : 0.f);
					float cyf = cy + (use_flow ? (dir > 0 ? fflow(cx,cy,ct,1) : bflow(cx,cy,ct,1)) : 0.f);
					cx0.push_back(std::max(0., std::min((double)sz.width  - 1, round(cxf))));
					cy0.push_back(std::max(0., std::min((double)sz.height - 1, round(cyf))));
				}
			}
			else
			{
				// only a single square search region
				cx0.push_back(px);
				cy0.push_back(py);
				sWx = params.sizeSearchWindow;
				sWy = params.sizeSearchWindow;
			}

//			printf("qt=%d - dt=%d - dir=%d - sWx=%d\n", qt, dt, dir, sWx);

			//! Search the square sub-regions
			std::vector<std::pair<float, unsigned>> distance_t;
			std::unordered_set<int> search_region;
			for (int i = 0; i < cx0.size(); i++)
			{
				int rangex[2];
				int rangey[2];
#ifdef CENTRED_SEARCH
				rangex[0] = std::max(0, cx0[i] - (sWx-1)/2);
				rangey[0] = std::max(0, cy0[i] - (sWy-1)/2);

				rangex[1] = std::min((int)sz.width  - sPx, cx0[i] + (sWx-1)/2);
				rangey[1] = std::min((int)sz.height - sPx, cy0[i] + (sWy-1)/2);
#else
				int shift_x = std::min(0, cx0[i] - (sWx-1)/2);
				int shift_y = std::min(0, cy0[i] - (sWy-1)/2);

				shift_x += std::max(0, cx0[i] + (sWx-1)/2 - (int)sz.width  + sPx);
				shift_y += std::max(0, cy0[i] + (sWy-1)/2 - (int)sz.height + sPx);

				rangex[0] = std::max(0, cx0[i] - (sWx-1)/2 - shift_x);
				rangey[0] = std::max(0, cy0[i] - (sWy-1)/2 - shift_y);

				rangex[1] = std::min((int)sz.width  - sPx, cx0[i] + (sWx-1)/2 - shift_x);
				rangey[1] = std::min((int)sz.height - sPx, cy0[i] + (sWy-1)/2 - shift_y);
#endif

//				printf("c0=[%d, %d] rx=[%d, %d] - ry=[%d, %d]\n",
//						cx0[i], cy0[i], rangex[0], rangex[1], rangey[0], rangey[1]);

				// add elements of the square search region to the union
				for (int qy = rangey[0]; qy <= rangey[1]; qy++)
				for (int qx = rangex[0]; qx <= rangex[1]; qx++)
				if (search_region.insert(sz.index(qx, qy, qt, 0)).second)
				{
					//! Squared L2 distance
					float dist = 0.f, dif;
					for (int c = 0; c < chnls; c++)
					for (int ht = 0; ht < sPt; ht++)
					for (int hy = 0; hy < sPx; hy++)
					for (int hx = 0; hx < sPx; hx++)
						dist += (dif = imBasic(px + hx, py + hy, pt + ht)
										 - imBasic(qx + hx, qy + hy, qt + ht)) * dif;

					//! Save distance and corresponding patch index
					distance_t.push_back(std::make_pair(dist, sz.index(qx, qy, qt, 0)));
				}
			}

			//! Keep only the N2 best similar patches
			std::partial_sort(distance_t.begin(), distance_t.begin() + nSimPFrame,
									distance_t.end(), comparaisonFirst);

			for (int i = 0; i < nSimPFrame; i++)
			{
				//! Store best matches in distance vector
				distance.push_back(distance_t[i]);

				//! Store coordinates of best matches for the trajetories
				unsigned cx, cy, ct, cc;
				sz.coords(distance_t[i].second, cx, cy, ct, cc);
				traj_cx[i][dt] = cx;
				traj_cy[i][dt] = cy;
				traj_ct[i][dt] = ct;
//
//				printf("\tc=[%d, %d, %d] ", cx, cy, ct);
//				printf("-> dist = %f\n", distance_t[i].first);
			}
		}

		//! Keep only the N2 best similar patches
		nSimP = std::min(params.nSimilarPatches, (unsigned)distance.size());
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
		                  distance.end(), comparaisonFirst);

//		for (int i = 0; i < nSimP; i++)
//		{
//				//! Store coordinates of best matches for the trajetories
//				unsigned cx, cy, ct, cc;
//				im.sz.coords(distance[i].second, cx, cy, ct, cc);
//				printf("final c=[%d, %d, %d] ", cx, cy, ct);
//				printf("-> dist = %f\n", distance[i].first);
//		}
//
//		while (1);

		//! Store indices of most similar patches
		for (unsigned n = 0; n < nSimP; n++)
			index[n] = distance[n].second;
	}
	else // nSimilarPatches == 1
		index[0] = pidx;

	//! Save similar patches into 3D groups
	const unsigned w   = sz.width;
	const unsigned wh  = sz.wh;
	const unsigned whc = sz.whc;
#ifdef USE_FFTW
	for (unsigned c  = 0, k = 0; c < chnls; c++)
	for (unsigned n  = 0; n < nSimP; n++)
	for (unsigned ht = 0; ht < sPt; ht++)
	for (unsigned hy = 0; hy < sPx; hy++)
	for (unsigned hx = 0; hx < sPx; hx++, k++)
	{
		groupNoisy[k] = imNoisy(c * wh + index[n] + ht * whc + hy * w + hx);
		groupBasic[k] = imBasic(c * wh + index[n] + ht * whc + hy * w + hx);
 #ifdef DEBUG_COMPUTE_GROUP_ERROR
		groupClean[k] = imClean(c * wh + index[n] + ht * whc + hy * w + hx);
 #endif
	}
#else
	for (unsigned c = 0, k = 0; c < chnls; c++)
	for (unsigned ht = 0; ht < sPt; ht++)
	for (unsigned hy = 0; hy < sPx; hy++)
	for (unsigned hx = 0; hx < sPx; hx++)
	for (unsigned n = 0 ; n < nSimP; n++, k++)
	{
		groupNoisy[k] = imNoisy(c * wh + index[n] + ht * whc + hy * w + hx);
		groupBasic[k] = imBasic(c * wh + index[n] + ht * whc + hy * w + hx);
 #ifdef DEBUG_COMPUTE_GROUP_ERROR
		groupClean[k] = imClean(c * wh + index[n] + ht * whc + hy * w + hx);
 #endif
	}
#endif

	return nSimP;
}
#else
/**
 * @brief Keep from all near patches the similar ones to the reference patch
 * for the second step.
 *
 * @param imNoisy: contains the original noisy image;
 * @param imBasic: contains the basic estimation;
 * @param fflow: forward optical flow (optional)
 * @param bflow: backward optical flow (optional)
 * @param groupNoisy: will contain similar patches for all channels of imNoisy;
 * @param groupBasic: will contain similar patches for all channels of i_imBasic;
 * @param index: will contain index of similar patches;
 * @param pidx: index of the reference patch;
 * @param params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep2(
	Video<float> const& imNoisy
,	Video<float> const& imBasic
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<float> &groupNoisy
,	std::vector<float> &groupBasic
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
#ifdef DEBUG_COMPUTE_GROUP_ERROR
,	Video<float> const& imClean
,	std::vector<float> & groupClean
#endif
){
	const VideoSize sz = imNoisy.sz;
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime;
	const int chnls = sz.channels;
	unsigned nSimP = 1;

	if (params.nSimilarPatches > 1)
	{
		bool use_flow = (fflow.sz.width > 0);

		//! Coordinates of center of search box
		unsigned px, py, pt, pc;
		sz.coords(pidx, px, py, pt, pc);

		//! Determine search range
		int sWx   = params.sizeSearchWindow;
		int sWy   = params.sizeSearchWindow;
		const int sWt_f = params.sizeSearchTimeRangeFwd;
		const int sWt_b = params.sizeSearchTimeRangeBwd;

		//! Temporal search range
		int ranget[2];
#ifdef CENTRED_SEARCH
		ranget[0] = std::max(0, (int)pt -  sWt_b);
		ranget[1] = std::min((int)sz.frames - sPt, (int)pt + sWt_f);
#else
		int shift_t = std::min(0, (int)pt -  sWt_b)  
		            + std::max(0, (int)pt +  sWt_f - (int)sz.frames + sPt); 

		ranget[0] = std::max(0, (int)pt - sWt_b - shift_t);
		ranget[1] = std::min((int)sz.frames - sPt, (int)pt +  sWt_f - shift_t);
#endif

		//! Redefine size of temporal search range
		int sWt = ranget[1] - ranget[0] + 1;

		//! Allocate vector of patch distances
		std::vector<std::pair<float, unsigned> > distance(sWx * sWy * sWt);

		//! Number of patches in search region
		int nsrch = 0;

		//! Schedule frames in temporal range
		std::vector<int> srch_ranget;
		srch_ranget.push_back(pt);
		for (int qt = pt+1; qt <= ranget[1]; ++qt) srch_ranget.push_back(qt);
		for (int qt = pt-1; qt >= ranget[0]; --qt) srch_ranget.push_back(qt);

		//! Trajectory of search center
		std::vector<int> cx(sWt,0), cy(sWt,0), ct(sWt,0);

		//! Search
		for (int ii = 0; ii < sWt; ++ii)
		{
			int qt = srch_ranget[ii]; // video frame number
			int dt = qt - ranget[0]; // search region frame number
			int dir = std::max(-1, std::min(1, qt - (int)pt)); // direction (forward or backwards from pt)

			//! Integrate optical flow to new center
			if (dir != 0)
			{
				int cx0 = cx[dt - dir];
				int cy0 = cy[dt - dir];
				int ct0 = ct[dt - dir];
//				if (cx0 > (int)bflow.sz.width - 1 || cx0 < 0 ||
//				    cy0 > (int)bflow.sz.height- 1 || cy0 < 0 ||
//				    ct0 > (int)bflow.sz.frames- 1 || ct0 < 0 )
//					printf("sz = %d,%d,%d ~ c = %d,%d,%d\n",bflow.sz.width, bflow.sz.height, bflow.sz.frames, cx0, cy0, ct0);
				float cx_f = cx0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,0) : bflow(cx0,cy0,ct0,0)) : 0.f);
				float cy_f = cy0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,1) : bflow(cx0,cy0,ct0,1)) : 0.f);
				cx[dt] = std::max(0., std::min((double)sz.width  - 1, round(cx_f)));
				cy[dt] = std::max(0., std::min((double)sz.height - 1, round(cy_f)));
				ct[dt] = qt;
			}
			else
			{
				cx[dt] = px;
				cy[dt] = py;
				ct[dt] = pt;
			}

			//! Spatial search range
			int rangex[2];
			int rangey[2];
#ifdef CENTRED_SEARCH
			rangex[0] = std::max(0, cx[dt] - (sWx-1)/2);
			rangey[0] = std::max(0, cy[dt] - (sWy-1)/2);

			rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2);
			rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2);
#else
			int shift_x = std::min(0, cx[dt] - (sWx-1)/2); 
			int shift_y = std::min(0, cy[dt] - (sWy-1)/2); 

			shift_x += std::max(0, cx[dt] + (sWx-1)/2 - (int)sz.width  + sPx); 
			shift_y += std::max(0, cy[dt] + (sWy-1)/2 - (int)sz.height + sPx); 

			rangex[0] = std::max(0, cx[dt] - (sWx-1)/2 - shift_x);
			rangey[0] = std::max(0, cy[dt] - (sWy-1)/2 - shift_y);

			rangex[1] = std::min((int)sz.width  - sPx, cx[dt] + (sWx-1)/2 - shift_x);
			rangey[1] = std::min((int)sz.height - sPx, cy[dt] + (sWy-1)/2 - shift_y);
#endif

			//! Compute distance between patches in search range
			for (int qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
			for (int qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
			{
				//! Squared L2 distance
				float dist = 0.f, dif;
				for (int c = 0; c < chnls; c++)
				for (int ht = 0; ht < sPt; ht++)
				for (int hy = 0; hy < sPx; hy++)
				for (int hx = 0; hx < sPx; hx++)
					dist += (dif = imBasic(px + hx, py + hy, pt + ht, c)
					             - imBasic(qx + hx, qy + hy, qt + ht, c) ) * dif;

				//! Save distance and corresponding patch index
				distance[nsrch++] = std::make_pair(dist, sz.index(qx, qy, qt, 0));
			}
		}

		distance.resize(nsrch);

		//! Keep only the nSimilarPatches best similar patches
		nSimP = std::min(params.nSimilarPatches, (unsigned)distance.size());
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
		                  distance.end(), comparaisonFirst);

		//! Add more patches if their distance is bellow the threshold
		const float threshold = (params.tau > distance[nSimP - 1].first ?
		                         params.tau : distance[nSimP - 1].first);
		nSimP = 0;
		for (unsigned n = 0; n < distance.size(); n++)
			if (distance[n].first <= threshold)
				index[nSimP++] = distance[n].second;

//		if (nSimP > params.nSimilarPatches)
//			printf("SR2 [%d,%d,%d] ~ nsim = %d ~ nsim ratio = %f\n", px,py,pt, nSimP, (float)nSimP/(float)(sWx*sWy*sWt));
	}
	else // nSimilarPatches == 1
		index[0] = pidx;

	//! Save similar patches into 3D groups
	const unsigned w   = sz.width;
	const unsigned wh  = sz.wh;
	const unsigned whc = sz.whc;
#ifdef USE_FFTW
	for (unsigned c  = 0, k = 0; c < chnls; c++)
	for (unsigned n  = 0; n < nSimP; n++)
	for (unsigned ht = 0; ht < sPt; ht++)
	for (unsigned hy = 0; hy < sPx; hy++)
	for (unsigned hx = 0; hx < sPx; hx++, k++)
	{
		groupNoisy[k] = imNoisy(c * wh + index[n] + ht * whc + hy * w + hx);
		groupBasic[k] = imBasic(c * wh + index[n] + ht * whc + hy * w + hx);
 #ifdef DEBUG_COMPUTE_GROUP_ERROR
		groupClean[k] = imClean(c * wh + index[n] + ht * whc + hy * w + hx);
 #endif
	}
#else
	for (unsigned c = 0, k = 0; c < chnls; c++)
	for (unsigned ht = 0; ht < sPt; ht++)
	for (unsigned hy = 0; hy < sPx; hy++)
	for (unsigned hx = 0; hx < sPx; hx++)
	for (unsigned n = 0 ; n < nSimP; n++, k++)
	{
		groupNoisy[k] = imNoisy(c * wh + index[n] + ht * whc + hy * w + hx);
		groupBasic[k] = imBasic(c * wh + index[n] + ht * whc + hy * w + hx);
 #ifdef DEBUG_COMPUTE_GROUP_ERROR
		groupClean[k] = imClean(c * wh + index[n] + ht * whc + hy * w + hx);
 #endif
	}
#endif

	return nSimP;
}
#endif

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group: contains for each channels values of similar patches. If an homogeneous area
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
	std::vector<std::vector<float> > &io_group
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
		stdDev += computeStdDeviation(io_group[c], p_sP * p_sP, p_nSimP, 1);

	//! If we are in an homogeneous area
	if (stdDev < p_threshold)
	{
		for (unsigned c = 0; c < p_imSize.channels; c++)
		{
			float mean = 0.f;

			for (unsigned k = 0; k < N; k++) mean += io_group[c][k];
			mean /= (float) N;

			for (unsigned k = 0; k < N; k++) io_group[c][k] = mean;
		}
		return 1;
	}
	else return 0;
}

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_groupNoisy: inputs values of similar patches for the noisy video;
 *                         if the area is classified as homogeneous, outputs the
 *                         average of all pixels in all patches.
 * @param i_groupBasic: contains values of similar patches for the basic video.
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_imSize: size of the video.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep2(
	std::vector<float> &io_groupNoisy
,	std::vector<float> const &i_groupBasic
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
){
	//! Parameters
	const unsigned sP2 = p_sP * p_sP;
	const unsigned sPC = sP2 * p_imSize.channels;

	//! Compute the standard deviation of the set of patches
	const float stdDev = computeStdDeviation(io_groupNoisy, sP2, p_nSimP, p_imSize.channels);

	//! If we are in an homogeneous area
	if (stdDev < p_threshold)
	{
		for (unsigned c = 0; c < p_imSize.channels; c++)
		{
				float mean = 0.f;
				for (unsigned n = 0; n < p_nSimP; n++)
				for (unsigned k = 0; k < sP2; k++)
					mean += i_groupBasic[n * sPC + c * sP2 + k];

				mean /= float(sP2 * p_nSimP);

				for (unsigned n = 0; n < p_nSimP; n++)
				for (unsigned k = 0; k < sP2; k++)
					io_groupNoisy[n * sPC + c * sP2 + k] = mean;
		}
		return 1;
	}
	else return 0;
}

#ifdef USE_FFTW
/**
 * @brief Implementation of computeBayesEstimateStep1 using an 
 * external basis provided by the user in the workspace i_mat.covEigVecs.
 *
 * See computeBayesEstimateStep1 for information about the arguments.
 **/
float computeBayesEstimateStep1_externalBasisFFTW(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	std::vector<std::vector<float> > &aggreWeights 
){
	//! Parameters initialization
	const float beta_sigma = p_params.beta * p_params.sigma;
	const float beta_sigma2 = beta_sigma * beta_sigma;
#ifndef USE_BETA_FOR_VARIANCE
	const float sigma2 = p_params.sigma * p_params.sigma;
#else
	const float sigma2 = beta_sigma2;
#endif
	const unsigned sPC = p_params.sizePatch * p_params.sizePatch
	                   * p_params.sizePatchTime;
	const unsigned r   = p_params.rank ? sPC : 0; // XXX FIXME TODO

	const float inSimP = 1.f/(float)p_nSimP;

	//! Variances
	float  rank_variance = 0.f;
	float total_variance = 0.f;

//	float entropy = 0;
//	float total_weights = 0;

	for (unsigned c = 0; c < io_group.size(); c++)
	{
#ifndef DCT_DONT_CENTER1
		//! Center 3D group
		i_mat.baricenter.assign(sPC, 0.f);
		for (int i = 0; i < p_nSimP; ++i)
		for (int k = 0; k < sPC; ++k)
			i_mat.baricenter[k] += io_group[c][i*sPC + k] * inSimP;

		for (int i = 0; i < p_nSimP; ++i)
		for (int k = 0; k < sPC; ++k)
			io_group[c][i*sPC + k] -= i_mat.baricenter[k];
#endif

		if (r > 0)
		{
			//! Forward DCT
			globalDCT.forward(io_group[c]);

			//! Compute variance over each component
			i_mat.covEigVals.assign(sPC, 0.f);
			for (int i = 0; i < p_nSimP; ++i)
			for (int k = 0; k < sPC; ++k)
				i_mat.covEigVals[k] += io_group[c][i*sPC + k]
				                     * io_group[c][i*sPC + k];

			for (int k = 0; k < sPC; ++k)
			{
				i_mat.covEigVals[k] *= inSimP;
				total_variance += i_mat.covEigVals[k]/(float)sPC/(float)io_group.size();
			}

//			// XXX DEBUG - compute entropy of variances
//			{ 
//				float normalization = 0;
//				for (int k = 0; k < sPC; ++k)
//					normalization += i_mat.covEigVals[k];
//
//				for (int k = 0; k < sPC; ++k)
//					entropy -= (i_mat.covEigVals[k]/normalization) * log(i_mat.covEigVals[k]/normalization);
//			}

			//! Compute aggregation weights
			if (p_params.aggreGammaPatch)
			{
				aggreWeights[c].assign(p_nSimP, 0.f);

				//! Mahalanobis distance
				float tmp;
				for (int n = 0; n < p_nSimP; ++n)
				for (int i = 0; i < sPC    ; ++i)
						aggreWeights[c][n] += (tmp = io_group[c][n*sPC + i])
						                    *  tmp / std::max(i_mat.covEigVals[i],sigma2);

				//! Exponentiate
				float igamma = 1.f/2.f/p_params.aggreGammaPatch;
				for(int n = 0; n < p_nSimP; ++n)
					aggreWeights[c][n] = std::exp(-igamma * aggreWeights[c][n]);
			}
			else
				aggreWeights[c].assign(p_nSimP, 1.f);

			//! Compute eigenvalues-based coefficients of Bayes' filter
			for (unsigned k = 0; k < r; ++k)
			{
				if (i_mat.covEigVals[k] == 0) continue;

				rank_variance += i_mat.covEigVals[k];
				float var = i_mat.covEigVals[k] - sigma2;

#if defined(THRESHOLD_WEIGHTS1)
				var = std::max(0.f, var);
#elif defined(LI_ZHANG_DAI1)
				var += sigma2; // add back sigma2
				var = (var < 4.f*sigma2) ? 0.f
				    : (var - 2.f*sigma2 + sqrtf(var*(var - 4.f*sigma2)))*.5f;
#endif

#ifndef LINEAR_THRESHOLDING1
				i_mat.covEigVals[k] = (fabs(var + beta_sigma2) > 1e-8f) ? 
					                   var / ( var + beta_sigma2 ) : 
											 0.f;
#else
				// this doesn't make sense for Li-Zhang-Dai's variance
				i_mat.covEigVals[k] = (var > 0.f) ? 1.f : 0.f;
#endif
			}

//			// XXX DEBUG compute total weights
//			for (int k = 0; k < sPC; ++k)
//				total_weights += i_mat.covEigVals[k]/(float)sPC/(float)io_group.size();

			//! W*Z
			for (int i = 0; i < p_nSimP; ++i)
			for (int k = 0; k < sPC; ++k)
				io_group[c][i*sPC + k] *= i_mat.covEigVals[k];

			//! Inverse DCT
			globalDCT.inverse(io_group[c]);

#ifndef DCT_DONT_CENTER1
			//! Add baricenter
			for (int i = 0; i < p_nSimP; ++i)
			for (int k = 0; k < sPC; ++k)
				io_group[c][i*sPC + k] += i_mat.baricenter[k];
#endif
		}
		else
		{
#ifndef DCT_DONT_CENTER1
			//! rank = 0: set as baricenter
			for (int i = 0; i < p_nSimP; ++i)
			for (int k = 0; k < sPC; ++k)
				io_group[c][i*sPC + k] = i_mat.baricenter[k];
#endif

			//! Avoid 0/0 in return statement
			total_variance = 1.f;
		}
	}

	// compute group aggregation weights
	if (p_params.aggreGammaGroup)
	{
		const float aggreSigma = p_params.aggreGammaGroup*p_params.sigma;

#ifdef GAUSSIAN_GROUP_AGGREGATION_WEIGHTS
		// Gaussian decay
		float tmp;
		const float group_weight = std::exp(-total_variance/2.f/aggreSigma/aggreSigma);
#else
		// logistic decay
		float total_stddev = sqrtf(total_variance);
		const float group_weight = 1.f - 1.f/(1.f + 1e-4f + std::exp(-(total_stddev - aggreSigma)));
#endif

		for (int c = 0; c < io_group.size(); ++c)
		for (int n = 0; n < p_nSimP; ++n)
			aggreWeights[c][n] *= group_weight;

		// XXX this is for visualization
		total_variance = group_weight;
	}



	// return percentage of captured variance
//	return rank_variance / total_variance;
//	return entropy;
//	return total_weights;
//	return sqrtf(total_variance);
	return total_variance;
}
#endif

/**
 * @brief Implementation of computeBayesEstimateStep1 using an 
 * external basis provided by the user in the workspace i_mat.covEigVecs.
 *
 * See computeBayesEstimateStep1 for information about the arguments.
 **/
float computeBayesEstimateStep1_externalBasis(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	std::vector<std::vector<float> > &aggreWeights
){
	//! Parameters initialization
	const float beta_sigma = p_params.beta * p_params.sigma;
	const float beta_sigma2 = beta_sigma * beta_sigma;
#ifndef USE_BETA_FOR_VARIANCE
	const float sigma2 = p_params.sigma * p_params.sigma;
#else
	const float sigma2 = beta_sigma2;
#endif
	const unsigned sPC = p_params.sizePatch * p_params.sizePatch
	                   * p_params.sizePatchTime;
	const unsigned r   = sPC; // p_params.rank; // XXX FIXME TODO

	//! Variances
	float  rank_variance = 0.f;
	float total_variance = 0.f;

	for (unsigned c = 0; c < io_group.size(); c++)
	{
#ifndef DCT_DONT_CENTER1
		//! Center 3D group
		centerData(io_group[c], i_mat.baricenter, p_nSimP, sPC);
#endif

		if (r > 0)
		{
			/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
			 * row a patch. Thus, in column-major storage it corresponds to X^T, where
			 * each column of X contains a centered data point.
			 *
			 * We need to compute the noiseless estimage hX as 
			 * hX = U * W * U' * X
			 * where U is the matrix with the eigenvectors and W is a diagonal matrix
			 * with the filter coefficients.
			 *
			 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
			 * we compute 
			 * hX' = X' * U * (W * U')
			 */

			//! Project over basis: Z' = X'*U
			productMatrix(i_mat.groupTranspose,
			              io_group[c],
			              i_mat.patch_basis,
			              p_nSimP, sPC, sPC,
			              false, false);

			//! Compute variance over each component
			//  TODO: compute r leading components
			i_mat.covEigVals.resize(sPC);
			for (int k = 0; k < sPC; ++k)
			{
				float  comp_k_var = 0.f;
				float *comp_k = i_mat.groupTranspose.data() + k * p_nSimP;

				for (int i = 0; i < p_nSimP; ++i)
					comp_k_var += comp_k[i] * comp_k[i];

				i_mat.covEigVals[k] = comp_k_var / (float)p_nSimP;
				total_variance += i_mat.covEigVals[k];
			}

			//! Compute aggregation weights
			if (p_params.aggreGammaPatch)
			{
				aggreWeights[c].assign(p_nSimP, 0.f);

				//! Mahalanobis distance
				float tmp;
				for (int i = 0; i < sPC; ++i)
				{
					const float iEigVal = 1.f/std::max(i_mat.covEigVals[i],sigma2);
					for (int n = 0; n < p_nSimP; ++n)
						aggreWeights[c][n] += (tmp = i_mat.groupTranspose[i*p_nSimP + n])
						                   *  tmp * iEigVal ;
				}

				//! Exponentiate
				float igamma = 1.f/2.f/p_params.aggreGammaPatch;
				for(int n = 0; n < p_nSimP; ++n)
					aggreWeights[c][n] = std::exp(-igamma * aggreWeights[c][n]);
			}

			//! Compute eigenvalues-based coefficients of Bayes' filter
			for (unsigned k = 0; k < r; ++k)
			{
				rank_variance  += i_mat.covEigVals[k];
				float var = i_mat.covEigVals[k] - sigma2;

#if defined(THRESHOLD_WEIGHTS1)
				var = std::max(0.f, var);
#elif defined(LI_ZHANG_DAI1)
				var += sigma2; // add back sigma2
				var = (var < 4.f*sigma2) ? 0.f
				    : (var - 2.f*sigma2 + sqrtf(var*(var - 4.f*sigma2)))*.5f;
#endif

#ifndef LINEAR_THRESHOLDING1
				i_mat.covEigVals[k] = (fabs(var + beta_sigma2) > 1e-8f) ? 
					                   var / ( var + beta_sigma2 ) : 
											 0.f;
#else
				// this doesn't make sense for Li-Zhang-Dai's variance
				i_mat.covEigVals[k] = (var > 0.f) ? 1.f : 0.f;
#endif
			}

			//! U * W
			i_mat.covEigVecs.resize(sPC* sPC);
			float *eigVecs = i_mat.covEigVecs .data();
			float *basis   = i_mat.patch_basis.data();
			for (unsigned k = 0; k < r  ; ++k)
			for (unsigned i = 0; i < sPC; ++i)
				*eigVecs++ = *basis++ * i_mat.covEigVals[k];

			//! hX' = Z'*(U*W)'
			productMatrix(io_group[c],
			              i_mat.groupTranspose,
			              i_mat.covEigVecs,
			              p_nSimP, sPC, r,
			              false, true);

#ifndef DCT_DONT_CENTER1
			//! Add baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					io_group[c][k] += i_mat.baricenter[j];
#endif
		}
		else
		{
#ifndef DCT_DONT_CENTER1
			//! rank = 0: set as baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					io_group[c][k] = i_mat.baricenter[j];
#endif

			//! Avoid 0/0 in return statement
			total_variance = 1.f;
		}
	}

	// return percentage of captured variance
	return rank_variance / total_variance;
}

/**
 * @brief Implementation of computeBayesEstimateStep1 using an 
 * external basis provided by the user in the workspace i_mat.covEigVecs.
 * This version implements the Bayesian hyper-prior on the mean patch.
 *
 * See computeBayesEstimateStep1 for information about the arguments.
 **/
float computeBayesEstimateStep1_externalBasisHyper(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float sigma  = p_params.beta * p_params.sigma;
	const float sigma2 = sigma * sigma;

	const float sigmaM2   = p_params.sigma    * p_params.sigma;
	const float betaMAPM2 = p_params.betaMean * p_params.betaMean;
	const float betaVARM2 = 1; //betaMAPM2;

	const unsigned sPC = p_params.sizePatch * p_params.sizePatch
	                   * p_params.sizePatchTime;
	const unsigned r   = sPC;

	//! Variances
	float  rank_variance = 0.f;
	float total_variance = 0.f;

	for (unsigned c = 0; c < io_group.size(); c++)
	{
		/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
		 * row a patch. Thus, in column-major storage it corresponds to X^T, where
		 * each column of X contains a centered data point.
		 *
		 * We need to compute the noiseless estimage hX as 
		 * hX = U * W * U' * X
		 * where U is the matrix with the eigenvectors and W is a diagonal matrix
		 * with the filter coefficients.
		 *
		 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
		 * we compute 
		 * hX' = X' * U * (W * U')
		 */

		//! Project data over basis: Z' = X'*U
		productMatrix(i_mat.groupTranspose,
		              io_group[c],
		              i_mat.patch_basis,
		              p_nSimP, sPC, sPC,
		              false, false);

		//! Compute baricenter and center data
		centerData(i_mat.groupTranspose, i_mat.baricenter, p_nSimP, sPC);

#ifdef MEAN_HYPERPRIOR_BM3D1
		//! Store currrent baricenter before update
		i_mat.tmpMat = i_mat.baricenter;

		//! Filter baricenter
		for (int i = 0; i < sPC; ++i)
		{
			float var = std::max(0.f,
			            (float)p_nSimP*i_mat.baricenter[i]*i_mat.baricenter[i] -
			            betaVARM2*sigmaM2);
			i_mat.baricenter[i] *= var / (var + betaMAPM2*sigmaM2);
		}

		//! Re-center data at filtered baricenter
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				i_mat.groupTranspose[k] -= (i_mat.baricenter[j] - i_mat.tmpMat[j]);
#else
		const int MAX_BARI_ITERS = 10;
		float bari_delta;
		i_mat.covMat = i_mat.baricenter;
		for (int iter = 0; iter < MAX_BARI_ITERS; ++iter)
		{
			//! Store currrent baricenter before update
			i_mat.tmpMat = i_mat.baricenter;

			//! Compute variance over each component
			i_mat.covEigVals.resize(sPC);
			for (int k = 0; k < sPC; ++k)
			{
				float  comp_k_var = 0.f;
				float *comp_k = i_mat.groupTranspose.data() + k * p_nSimP;

				for (int i = 0; i < p_nSimP; ++i)
					comp_k_var += comp_k[i] * comp_k[i];

				i_mat.covEigVals[k] = comp_k_var / (float)p_nSimP;
				total_variance += i_mat.covEigVals[k];
			}

			//! Filter baricenter
			for (int i = 0; i < sPC; ++i)
			{
				float var = std::max(0.f,
				            (float)p_nSimP*i_mat.baricenter[i]*i_mat.baricenter[i] - 
				            betaVARM2 * (sigmaM2 + i_mat.covEigVals[i]));
				i_mat.baricenter[i] = i_mat.covMat[i] * var
				                    / (var + betaMAPM2 * (sigmaM2 + i_mat.covEigVals[i]));
			}

			//! Re-center data at filtered baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					i_mat.groupTranspose[k] -= (i_mat.baricenter[j] - i_mat.tmpMat[j]);

			bari_delta = 0.f;
			for (unsigned j = 0; j < sPC; j++)
				bari_delta += (i_mat.baricenter[j] - i_mat.tmpMat[j])
				            * (i_mat.baricenter[j] - i_mat.tmpMat[j]);
			
			printf("%g ", bari_delta);
		}

		printf("\n");
#endif
		//! Update variance over each component
		i_mat.covEigVals.resize(sPC);
		for (int k = 0; k < sPC; ++k)
		{
			float  comp_k_var = 0.f;
			float *comp_k = i_mat.groupTranspose.data() + k * p_nSimP;

			for (int i = 0; i < p_nSimP; ++i)
				comp_k_var += comp_k[i] * comp_k[i];

			i_mat.covEigVals[k] = comp_k_var / (float)p_nSimP;
			total_variance += i_mat.covEigVals[k];
		}


		//! Compute filter coefficients
		for (int i = 0; i < r; ++i)
		{
			float var = i_mat.covEigVals[i] - sigma2;
#ifdef THRESHOLD_WEIGHTS1
			var = std::max(0.f, var);
#endif
#ifndef LINEAR_THRESHOLDING1
			i_mat.covEigVals[i] = var / ( var + sigma2 );
#else
			i_mat.covEigVals[i] = (var > 0.f) ? 1.f : 0.f;
#endif
			rank_variance += var;
		}

		for (int i = r; i < sPC; ++i)
			i_mat.covEigVals[i] = 0.f;

		//! Z' * W
		float *comp = i_mat.groupTranspose.data();
		for (unsigned k = 0; k < sPC; ++k)
			for (unsigned i = 0; i < p_nSimP; ++i)
				*comp++ *= i_mat.covEigVals[k];


		//! Add baricenter
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				i_mat.groupTranspose[k] += i_mat.baricenter[j];

		//! hX' = Z'*W*U'
		productMatrix(io_group[c],
		              i_mat.groupTranspose,
		              i_mat.patch_basis,
		              p_nSimP, sPC, r,
		              false, true);
	}

	// return percentage of captured variance
	return rank_variance / total_variance;
}

/**
 * @brief Implementation of computeBayesEstimateStep1 using an 
 * external basis provided by the user in the workspace i_mat.covEigVecs.
 * This version considers different non-linear thresholding operators.
 *
 * See computeBayesEstimateStep1 for information about the arguments.
 **/
float computeBayesEstimateStep1_externalBasisTh(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float beta_sigma = p_params.beta * p_params.sigma;
	const float beta_sigma2 = beta_sigma * beta_sigma;
#ifndef USE_BETA_FOR_VARIANCE
	const float sigma2 = p_params.sigma * p_params.sigma;
#else
	const float sigma2 = beta_sigma2;
#endif
	const unsigned sPC = p_params.sizePatch * p_params.sizePatch
	                   * p_params.sizePatchTime;
	const unsigned r   = sPC; // p_params.rank; // XXX FIXME TODO

	//! Variances
	float  rank_variance = 0.f;
	float total_variance = 0.f;

	for (unsigned c = 0; c < io_group.size(); c++)
	{
		if (r > 0)
		{
			/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
			 * row a patch. Thus, in column-major storage it corresponds to X^T, where
			 * each column of X contains a centered data point.
			 *
			 * We need to compute the noiseless estimage hX as 
			 * hX = U * W * U' * X
			 * where U is the matrix with the eigenvectors and W is a diagonal matrix
			 * with the filter coefficients.
			 *
			 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
			 * we compute 
			 * hX' = X' * U * (W * U')
			 */

			//! Project over basis: Z' = X'*U
			productMatrix(i_mat.groupTranspose,
			              io_group[c],
			              i_mat.patch_basis,
			              p_nSimP, sPC, sPC,
			              false, false);

#ifndef DCT_DONT_CENTER1
			//! Center 3D group
			centerData(i_mat.groupTranspose, i_mat.baricenter, p_nSimP, sPC);
#endif

#if defined(SOFT_THRESHOLD1_BAYES)
			//! Compute variance over each component
			i_mat.covEigVals.resize(sPC);
			for (int k = 0; k < sPC; ++k)
			{
				float  comp_k_var = 0.f;
				float *comp_k = i_mat.groupTranspose.data() + k * p_nSimP;

				for (int i = 0; i < p_nSimP; ++i)
					comp_k_var += comp_k[i] * comp_k[i];

				i_mat.covEigVals[k] = comp_k_var / (float)p_nSimP;
				total_variance += i_mat.covEigVals[k];
			}

			//! Substract sigma2 and compute variance captured by the r leading eigenvectors
			for (int i = 0; i < r; ++i)
			{
				float tmp = std::max(i_mat.covEigVals[i] - sigma2, 0.f);
				i_mat.covEigVals[i] = (tmp < 1e-6) ? FLT_MAX
					                 : sqrt(2.f) * beta_sigma2 / sqrt(tmp);
				rank_variance += tmp;
			}
#endif

			//! Thresholding
			float *z = i_mat.groupTranspose.data();
			for (unsigned k = 0; k < r  ; ++k)
			for (unsigned i = 0; i < p_nSimP; ++i, ++z)
#if defined(SOFT_THRESHOLD1)
				*z = *z > 0 ? std::max(*z - beta_sigma, 0.f)
				            : std::min(*z + beta_sigma, 0.f);
#elif defined(SOFT_THRESHOLD1_BAYES)
				*z = *z > 0 ? std::max(*z - i_mat.covEigVals[k], 0.f)
				            : std::min(*z + i_mat.covEigVals[k], 0.f);
#else
				*z = (*z * *z) > beta_sigma2 ? *z : 0.f;
#endif

#ifdef MEAN_HYPERPRIOR_BM3D1
			z = i_mat.baricenter.data();
			for (unsigned k = 0; k < r  ; ++k, ++z)// if (k != 0)
 #if defined(SOFT_THRESHOLD1)
				*z = *z > 0 ? std::max(*z - beta_sigma/sqrtf((float)p_nSimP), 0.f)
				            : std::min(*z + beta_sigma/sqrtf((float)p_nSimP), 0.f);
 #elif defined(SOFT_THRESHOLD1_BAYES)
				*z = *z > 0 ? std::max(*z - i_mat.covEigVals[k]/sqrtf((float)p_nSimP), 0.f)
				            : std::min(*z + i_mat.covEigVals[k]/sqrtf((float)p_nSimP), 0.f);
 #else
				*z = (*z * *z) > 0.1*beta_sigma2/(float)p_nSimP ? *z : 0.f;
 #endif
#endif

#ifndef DCT_DONT_CENTER1
			//! Add baricenter
			z = i_mat.groupTranspose.data();
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++, ++z)
					*z += i_mat.baricenter[j];
#endif

			//! hX' = Z'*U'
			productMatrix(io_group[c],
			              i_mat.groupTranspose,
//			              i_mat.patch_basis,
			              i_mat.patch_basis_inv,
			              p_nSimP, sPC, r,
			              false, true);
		}
		else
		{
			//! Center 3D group
			centerData(io_group[c], i_mat.baricenter, p_nSimP, sPC);

			//! rank = 0: set as baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					io_group[c][k] = i_mat.baricenter[j];

			//! Avoid 0/0 in return statement
			total_variance = 1.f;
		}
	}

	// return percentage of captured variance
	return rank_variance / total_variance;
}

/**
 * @brief Compute the Bayes estimation assuming a low rank covariance matrix.
 *
 * @param io_group: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *		- groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 *		- baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 *		- covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *		- covMatTmp: allocated memory. Used to process the Bayes estimate;
 *		- tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 * @param aggreWeights: output aggregation weights.
 *
 * @return none.
 **/
float computeBayesEstimateStep1(
	std::vector<std::vector<float> > &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	std::vector<std::vector<float> > &aggreWeights 
){
	/* In some images might have a dark noiseless frame (for instance the image
	 * might be the result of a transformation, and a black region correspond to
	 * the pixels that where before the frame before the registration). A group
	 * of patches with 0 variance create problems afterwards. Such groups have to
	 * be avoided.
	 * */

	/* TODO
	 * This fix only addresses the problem in which all patches are equal. A more
	 * subtle situation is when some patches are equal. These patches shouldn't be
	 * denoised, since it is extremely unlekely that such an event occurs.
	 * */

	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime;

	for (unsigned c = 0; c < io_group.size(); c++)
	for (unsigned j = 0; j < sPC; j++)
	{
#ifdef USE_FFTW
		float v = io_group[c][j];
		for (unsigned i = 1; i < p_nSimP; i++)
			if (v != io_group[c][j + i * sPC])
				goto not_equal;
#else
		float v = io_group[c][j * p_nSimP];
		for (unsigned i = 1; i < p_nSimP; i++)
			if (v != io_group[c][j * p_nSimP + i])
				goto not_equal;
#endif
	}

	//! All patches are equal ~ do nothing
	return 0.f;

not_equal:
	//! Not all patches are equal ~ denoise
	return
#if (defined(THRESHOLDING1))
		computeBayesEstimateStep1_externalBasisTh(io_group, i_mat,
			io_nInverseFailed, p_params, p_nSimP);
#elif (defined(MEAN_HYPERPRIOR1))
		computeBayesEstimateStep1_externalBasisHyper(io_group, i_mat,
			io_nInverseFailed, p_params, p_nSimP);
#else
  #ifdef USE_FFTW
		computeBayesEstimateStep1_externalBasisFFTW(io_group, i_mat,
			io_nInverseFailed, p_params, p_nSimP, aggreWeights);
  #else
		computeBayesEstimateStep1_externalBasis(io_group, i_mat,
			io_nInverseFailed, p_params, p_nSimP, aggreWeights);
  #endif
#endif

}

#ifdef USE_FFTW
/**
 * @brief Implementation of computeBayesEstimateStep2 using an 
 * external basis provided by the user in the workspace i_mat.covEigVecs.
 *
 * See computeBayesEstimateStep2 for information about the arguments.
 **/
float computeBayesEstimateStep2_externalBasisFFTW(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	std::vector<float>  &aggreWeights
){
	//! Parameters initialization
	const float beta_sigma = p_params.beta * p_params.sigma;
	const float beta_sigma2 = beta_sigma * beta_sigma;
#ifndef USE_BETA_FOR_VARIANCE
	const float sigma2 = p_params.sigma * p_params.sigma;
#else
	const float sigma2 = beta_sigma2;
#endif
	const unsigned sPC = p_params.sizePatch * p_params.sizePatch
	                   * p_params.sizePatchTime;
	const unsigned r   = sPC; // p_params.rank; // XXX FIXME TODO

	const float inSimP = 1.f/(float)p_nSimP;

	//! Variances
	float r_variance = 0.f;
	float total_variance = 0.f;
	aggreWeights.assign(p_nSimP, p_params.aggreGammaPatch ? 0.f : 1.f);

	for (unsigned c = 0; c < p_imSize.channels; c++)
	{
		std::vector<float> groupNoisy_c(io_groupNoisy.begin() + sPC*p_nSimP * c   ,
		                                io_groupNoisy.begin() + sPC*p_nSimP *(c+1));

		std::vector<float> groupBasic_c( i_groupBasic.begin() + sPC*p_nSimP * c   ,
		                                 i_groupBasic.begin() + sPC*p_nSimP *(c+1));

#ifndef DCT_DONT_CENTER2
		//! Center 3D groups around their baricenter
 #ifdef BARICENTER_BASIC
		i_mat.baricenter.assign(sPC, 0.f);
		for (int i = 0; i < p_nSimP; ++i)
		for (int k = 0; k < sPC; ++k)
			i_mat.baricenter[k] += groupBasic_c[i*sPC + k] * inSimP;

		for (int i = 0; i < p_nSimP; ++i)
		for (int k = 0; k < sPC; ++k)
		{
			groupBasic_c[i*sPC + k] -= i_mat.baricenter[k];
			groupNoisy_c[i*sPC + k] -= i_mat.baricenter[k];
		}

 #else //BARICENTER_BASIC
		i_mat.baricenter.assign(sPC, 0.f);
		for (int i = 0; i < p_nSimP; ++i)
		for (int k = 0; k < sPC; ++k)
			i_mat.baricenter[k] += groupBasic_c[i*sPC + k] * inSimP;

		for (int i = 0; i < p_nSimP; ++i)
		for (int k = 0; k < sPC; ++k)
			groupBasic_c[i*sPC + k] -= i_mat.baricenter[k];

		i_mat.baricenter.assign(sPC, 0.f);
		for (int i = 0; i < p_nSimP; ++i)
		for (int k = 0; k < sPC; ++k)
			i_mat.baricenter[k] += groupNoisy_c[i*sPC + k] * inSimP;

		for (int i = 0; i < p_nSimP; ++i)
		for (int k = 0; k < sPC; ++k)
			groupNoisy_c[i*sPC + k] -= i_mat.baricenter[k];
 #endif//BARICENTER_BASIC
#endif//DCT_DONT_CENTER2

		if (r > 0)
		{
			//! Forward DCT
			globalDCT.forward(groupBasic_c);
			globalDCT.forward(groupNoisy_c);

			//! Compute variance over each component
			i_mat.covEigVals.assign(sPC, 0.f);
			for (int i = 0; i < p_nSimP; ++i)
			for (int k = 0; k < sPC; ++k)
				i_mat.covEigVals[k] += groupBasic_c[i*sPC + k]
				                     * groupBasic_c[i*sPC + k];

			for (int k = 0; k < sPC; ++k)
			{
				i_mat.covEigVals[k] *= inSimP;
				total_variance += i_mat.covEigVals[k]/(float)sPC/(float)p_imSize.channels;
			}

			//! Compute Mahalanobis distance for patch aggregation weights
			if (p_params.aggreGammaPatch)
			{
				float tmp;
				for (int n = 0; n < p_nSimP; ++n)
				for (int i = 0; i < sPC    ; ++i)
						aggreWeights[n] += (tmp = groupBasic_c[n*sPC + i])
						                 *  tmp / std::max(i_mat.covEigVals[i],sigma2);
			}

			//! Compute eigenvalues-based coefficients of Bayes' filter
			for (unsigned k = 0; k < r; ++k)
			{
				if (i_mat.covEigVals[k] == 0) continue;

				r_variance += i_mat.covEigVals[k];
				float var = i_mat.covEigVals[k];

#ifdef NOISY_COVARIANCE2
 #if defined(THRESHOLD_WEIGHTS2)
				var -= std::min(var, sigma2);
 #elif defined(LI_ZHANG_DAI2)
				var = (var < 4.f*sigma2) ? 0.f
				    : (var - 2.f*sigma2 + sqrtf(var*(var - 4.f*sigma2)))*.5f;
 #endif
#endif

#if defined(LINEAR_HARD_THRESHOLDING2)
				i_mat.covEigVals[k] = var > beta_sigma2 ? 1.f : 0.f;
#elif defined(LINEAR_SOFT_THRESHOLDING2)
				i_mat.covEigVals[k] = var > beta_sigma2
				                    ? 1.f - beta_sigma2/var : 0.f;
#else
				i_mat.covEigVals[k] = var / ( var + beta_sigma2);
#endif
			}

			//! W*Z
			for (int i = 0; i < p_nSimP; ++i)
			for (int k = 0; k < sPC; ++k)
				groupNoisy_c[i*sPC + k] *= i_mat.covEigVals[k];

			//! Inverse DCT
			globalDCT.inverse(groupNoisy_c);

#ifndef DCT_DONT_CENTER2
			//! Add baricenter
			for (int i = 0; i < p_nSimP; ++i)
			for (int k = 0; k < sPC; ++k)
				groupNoisy_c[i*sPC + k] += i_mat.baricenter[k];
#endif
		}
		else
		{
#ifndef DCT_DONT_CENTER2
			//! rank = 0: set as baricenter
			for (int i = 0; i < p_nSimP; ++i)
			for (int k = 0; k < sPC; ++k)
				groupNoisy_c[i*sPC + k] = i_mat.baricenter[k];
#endif
			//! Avoid 0/0 in return statement
			total_variance = 1.f;
		}

		//! Copy channel back into vector
		std::copy(groupNoisy_c.begin(), groupNoisy_c.end(),
		          io_groupNoisy.begin() + sPC*p_nSimP*c);
	}


	//! Exponentiate Mahalanobis distance for patch aggregation weights
	if (p_params.aggreGammaPatch)
	{
		float igamma = 1.f/2.f/p_params.aggreGammaPatch;
		for(int n = 0; n < p_nSimP; ++n)
			aggreWeights[n] = std::exp(-igamma * aggreWeights[n]);
	}

	// compute group aggregation weights
	if (p_params.aggreGammaGroup)
	{
		const float aggreSigma = p_params.aggreGammaGroup*p_params.sigma;

#ifdef GAUSSIAN_GROUP_AGGREGATION_WEIGHTS
		// Gaussian decay
		float tmp;
		const float group_weight = std::exp(-total_variance/2.f/aggreSigma/aggreSigma);
#else
		// logistic decay
		float total_stddev = sqrtf(total_variance);
		const float group_weight = 1.f - 1.f/(1.f + 1e-4f + std::exp(-(total_stddev - aggreSigma)));
#endif

		for (int n = 0; n < p_nSimP; ++n)
			aggreWeights[n] *= group_weight;

		// XXX this is for visualization
		total_variance = group_weight;
	}

	// return percentage of captured variance
	return total_variance;
}
#endif

/**
 * @brief Implementation of computeBayesEstimateStep2 computing the
 * principal directions of the a priori covariance matrix. This functions
 * computes the eigenvectors/values of the data covariance matrix using LAPACK.
 *
 * See computeBayesEstimateStep2 for information about the arguments.
 **/
float computeBayesEstimateStep2_externalBasis(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	std::vector<float>  &aggreWeights
){
	//! Parameters initialization
	const float beta_sigma = p_params.beta * p_params.sigma;
	const float beta_sigma2 = beta_sigma * beta_sigma;
#ifndef USE_BETA_FOR_VARIANCE
	const float sigma2 = p_params.sigma * p_params.sigma;
#else
	const float sigma2 = beta_sigma2;
#endif
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime;
	const unsigned r    = sPC; // p_params.rank; // XXX FIXME TODO

	float r_variance = 0.f;
	float total_variance = 1.f;
//	float total_weights = 0.f;
	aggreWeights.assign(p_nSimP, p_params.aggreGammaPatch ? 0.f : 1.f);

	for (unsigned c = 0; c < p_imSize.channels; c++)
	{
		std::vector<float> groupNoisy_c(io_groupNoisy.begin() + sPC*p_nSimP * c   ,
		                                io_groupNoisy.begin() + sPC*p_nSimP *(c+1));

		std::vector<float> groupBasic_c( i_groupBasic.begin() + sPC*p_nSimP * c   ,
		                                 i_groupBasic.begin() + sPC*p_nSimP *(c+1));

		if (r > 0)
		{
			/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
			 * row a patch. Thus, in column-major storage it corresponds to X^T, where
			 * each column of X contains a centered data point.
			 *
			 * We need to compute the noiseless estimage hX as 
			 * hX = U * W * U' * X
			 * where U is the matrix with the eigenvectors and W is a diagonal matrix
			 * with the filter coefficients.
			 *
			 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
			 * we compute 
			 * hX' = X' * U * (W * U')
			 */

			//! Project noisy patches over basis: Z' = X'*U (to compute variances)
			productMatrix(i_mat.groupTransposeNoisy,
			              groupNoisy_c,
			              i_mat.patch_basis,
			              p_nSimP, sPC, sPC,
			              false, false);

			//! Project basic patches over basis: Z' = X'*U (to compute variances)
#ifndef NOISY_COVARIANCE2
			productMatrix(i_mat.groupTranspose,
			              groupBasic_c,
			              i_mat.patch_basis,
			              p_nSimP, sPC, sPC,
			              false, false);
#else
			i_mat.groupTranspose = i_mat.groupTransposeNoisy;
#endif

#ifndef DCT_DONT_CENTER2
			//! Center 3D groups around their baricenter
 #ifdef BARICENTER_BASIC

			//! Center basic and noisy patches using basic's baricenter
			centerData(i_mat.groupTranspose, i_mat.baricenter, p_nSimP, sPC);
			for (unsigned j = 0, k = 0; j < sPC; j++)
			{
				i_mat.baricenterNoisy[j] = i_mat.baricenter[j];
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					i_mat.groupTransposeNoisy[k] -= i_mat.baricenter[j];
			}

 #else //BARICENTER_BASIC

			//! Center basic noisy patches each with its own baricenter
			centerData(i_mat.groupTranspose     , i_mat.baricenter     , p_nSimP, sPC);
			centerData(i_mat.groupTransposeNoisy, i_mat.baricenterNoisy, p_nSimP, sPC);

 #endif//BARICENTER_BASIC
#endif//DCT_DONT_CENTER2


			//! Compute variance over each component
			//  TODO: compute r leading components
			i_mat.covEigVals.resize(sPC);
			for (int k = 0; k < sPC; ++k)
			{
				float  comp_k_var = 0.f;
				float *comp_k = i_mat.groupTranspose.data() + k * p_nSimP;

				for (int i = 0; i < p_nSimP; ++i)
					comp_k_var += comp_k[i] * comp_k[i];

				i_mat.covEigVals[k] = comp_k_var / (float)p_nSimP;
				total_variance += i_mat.covEigVals[k]/(float)sPC/(float)p_imSize.channels;
			}

			//! Compute Mahalanobis distance for patch aggregation weights
			if (p_params.aggreGammaPatch)
			{
				float tmp;
				for (int i = 0; i < sPC; ++i)
				{
					const float iEigVal = 1.f/i_mat.covEigVals[i];
					for (int n = 0; n < p_nSimP; ++n)
						aggreWeights[n] += (tmp = i_mat.groupTranspose[i*p_nSimP + n])
						                 *  tmp * iEigVal ;
				}
			}

			//! Compute eigenvalues-based coefficients of Bayes' filter
			for (int i = 0; i < r; ++i)
			{
				r_variance  += i_mat.covEigVals[i];
				float var = i_mat.covEigVals[i];

#ifdef NOISY_COVARIANCE2
 #if defined(THRESHOLD_WEIGHTS2)
				var -= std::min(var, sigma2);
 #elif defined(LI_ZHANG_DAI2)
				var = (var < 4.f*sigma2) ? 0.f
				    : (var - 2.f*sigma2 + sqrtf(var*(var - 4.f*sigma2)))*.5f;
 #endif
#endif

#if defined(LINEAR_HARD_THRESHOLDING2)
				i_mat.covEigVals[i] = var > beta_sigma2 ? 1.f : 0.f;
#elif defined(LINEAR_SOFT_THRESHOLDING2)
				i_mat.covEigVals[i] = var > beta_sigma2
				                    ? 1.f - beta_sigma2/var : 0.f;
#else
				i_mat.covEigVals[i] = var / ( var + beta_sigma2);
#endif
			}

//			for (int k = 0; k < sPC; ++k)
//				total_weights += i_mat.covEigVals[k]/(float)sPC/(float)p_imSize.channels;

			//! U * W
			i_mat.covEigVecs.resize(sPC*sPC);
			float *basis   = i_mat.patch_basis.data();
			float *eigVecs = i_mat.covEigVecs .data();
			for (unsigned k = 0; k < r  ; ++k)
			for (unsigned i = 0; i < sPC; ++i)
				*eigVecs++ = *basis++ * i_mat.covEigVals[k];

			//! hX' = Z'*(U*W)'
			productMatrix(groupNoisy_c,
			              i_mat.groupTransposeNoisy,
			              i_mat.covEigVecs,
			              p_nSimP, sPC, r,
			              false, true);

#ifndef DCT_DONT_CENTER2
			//! Add baricenter

#ifdef MEAN_HYPERPRIOR_BM3D2
			float *z = i_mat.baricenterNoisy.data();
			float *y = i_mat.baricenter.data();
			for (unsigned k = 0; k < r  ; ++k, ++z, ++y)// if (k != 0)
				*z *= (*y**y)/(*y**y + beta_sigma2/(float)p_nSimP);
#endif

			//! invert dct
			productMatrix(i_mat.baricenter,
			              i_mat.baricenterNoisy,
			              i_mat.patch_basis,
			              1, sPC, r,
			              false, true);

			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					groupNoisy_c[k] += i_mat.baricenter[j];
#endif

//	printMatrix(i_mat.groupTranspose, sPC, p_nSimP, "/tmp/z.asc");
//	printMatrix(io_groupNoisy       , sPC, p_nSimP, "/tmp/x_filtered.asc");
//
//	if (1)
//	{
//		// stop here
//		printf("stopped: PRESS Ctrl-C\n");
//		while (1) int a = 1;
//	}
		}
		else
		{
			//! r = 0: set all patches as baricenter
 #ifdef BARICENTER_BASIC
			centerData(groupBasic_c, i_mat.baricenter, p_nSimP, sPC);
 #else
			centerData(groupNoisy_c, i_mat.baricenter, p_nSimP, sPC);
 #endif
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					groupNoisy_c[k] = i_mat.baricenter[j];
		}

		//! Copy channel back into vector
		std::copy(groupNoisy_c.begin(), groupNoisy_c.end(),
		          io_groupNoisy.begin() + sPC*p_nSimP*c);

	}

	//! Exponentiate Mahalanobis distance
	if (p_params.aggreGammaPatch)
	{
		float igamma = 1.f/2.f/p_params.aggreGammaPatch;
		for(int n = 0; n < p_nSimP; ++n)
			aggreWeights[n] = std::exp(-igamma * aggreWeights[n]);
	}

	// compute group aggregation weights
	if (p_params.aggreGammaGroup)
	{
		const float aggreSigma = p_params.aggreGammaGroup*p_params.sigma;

#ifdef GAUSSIAN_GROUP_AGGREGATION_WEIGHTS
		// Gaussian decay
		float tmp;
		const float group_weight = std::exp(-total_variance/2.f/aggreSigma/aggreSigma);
#else
		// logistic decay
		float total_stddev = sqrtf(total_variance);
		const float group_weight = 1.f - 1.f/(1.f + std::exp(-(total_stddev - aggreSigma)));
#endif

		for (int n = 0; n < p_nSimP; ++n)
			aggreWeights[n] *= group_weight;

		// XXX this is for visualization
		total_variance = group_weight;
	}

	// return percentage of captured variance
//	return r_variance / total_variance;
//	return total_weights;
//	return sqrtf(total_variance);
	return total_variance;
}

/**
 * @brief Implementation of computeBayesEstimateStep2 using an 
 * external basis provided by the user in the workspace i_mat.covEigVecs.
 * This version implements the Bayesian hyper-prior on the mean patch.
 *
 * See computeBayesEstimateStep2 for information about the arguments.
 **/
float computeBayesEstimateStep2_externalBasisHyper(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float sigma  = p_params.beta * p_params.sigma;
	const float sigma2 = sigma * sigma;

	const float sigmaM2 = p_params.sigma * p_params.sigma;
	const float betaM2  = p_params.betaMean * p_params.betaMean;

	const unsigned sPC = p_params.sizePatch * p_params.sizePatch
	                   * p_params.sizePatchTime;
	const unsigned r   = sPC; // p_params.rank; // XXX FIXME TODO

	//! Variances
	float r_variance = 0.f;
	float total_variance = 1.f;

	for (unsigned c = 0; c < p_imSize.channels; c++)
	{
		std::vector<float> groupNoisy_c(io_groupNoisy.begin() + sPC*p_nSimP * c   ,
		                                io_groupNoisy.begin() + sPC*p_nSimP *(c+1));

		std::vector<float> groupBasic_c( i_groupBasic.begin() + sPC*p_nSimP * c   ,
		                                 i_groupBasic.begin() + sPC*p_nSimP *(c+1));

		/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
		 * row a patch. Thus, in column-major storage it corresponds to X^T, where
		 * each column of X contains a centered data point.
		 *
		 * We need to compute the noiseless estimage hX as 
		 * hX = U * W * U' * X
		 * where U is the matrix with the eigenvectors and W is a diagonal matrix
		 * with the filter coefficients.
		 *
		 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
		 * we compute 
		 * hX' = X' * U * (W * U')
		 */

		//! Learn priors from basic patches -------------------------------

		//! Project basic data over basis: Z' = X'*U
		productMatrix(i_mat.groupTranspose,
		              groupBasic_c,
		              i_mat.patch_basis,
		              p_nSimP, sPC, sPC,
		              false, false);

		//! Compute baricenter and center data
		centerData(i_mat.groupTranspose, i_mat.baricenter, p_nSimP, sPC);

		//! Compute prior variance over each component
		i_mat.covEigVals.resize(sPC);
		for (int k = 0; k < sPC; ++k)
		{
			float  comp_k_var = 0.f;
			float *comp_k = i_mat.groupTranspose.data() + k * p_nSimP;

			for (int i = 0; i < p_nSimP; ++i)
				comp_k_var += comp_k[i] * comp_k[i];

			i_mat.covEigVals[k] = comp_k_var / (float)p_nSimP;
			total_variance += i_mat.covEigVals[k];
		}

		//! Compute prior variances for baricenter
		i_mat.covMat.resize(sPC);
		for (int k = 0; k < sPC; ++k)
#ifndef MEAN_HYPERPRIOR_BM3D2
			i_mat.covMat[k] = std::max(0.f,
			      (float)p_nSimP * i_mat.baricenter[k] * i_mat.baricenter[k]
			      - i_mat.covEigVals[k]);
#else
			i_mat.covMat[k] =
			      (float)p_nSimP * i_mat.baricenter[k] * i_mat.baricenter[k];
#endif

		//! Filter noisypatches -------------------------------------------

		//! Project noisy data over basis: Z' = X'*U
		productMatrix(i_mat.groupTranspose,
		              groupNoisy_c,
		              i_mat.patch_basis,
		              p_nSimP, sPC, sPC,
		              false, false);

		//! Compute baricenter and center data
		centerData(i_mat.groupTranspose, i_mat.baricenter, p_nSimP, sPC);

		//! Store currrent baricenter before update
		i_mat.tmpMat = i_mat.baricenter;

		//! Filter baricenter
		for (int i = 1; i < sPC; ++i) //NOTE: we dont't filter baricenters DC!!
#ifndef MEAN_HYPERPRIOR_BM3D2
			i_mat.baricenter[i] *= (i_mat.covMat[i] != 0) ? i_mat.covMat[i] /
			        (i_mat.covMat[i] + betaM2 * (sigmaM2 + i_mat.covEigVals[i])) : 0;
#else
			i_mat.baricenter[i] /= 1 + betaM2 * sigmaM2 / i_mat.covMat[i];
#endif

		//! Re-center data at filtered baricenter
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				i_mat.groupTranspose[k] -= (i_mat.baricenter[j] - i_mat.tmpMat[j]);

		//! Compute filter coefficients
		for (int i = 0; i < r; ++i)
		{
			float var = i_mat.covEigVals[i];
#if defined(LINEAR_HARD_THRESHOLDING2)
			i_mat.covEigVals[k] = (var > sigma2) ? 1.f : 0.f;
#elif defined(LINEAR_SOFT_THRESHOLDING2)
			i_mat.covEigVals[k] = (var > sigma2) ? 1.f - sigma2/var : 0.f;
#else
			i_mat.covEigVals[i] = var / ( var + sigma2 );
#endif
			r_variance += var;
		}

		for (int i = r; i < sPC; ++i)
			i_mat.covEigVals[i] = 0.f;

		//! Z' * W
		float *comp = i_mat.groupTranspose.data();
		for (unsigned k = 0; k < sPC; ++k)
			for (unsigned i = 0; i < p_nSimP; ++i)
				*comp++ *= i_mat.covEigVals[k];

		//! Add baricenter
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				i_mat.groupTranspose[k] += i_mat.baricenter[j];

		//! hX' = Z'*W*U'
		productMatrix(groupNoisy_c,
		              i_mat.groupTranspose,
		              i_mat.patch_basis,
		              p_nSimP, sPC, r,
		              false, true);

		//! Copy channel back into vector
		std::copy(groupNoisy_c.begin(), groupNoisy_c.end(),
		          io_groupNoisy.begin() + sPC*p_nSimP*c);
	}

	// return percentage of captured variance
	return r_variance / total_variance;
}

/**
 * @brief Implementation of computeBayesEstimateStep2 computing the
 * principal directions of the a priori covariance matrix. This functions
 * computes the eigenvectors/values of the data covariance matrix using LAPACK.
 *
 * See computeBayesEstimateStep2 for information about the arguments.
 **/
float computeBayesEstimateStep2_externalBasisTh(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
){
	//! Parameters initialization
	const float sigma2 = p_params.beta * p_params.sigma * p_params.sigma;
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime;
	const unsigned r    = sPC; // p_params.rank; // XXX FIXME TODO

	float r_variance = 0.f;
	float total_variance = 1.f;

	for (unsigned c = 0; c < p_imSize.channels; c++)
	{
		std::vector<float> groupNoisy_c(io_groupNoisy.begin() + sPC*p_nSimP * c   ,
		                                io_groupNoisy.begin() + sPC*p_nSimP *(c+1));

		std::vector<float> groupBasic_c( i_groupBasic.begin() + sPC*p_nSimP * c   ,
		                                 i_groupBasic.begin() + sPC*p_nSimP *(c+1));

#ifndef DCT_DONT_CENTER2
		//! Center 3D groups around their baricenter
 #ifdef BARICENTER_BASIC

		//! Center basic and noisy patches using basic's baricenter
		centerData(groupBasic_c, i_mat.baricenter, p_nSimP, sPC);
		for (unsigned j = 0, k = 0; j < sPC; j++)
			for (unsigned i = 0; i < p_nSimP; i++, k++)
				groupNoisy_c[k] -= i_mat.baricenter[j];

 #else //BARICENTER_BASIC

		//! Center basic noisy patches each with its own baricenter
		centerData(groupNoisy_c, i_mat.baricenter, p_nSimP, sPC);
		centerData(groupBasic_c, i_mat.baricenter, p_nSimP, sPC);

 #endif//BARICENTER_BASIC
#endif//DCT_DONT_CENTER2

		if (r > 0)
		{
			/* NOTE: io_groupNoisy, if read as a column-major matrix, contains in each
			 * row a patch. Thus, in column-major storage it corresponds to X^T, where
			 * each column of X contains a centered data point.
			 *
			 * We need to compute the noiseless estimage hX as 
			 * hX = U * W * U' * X
			 * where U is the matrix with the eigenvectors and W is a diagonal matrix
			 * with the filter coefficients.
			 *
			 * Matrix U is stored (column-major) in i_mat.covEigVecs. Since we have X^T
			 * we compute 
			 * hX' = X' * U * (W * U')
			 */

			//! Project basic patches over basis: Z' = X'*U (to compute variances)
#ifndef NOISY_COVARIANCE2
			productMatrix(i_mat.groupTranspose,
			              groupBasic_c,
			              i_mat.patch_basis,
			              p_nSimP, sPC, sPC,
			              false, false);
#else
			productMatrix(i_mat.groupTranspose,
			              groupNoisy_c,
			              i_mat.patch_basis,
			              p_nSimP, sPC, sPC,
			              false, false);
#endif

#ifdef SOFT_THRESHOLD2_BAYES
			//! Compute variance over each component
			//  TODO: compute r leading components
			i_mat.covEigVals.resize(sPC);
			for (int k = 0; k < sPC; ++k)
			{
				float  comp_k_var = 0.f;
				float *comp_k = i_mat.groupTranspose.data() + k * p_nSimP;

				for (int i = 0; i < p_nSimP; ++i)
					comp_k_var += comp_k[i] * comp_k[i];

				i_mat.covEigVals[k] = comp_k_var / (float)p_nSimP;
				total_variance += i_mat.covEigVals[k];
			}

			//! Substract sigma2 and compute variance captured by the r leading eigenvectors
			for (int i = 0; i < r; ++i)
			{
 #ifdef NOISY_COVARIANCE2
				float tmp = std::max(i_mat.covEigVals[i] - sigma2, 0.f);
 #else
				float tmp = i_mat.covEigVals[i];
 #endif
				i_mat.covEigVals[i] = (tmp > 1e-6) ? sqrt(2.f) * sigma2 / sqrt(tmp) : FLT_MAX;
			}
#endif

			//! Thresholding
			float *z = i_mat.groupTranspose.data();
			for (unsigned k = 0; k < r  ; ++k)
			for (unsigned i = 0; i < p_nSimP; ++i, ++z)
 #if defined(SOFT_THRESHOLD2)
				*z = *z > 0 ? std::max(*z - sigma, 0.f) : std::min(*z + sigma, 0.f);
 #elif defined(SOFT_THRESHOLD2_BAYES)
				*z = *z > 0 ? std::max(*z - i_mat.covEigVals[k], 0.f)
				            : std::min(*z + i_mat.covEigVals[k], 0.f);
 #else
				*z = (*z * *z) > sigma2 ? *z : 0.f;
 #endif

			//! hX' = Z'*U'
			productMatrix(groupNoisy_c,
			              i_mat.groupTranspose,
			              i_mat.patch_basis,
			              p_nSimP, sPC, r,
			              false, true);

#ifndef DCT_DONT_CENTER2
			//! Add baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					groupNoisy_c[k] += i_mat.baricenter[j];
#endif
		}
		else
			//! r = 0: set all patches as baricenter
			for (unsigned j = 0, k = 0; j < sPC; j++)
				for (unsigned i = 0; i < p_nSimP; i++, k++)
					groupNoisy_c[k] = i_mat.baricenter[j];

		//! Copy channel back into vector
		std::copy(groupNoisy_c.begin(), groupNoisy_c.end(),
		          io_groupNoisy.begin() + sPC*p_nSimP*c);

	}


	// return percentage of captured variance
	return r_variance / total_variance;

}

/**
 * @brief Compute the Bayes estimation assuming a low rank covariance matrix.
 *
 * @param io_groupNoisy: inputs all similar patches in the noisy image,
 *                         outputs their denoised estimates.
 * @param i_grougroup contains all similar patches in the basic image.
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
 * @param aggreWeights: output aggregation weights.
 *
 * @return estimate of kept variance.
 **/
float computeBayesEstimateStep2(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_size
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	std::vector<float> &aggreWeights
){
	/* In some images might have a dark noiseless frame (for instance the image
	 * might be the result of a transformation, and a black region correspond to
	 * the pixels that where before the frame before the registration). A group
	 * of patches with 0 variance create problems afterwards. Such groups have to
	 * be avoided.
	 * */

	/* TODO
	 * This fix only addresses the problem in which all patches are equal. A more
	 * subtle situation is when some patches are equal. These patches shouldn't be
	 * denoised, since it is extremely unlekely that such an event occurs.
	 * */

	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime * p_size.channels;

	for (unsigned c = 0; c < io_groupNoisy.size(); c++)
	for (unsigned j = 0; j < sPC; j++)
	{
		float v = io_groupNoisy[j * p_nSimP];
		for (unsigned i = 1; i < p_nSimP; i++)
			if (v != io_groupNoisy[j * p_nSimP + i])
				goto not_equal;
	}

	//! All patches are equal ~ do nothing
	return 0.f;

not_equal:
	//! Not all patches are equal ~ denoise
	return
#if defined(THRESHOLDING2)
		computeBayesEstimateStep2_externalBasisTh(io_groupNoisy, i_groupBasic,
			i_mat, io_nInverseFailed, p_size, p_params, p_nSimP);
#elif defined(MEAN_HYPERPRIOR2)
		computeBayesEstimateStep2_externalBasisHyper(io_groupNoisy, i_groupBasic,
			i_mat, io_nInverseFailed, p_size, p_params, p_nSimP);
#else
  #ifdef USE_FFTW
		computeBayesEstimateStep2_externalBasisFFTW(io_groupNoisy, i_groupBasic, i_mat,
			io_nInverseFailed, p_size, p_params, p_nSimP, aggreWeights);
  #else
		computeBayesEstimateStep2_externalBasis(io_groupNoisy, i_groupBasic, i_mat,
			io_nInverseFailed, p_size, p_params, p_nSimP, aggreWeights);
  #endif
#endif

}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param aggreWeights: input aggregation weights.
 * @param i_index: contains index of all similar patches contained in i_group;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 *
 * @return masked: number of processable pixels that were flaged non-processable.
 **/
int computeAggregationStep1(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<std::vector<float> > const& i_group
,	std::vector<std::vector<float> > const& aggreWeights
,	std::vector<float> const& aggreWindow
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
	const unsigned sPC = sPx*sPx*sPt;

	int masked = 0;

#ifndef DEBUG_SHOW_PATCH_GROUPS
	//! Aggregate estimates
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / whc) * wh + ind % wh;
		for (unsigned pt = 0, i = 0; pt < sPt; pt++)
		for (unsigned py = 0       ; py < sPx; py++)
		for (unsigned px = 0       ; px < sPx; px++, i++)
		{
			for (unsigned c = 0; c < chnls; c++)
			{
				const unsigned ij = ind  + c * wh;
				// XXX NOTE: we only use aggregWeights[0] to avoid color artifacts XXX
#ifdef USE_FFTW
				io_im(ij + pt * whc + py * w + px) += 
					aggreWeights[0][n] * aggreWindow[i] * i_group[c][n*sPC + i];
#else
				io_im(ij + pt * whc + py * w + px) += 
					aggreWeights[0][n] * aggreWindow[i] * i_group[c][i * p_nSimP + n];
#endif
			}
			io_weight(ind1 + pt * wh + py * w + px) += aggreWeights[0][n] * aggreWindow[i];
		}

		//! Use Paste Trick
		unsigned px, py, pt;
		io_mask.sz.coords(ind1, px, py, pt);

		// TODO Modify this part so that if the weight are small enough, the
		//      patch is not deactivated
		if (p_params.doPasteBoost)
		{
			if (io_mask(ind1)) masked++;
			io_mask(ind1) = false;

			if ((py >     2*sPx) && io_mask(ind1 - w)) masked++;
			if ((py < h - 2*sPx) && io_mask(ind1 + w)) masked++;
			if ((px >     2*sPx) && io_mask(ind1 - 1)) masked++;
			if ((px < w - 2*sPx) && io_mask(ind1 + 1)) masked++;

			if (py >     2*sPx) io_mask(ind1 - w) = false;
			if (py < h - 2*sPx) io_mask(ind1 + w) = false;
			if (px >     2*sPx) io_mask(ind1 - 1) = false;
			if (px < w - 2*sPx) io_mask(ind1 + 1) = false;
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
					aggreWeights[c][n] * i_group[c][(pt * sPx*sPx + py * sPx + px) * p_nSimP + n];
			}
			io_weight(ind1 + pt * wh + py * w + px) += aggreWeights[c][n];
		}
	}

	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		const unsigned ind1 = (ind / whc) * wh + ind % wh;

		//! Use Paste Trick
		unsigned px, py, pt;
		io_mask.sz.coords(ind1, px, py, pt);

		if (p_params.doPasteBoost)
		{
			if (io_mask(ind1)) masked++;
			io_mask(ind1) = false;

			if ((py >     2*sPx) && io_mask(ind1 - w)) masked++;
			if ((py < h - 2*sPx) && io_mask(ind1 + w)) masked++;
			if ((px >     2*sPx) && io_mask(ind1 - 1)) masked++;
			if ((px < w - 2*sPx) && io_mask(ind1 + 1)) masked++;

			if (py >     2*sPx) io_mask(ind1 - w) = false;
			if (py < h - 2*sPx) io_mask(ind1 + w) = false;
			if (px >     2*sPx) io_mask(ind1 - 1) = false;
			if (px < w - 2*sPx) io_mask(ind1 + 1) = false;
		}
	}
#endif

	if (!p_params.doPasteBoost) 
	{
		masked++;
		io_mask(i_index[0]) = false;
	}

	return masked;
}

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the image with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param aggreWeights: input aggregation weights.
 * @param i_index: contains index of all similar patches contained in i_group;
 * @param p_imSize: size of io_im;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return masked: number of processable pixels that were flaged non-processable.
 **/
int computeAggregationStep2(
	Video<float> &io_im
,	Video<float> &io_weight
,	Video<char>  &io_mask
,	std::vector<float> const& i_group
,	std::vector<float> const& aggreWeights
,	std::vector<float> const& aggreWindow
,	Video<float> &variance
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
	const unsigned sPC = sPx*sPx*sPt;

	unsigned nAgg = p_nSimP;

	//! Aggregate estimates
	int masked = 0;
	for (unsigned n = 0; n < nAgg; n++)
	{
#ifdef USE_FFTW
		const unsigned ind = i_index[n];
		for (unsigned c = 0; c < chnls; c++)
		{
			const unsigned ij = ind + c * wh;
			for (unsigned pt = 0, k = 0; pt < sPt; pt++)
			for (unsigned py = 0       ; py < sPx; py++)
			for (unsigned px = 0       ; px < sPx; px++, k++)
				io_im(ij + pt * whc + py * w + px) +=
					aggreWeights[n] * aggreWindow[i] * i_group[k + (n + c * p_nSimP) * sPC];
#else
		const unsigned ind = i_index[n];
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			const unsigned ij = ind + c * wh;
			for (unsigned pt = 0, i = 0; pt < sPt; pt++)
			for (unsigned py = 0       ; py < sPx; py++)
			for (unsigned px = 0       ; px < sPx; px++, k++, i++)
				io_im(ij + pt * whc + py * w + px) +=
					aggreWeights[n] * aggreWindow[i] * i_group[k * p_nSimP + n];
#endif
		}

		const unsigned ind1 = (ind / whc) * wh + ind % wh;
		for (unsigned pt = 0, i = 0; pt < sPt; pt++)
		for (unsigned py = 0       ; py < sPx; py++)
		for (unsigned px = 0       ; px < sPx; px++, i++)
			io_weight(ind1 + pt * wh + py * w + px) += aggreWeights[n] * aggreWindow[i];

		//! Apply Paste Trick
		unsigned px, py, pt;
		io_mask.sz.coords(ind1, px, py, pt);

		// TODO Modify this part so that if the weight are small enough, the
		//      patch is not deactivated
		if (p_params.doPasteBoost)
		{
			if (io_mask(ind1)) masked++;
			io_mask(ind1) = false;

			if ((py >     2*sPx) && io_mask(ind1 - w)) masked++;
			if ((py < h - 2*sPx) && io_mask(ind1 + w)) masked++;
			if ((px >     2*sPx) && io_mask(ind1 - 1)) masked++;
			if ((px < w - 2*sPx) && io_mask(ind1 + 1)) masked++;

			if (py >     2*sPx) io_mask(ind1 - w) = false;
			if (py < h - 2*sPx) io_mask(ind1 + w) = false;
			if (px >     2*sPx) io_mask(ind1 - 1) = false;
			if (px < w - 2*sPx) io_mask(ind1 + 1) = false;

#ifdef DEBUG_SHOW_WEIGHT
			if (py >     2*sPx) variance(ind1 - w) = variance(ind1);
			if (py < h - 2*sPx) variance(ind1 + w) = variance(ind1);
			if (px >     2*sPx) variance(ind1 - 1) = variance(ind1);
			if (px < w - 2*sPx) variance(ind1 + 1) = variance(ind1);
#endif
		}
	}

	if (!p_params.doPasteBoost) 
	{
		masked++;
		io_mask(i_index[0]) = false;
	}

	return masked;
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
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group;
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
,	std::vector<std::vector<float> > const& i_group
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
				io_im(ij + p * width + q) += i_group[c][(p * sP + q) * p_nSimP + n];
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
				io_im(ij + p * width + q) += i_group[c][(p * sP + q) * p_nSimP + n];
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
 * @param i_group: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group;
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
,	std::vector<float> const& i_group
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
				io_im(ij + p * width + q) += i_group[k * p_nSimP + n];
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
