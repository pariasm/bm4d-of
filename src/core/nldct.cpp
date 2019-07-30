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



/* The parameter beta is used as a noise correction factor. It provides a way 
 * to control the thresholding/filtering strength of the algorithm, by modifying
 * sigma as beta * sigma. By default, this modifications should only be applied 
 * to the sigma in the threshold/filter operator. If the following flag is defined,
 * the modified sigma is also used to estimate the a priori variances in the first
 * step. This option causes a 0.05dB increase in the results for a particular value
 * of beta. However, the resultng method seems much more sensitive to beta. 
 */
//#define USE_BETA_FOR_VARIANCE

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
	nlbParams &params
,	const unsigned step
,	const float sigma
,	const VideoSize &size
,	const bool verbose
){
	const bool s1 = (step == 1);

	if (size.frames == 1)
	{
		// default parameters for BM3D
		if (sigma <= 40) // normal profile
		{
			params.sizePatch = 8; 
			params.nSimilarPatches = s1 ? 16 : 32;
			params.offSet = 3;
			params.beta = s1 ? 2.7 : 1.; // FIXME: beta here means the threshold
			params.tau = s1 ? 54.8 : 20.;
		}
		else // vn profile
		{
			params.sizePatch = s1 ? 8 : 11; 
			params.nSimilarPatches = 32;
			params.offSet = s1 ? 4 : 6;
			params.beta = s1 ? 2.8 : 1.;
			params.tau = s1 ? 158. : 59.2;
		}

		params.sizePatchTime = 1;
		params.sizeSearchTimeRangeFwd = 0;
		params.sizeSearchTimeRangeBwd = 0;
		params.sizeSearchWindowPred = 0;
		params.offSetTime = 1;
		params.dsub = 0;
		params.nSimilarPatchesPred = params.nSimilarPatches;
		params.sizeSearchWindow = 39;
	}
	else
	{
		// default parameters for VBM3D (only high noise profile)
		params.transform = s1 ? bior1_5 : dct;
		params.sizePatch = (s1 || sigma >= 30) ? 8 : 7; 
		params.sizePatchTime = 1;
		params.nSimilarPatches = 8;
		params.sizeSearchWindow = 7;
		params.sizeSearchWindowPred = 5;
		params.sizeSearchTimeRangeFwd = 4;
		params.sizeSearchTimeRangeBwd = 4;
		params.nSimilarPatchesPred = 2;
		params.offSet = s1 ? 6 : 4;
		params.offSetTime = 1;
		params.beta = s1 ? 2.7 : 1.;
		params.tau = s1 ? 67.1 : 54.8; // TODO this is only true for high noise
		params.dsub = s1 ? 7 : 3;
	}

	// these params are common for vbm3d and bm3d
	params.orderInvariance = false;
	params.isFirstStep = s1;
	params.sigma = sigma;
	params.verbose = verbose;
	params.colorSpace = YUV;
	params.doPasteBoost = false;
	params.agg_window = true;
	params.boundary = 2*(params.sizeSearchWindow/2) + (params.sizePatch - 1);
	params.transform = s1 ? bior1_5 : dct;
}

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

	if (prms.sizePatch != 8) prms.transform = dct;

	if (size.frames > 1)
		prms.offSet = (prms.isFirstStep) ? (float)sizePatch/4.*3. : sizePatch/2;
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

	if (prms.sizeSearchTimeRangeFwd + prms.sizeSearchTimeRangeFwd == 0)
		prms.nSimilarPatchesPred = prms.nSimilarPatches;

}

/**
 * @brief Display parameters of the NL-Bayes algorithm.
 *
 * @param i_params : nlbParams for first or second step of the algorithm;
 *
 * @return none.
 **/
void printNlbParameters(
	const nlbParams &p
){
	int px = p.sizePatch, pt = p.sizePatchTime;
	int wtb = p.sizeSearchTimeRangeBwd, wtf = p.sizeSearchTimeRangeBwd;
	int wx  = p.sizeSearchWindow;

	printf("\x1b[37;01m" "Parameters for step %d:" ANSI_RST "\n" , p.isFirstStep ? 1 : 2);
	printf("\tPatch search:\n");
	printf("\t\tPatch size                  = %dx%dx%d\n" , px, px, pt);
	printf("\t\tNumber of patches           = %d\n"       , p.nSimilarPatches);
	printf("\t\tSpatial search window       = %dx%d\n"    , wx, wx);
	printf("\t\tTemporal search window      = [-%d,%d]\n" , wtb, wtf);
	printf("\t\tDistance threshold (tau)    = %g\n"       , p.tau);
#ifdef VBM3D_SEARCH
	int wxpred = p.sizeSearchWindowPred;
	printf("\t\tPred. spatial search window = %dx%d\n"    , wxpred, wxpred);
	printf("\t\tPred. distance bias (dsub)  = %g\n"       , p.dsub);
	printf("\t\tSimilar patches per frame   = %d\n"       , p.nSimilarPatchesPred);
#endif
	printf("\tGroup filtering:\n");
	printf("\t\tTransform                   = %s\n"       , p.transform == dct ? "dct" : "bior1.5");
	printf("\t\tBeta                        = %g\n"       , p.beta);
	printf("\t\tOrder invariance            = %s\n"       , p.orderInvariance ? "yes" : "no");
	printf("\tAggregation:\n");
	printf("\t\tOffset                      = %d\n"       , p.offSet);
	printf("\t\tPasteBoost                  = %s\n"       , p.doPasteBoost
			? "active" : "inactive");
	printf("\t\tAggregation window          = %s\n"       , p.agg_window ? "yes" : "no");
	printf("\n");
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
,	Video<float> const& i_fflow
,	Video<float> const& i_bflow
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
#ifdef VBM3D_SEARCH
		printf(ANSI_BCYN "VBM3D_SEARCH > Using VBM3D predictive patch search\n" ANSI_RST);
#endif
#ifdef VBM3D_HAAR_TRANSFORM
		printf(ANSI_BCYN "VBM3D_HAAR_TRANSFORM > Using VBM3D 3D transform\n" ANSI_RST);
#endif
#ifdef USE_BETA_FOR_VARIANCE
		printf(ANSI_BCYN "USE_BETA_FOR_VARIANCE > Noise correction in step 1 applied both MAP and variances.\n" ANSI_RST);
#endif
	}


	//! Number of available cores
	int nThreads = 1;
#ifdef _OPENMP
	nThreads = omp_get_max_threads();
	if (p_prms1.verbose) printf(ANSI_CYN "OpenMP is using %d threads\n" ANSI_RST, nThreads);
#endif
	const int nParts = 2 * nThreads;

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
				processNlBayes(imNoisySub[n], fflowSub[n], bflowSub[n],
				               imBasicSub[n], imFinalSub[n], p_prms1, imCrops[n]);

		//! Get the basic estimate
		VideoUtils::subBuildTight(imBasicSub, o_imBasic, p_prms1.boundary);

		//! YUV to RGB
		if (p_prms1.colorSpace == YUV)
			VideoUtils::transformColorSpace(o_imBasic, false);

		for (int n = 0; n < (int)nParts; n++)
			groupsRatio[0] += 100.f * (float)groupsProcessedSub[n]/(float)imSize.whf;

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
		for (int n = 0; n < (int)nParts; n++)
			groupsProcessedSub[n] =
				processNlBayes(imNoisySub[n], fflowSub[n], bflowSub[n],
				               imBasicSub[n], imFinalSub[n], p_prms2, imCrops[n]);

		//! Get the final result
		VideoUtils::subBuildTight(imFinalSub, o_imFinal, p_prms2.boundary);

		//! Undo color transform
		if (p_prms2.colorSpace == YUV)
		{
			VideoUtils::transformColorSpace(o_imBasic, false);
			VideoUtils::transformColorSpace(o_imFinal, false);
		}

		for (int n = 0; n < (int)nParts; n++)
			groupsRatio[1] += 100.f * (float)groupsProcessedSub[n]/(float)imSize.whf;

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
		if ( (dy % stepy == 0) ||
		     (!border_y1 && y == end_y - 1) ||
		     (!border_y0 && y == ori_y    ) )
		if ( (dx % stepx == 0) ||
		     (!border_x1 && x == end_x - 1) ||
		     (!border_x0 && x == ori_x    ) )
		{
			mask(x,y,f) = true;
			n_groups++;
		}
	}


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
		if (p_params.agg_window && sPx == 8) // use BM3D's Kaiser window
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
		else if (p_params.agg_window && sPx == 10)
		{
			float *z = mat.agg_window.data();

			for (unsigned t = 0; t < sPt; t++)
			{
				*z++ = 0.1924; *z++ = 0.2763; *z++ = 0.3503; *z++ = 0.4055; *z++ = 0.4349; *z++ = 0.4349; *z++ = 0.4055; *z++ = 0.3503; *z++ = 0.2763; *z++ = 0.1924;		  
				*z++ = 0.2763; *z++ = 0.3967; *z++ = 0.5030; *z++ = 0.5822; *z++ = 0.6245; *z++ = 0.6245; *z++ = 0.5822; *z++ = 0.5030; *z++ = 0.3967; *z++ = 0.2763;		  
				*z++ = 0.3503; *z++ = 0.5030; *z++ = 0.6377; *z++ = 0.7381; *z++ = 0.7917; *z++ = 0.7917; *z++ = 0.7381; *z++ = 0.6377; *z++ = 0.5030; *z++ = 0.3503;		  
				*z++ = 0.4055; *z++ = 0.5822; *z++ = 0.7381; *z++ = 0.8544; *z++ = 0.9164; *z++ = 0.9164; *z++ = 0.8544; *z++ = 0.7381; *z++ = 0.5822; *z++ = 0.4055;		  
				*z++ = 0.4349; *z++ = 0.6245; *z++ = 0.7917; *z++ = 0.9164; *z++ = 0.9829; *z++ = 0.9829; *z++ = 0.9164; *z++ = 0.7917; *z++ = 0.6245; *z++ = 0.4349;		  
				*z++ = 0.4349; *z++ = 0.6245; *z++ = 0.7917; *z++ = 0.9164; *z++ = 0.9829; *z++ = 0.9829; *z++ = 0.9164; *z++ = 0.7917; *z++ = 0.6245; *z++ = 0.4349;		  
				*z++ = 0.4055; *z++ = 0.5822; *z++ = 0.7381; *z++ = 0.8544; *z++ = 0.9164; *z++ = 0.9164; *z++ = 0.8544; *z++ = 0.7381; *z++ = 0.5822; *z++ = 0.4055;		  
				*z++ = 0.3503; *z++ = 0.5030; *z++ = 0.6377; *z++ = 0.7381; *z++ = 0.7917; *z++ = 0.7917; *z++ = 0.7381; *z++ = 0.6377; *z++ = 0.5030; *z++ = 0.3503;		  
				*z++ = 0.2763; *z++ = 0.3967; *z++ = 0.5030; *z++ = 0.5822; *z++ = 0.6245; *z++ = 0.6245; *z++ = 0.5822; *z++ = 0.5030; *z++ = 0.3967; *z++ = 0.2763;		  
				*z++ = 0.1924; *z++ = 0.2763; *z++ = 0.3503; *z++ = 0.4055; *z++ = 0.4349; *z++ = 0.4349; *z++ = 0.4055; *z++ = 0.3503; *z++ = 0.2763; *z++ = 0.1924;		  
			}
		}
		else
			for (unsigned i = 0; i < sPt*sPx*sPx; i++)
				mat.agg_window[i] = 1.;
	}


	//! Variance captured by the principal components
	Video<float> variance(mask.sz);


	//! Total number of groups of similar patches processed
	unsigned group_counter = 0;

	if (step1)
	{
		//! Allocate Sizes
		io_imBasic.resize(sz);

		//! Matrices used for Bayes' estimate
		vector<float> groupNoisy(            patch_num * patch_dim * sz.channels);
		vector<float> groupBasic(step1 ? 0 : patch_num * patch_dim * sz.channels);
		vector<float> aggreWeights(patch_num, 1.f);

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
				unsigned nSimP = estimateSimilarPatches(i_imNoisy, io_imBasic,
						i_fflow, i_bflow, groupNoisy, groupBasic, index, ij3, p_params);

				//! Bayes' estimate
				variance(ij) = computeBayesEstimateStep1(groupNoisy, mat,
						nInverseFailed, sz, p_params, nSimP, aggreWeights);

				//! Aggregation
				remaining_groups -=
					computeAggregationStep1(io_imBasic, weight, mask, groupNoisy,
							aggreWeights, mat.agg_window, index, p_params, nSimP);
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
				unsigned nSimP = estimateSimilarPatches(i_imNoisy, io_imBasic,
						i_fflow, i_bflow, groupNoisy, groupBasic, index, ij3, p_params);

				//! Bayes' estimate
				variance(ij) = computeBayesEstimateStep2(groupNoisy, groupBasic,
						mat, nInverseFailed, sz, p_params, nSimP, aggreWeights);

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

	return group_counter;
}

// Comparison function for sorting
bool compareFirst(const std::vector<float> &a, const std::vector<float> &b){
	return a[0] < b[0];
}

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
unsigned estimateSimilarPatches(
	Video<float> const& imNoisy
,	Video<float> const& imBasic
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<float> &groupNoisy
,	std::vector<float> &groupBasic
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
){
	bool step1 = params.isFirstStep;
	const VideoSize sz = imNoisy.sz;
	const int sPx = params.sizePatch;
	const int sPt = params.sizePatchTime;
	const int dist_chnls = step1 ? 1 : sz.channels;
	unsigned nSimP = 1;

	if (params.nSimilarPatches > 1)
	{
		const unsigned nSimPFrame = params.nSimilarPatchesPred;
		const float dsub = params.dsub * params.dsub * 255;
		const float tau_match = params.tau * params.tau * sPx * sPx * sPt;

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
				sWx = params.sizeSearchWindowPred;
				sWy = params.sizeSearchWindowPred;
				for (int i = 0; i < nSimPFrame; i++)
				{
					int cx = traj_cx[i][dt - dir];
					int cy = traj_cy[i][dt - dir];
					int ct = traj_ct[i][dt - dir];
					float cxf = cx + (use_flow ? (dir > 0 ? fflow(cx,cy,ct,0) : bflow(cx,cy,ct,0)) : 0.f);
					float cyf = cy + (use_flow ? (dir > 0 ? fflow(cx,cy,ct,1) : bflow(cx,cy,ct,1)) : 0.f);
					cx0.push_back(std::max(0.f, std::min(sz.width  - 1.f, round(cxf))));
					cy0.push_back(std::max(0.f, std::min(sz.height - 1.f, round(cyf))));
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
//
				// Pointer to video used for distance computation
				const Video<float> *p_im = step1 ? &imNoisy : &imBasic;

				// add elements of the square search region to the union
				for (int qy = rangey[0]; qy <= rangey[1]; qy++)
				for (int qx = rangex[0]; qx <= rangex[1]; qx++)
				if (search_region.insert(sz.index(qx, qy, qt, 0)).second)
				{
					//! Squared L2 distance
					float dist = 0.f, dif;
					for (int c = 0; c < dist_chnls; c++)
					for (int ht = 0; ht < sPt; ht++)
					for (int hy = 0; hy < sPx; hy++)
					for (int hx = 0; hx < sPx; hx++)
						dist += (dif = (*p_im)(px + hx, py + hy, pt + ht, c)
						             - (*p_im)(qx + hx, qy + hy, qt + ht, c)) * dif;

					//! Small bias towards the center of the region
					dist -= (qy == py && qx == px)? dsub : 0;

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

		//! Remove patches with a distance larger than tau_match
		{
			int n = 1;
			while (2*n <= nSimP && distance[n].first < tau_match) n *= 2;
			nSimP = n;
		}

		//! Store indices of most similar patches
		for (unsigned n = 0; n < nSimP; n++)
			index[n] = distance[n].second;
	}
	else // nSimilarPatches == 1
		index[0] = pidx;

	//! Stack selected patches into the 3D group
	const unsigned w   = sz.width;
	const unsigned wh  = sz.wh;
	const unsigned whc = sz.whc;
	for (unsigned c = 0, k = 0; c < sz.channels; c++)
	for (unsigned ht = 0; ht < sPt; ht++)
	for (unsigned hy = 0; hy < sPx; hy++)
	for (unsigned hx = 0; hx < sPx; hx++)
	for (unsigned n = 0 ; n < nSimP; n++, k++)
	{
		groupNoisy[k] = imNoisy(c * wh + index[n] + ht * whc + hy * w + hx);
		if (!step1)
			groupBasic[k] = imBasic(c * wh + index[n] + ht * whc + hy * w + hx);
	}

	return nSimP;
}
#else
 #ifdef MC_PATCHES
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
unsigned estimateSimilarPatches(
	Video<float> const& imNoisy
,	Video<float> const& imBasic
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<float> &groupNoisy
,	std::vector<float> &groupBasic
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
){
	bool mc_patches = true;

	bool step1 = params.isFirstStep;
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime;
	const VideoSize sz = imNoisy.sz;
	const int dist_chnls = step1 ? 1 : sz.channels;
	unsigned nSimP = 1;

	//! Coordinates of center of search box
	unsigned px, py, pt, pc;
	sz.coords(pidx, px, py, pt, pc);

	if (params.nSimilarPatches > 1)
	{
		const float tau = params.tau * params.tau * sPx * sPx * sPt;
		bool use_flow = (fflow.sz.width > 0);

		//! Determine search range
		int sWx = params.sizeSearchWindow;
		int sWy = params.sizeSearchWindow;
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
		using std::vector;
		vector<vector<float> > distance(sWx*sWy*sWt, vector<float>(1 + 3*sPt));

		//! Number of patches in search region
		int nsrch = 0;

		//! Schedule frames in temporal range
		std::vector<int> srch_ranget;
		srch_ranget.push_back(pt);
		for (int qt = pt+1; qt <= ranget[1]; ++qt) srch_ranget.push_back(qt);
		for (int qt = pt-1; qt >= ranget[0]; --qt) srch_ranget.push_back(qt);

		// Pointer to video used for distance computation
		const Video<float> *p_im = step1 ? &imNoisy : &imBasic;

		//! Store reference patch
		std::vector<float> refpatch(sPx*sPx*sPt*dist_chnls, 0.f);
		{
			//! Initial temporal slice
			int k = 0;
			for (int hy = 0; hy < sPx; hy++)
			for (int hx = 0; hx < sPx; hx++)
			for (int c  = 0; c < dist_chnls; ++c, k++)
				refpatch[k] = (*p_im)(px + hx, py + hy, pt, c);

			//! Remaninig temporal slices
			float cx = px + sPx/2, cy = py + sPx/2;
			int cx0 = round(cx), cy0 = round(cy);
			for (int ht = 1; ht < sPt; ht++)
			{
				//! Integrate optical flow to new patch center
				cx += (mc_patches ? fflow(cx0, cy0, pt + ht - 1, 0) : 0.f);
				cy += (mc_patches ? fflow(cx0, cy0, pt + ht - 1, 1) : 0.f);
				cx = std::max(float(sPx/2), std::min(float(sz.width  + sPx/2 - sPx), cx));
				cy = std::max(float(sPx/2), std::min(float(sz.height + sPx/2 - sPx), cy));

				cx0 = round(cx); cy0 = round(cy);
				for (int hy = 0; hy < sPx; hy++)
				for (int hx = 0; hx < sPx; hx++)
				for (int c  = 0; c < dist_chnls; ++c, k++)
					refpatch[k] = (*p_im)(cx0 - sPx/2 + hx, cy0 - sPx/2 + hy, pt + ht, c);
			}
		}

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
				cx[dt] = std::max(0.f, std::min(sz.width  - 1.f, round(cx_f)));
				cy[dt] = std::max(0.f, std::min(sz.height - 1.f, round(cy_f)));
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
				float dist = 0.f, dif;

				//! Initial temporal slice
				std::vector<float>::const_iterator p_refpatch = refpatch.begin();
				for (int hy = 0; hy < sPx; hy++)
				for (int hx = 0; hx < sPx; hx++)
				for (int c = 0; c < dist_chnls; c++, ++p_refpatch)
					dist += (dif = *p_refpatch - (*p_im)(qx + hx, qy + hy, qt, c)) * dif;

				//! patch trajectory and patch distance
				std::vector<float> traj(1 + 3 * sPt);

				// Initial point in trajectory
				traj[1] = qx; // coordinates are stored with an offset
				traj[2] = qy; // since traj[0] is for the patch distance
				traj[3] = qt;

				//! Remaninig temporal slices
				float cx = qx + sPx/2, cy = qy + sPx/2;
				int cx0 = qx, cy0 = qy;
				for (int ht = 1; ht < sPt; ht++)
				{
					if (mc_patches)
					{
						//! Integrate optical flow to new patch center
						cx += fflow(cx0 + sPx/2, cy0 + sPx/2, qt + ht - 1, 0);
						cy += fflow(cx0 + sPx/2, cy0 + sPx/2, qt + ht - 1, 1);
						cx = std::max(float(sPx/2), std::min(float(sz.width  + sPx/2 - sPx), cx));
						cy = std::max(float(sPx/2), std::min(float(sz.height + sPx/2 - sPx), cy));
						cx0 = round(cx) - sPx/2; cy0 = round(cy) - sPx/2;
					}
					for (int hy = 0; hy < sPx; hy++)
					for (int hx = 0; hx < sPx; hx++)
					for (int c = 0; c < dist_chnls; c++, ++p_refpatch)
					{
//						if (cx0 - sPx/2 + hx >= sz.width ||
//							 cy0 - sPx/2 + hy >= sz.height ||
//							 qt + ht          >= sz.frames)
//							printf("c0 = [%d,%d,%d] - h=[%d,%d,%d] - [%d,%d,%d]\n",
//									cx0, cy0, qt, hx,hy,ht, cx0 - sPx/2 + hx, cy0 - sPx/2 + hy, qt + ht);
						dist += (dif = *p_refpatch - (*p_im)(cx0 + hx, cy0 + hy, qt + ht)) * dif;
					}

					traj[ht*3 + 1] = cx0;
					traj[ht*3 + 2] = cy0;
					traj[ht*3 + 3] = qt + ht;
				}

				//! Save distance and corresponding patch index
				traj[0] = dist;
				distance[nsrch++] = traj;
			}
		}

		distance.resize(nsrch);

		//! Keep only the nSimilarPatches best similar patches
		nSimP = std::min(params.nSimilarPatches, (unsigned)distance.size());
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
		                  distance.end(), compareFirst);

		//! Remove patches with a distance larger than tau_match
		{
			int n = 1;
			while (2*n <= nSimP && distance[n][0] < tau) n *= 2;
			nSimP = n;
		}

		//! Store trajectories of most similar patches
		index.resize(nSimP*3*sPt);
		for (int n = 0; n < nSimP; ++n)
		for (int i = 0; i < sPt*3; ++i)
			// shift 1 distance vector since distance[n][0] is the distance value
			index[n*sPt*3 + i] = (unsigned)distance[n][i+1];
	}
	else // nSimilarPatches == 1
	{
		//! Initial temporal slice
		float cx = px + sPx/2, cy = py + sPx/2;
		int cx0 = round(cx), cy0 = round(cy);
		int k = 0;
		index[k++] = cx0;
		index[k++] = cy0;
		index[k++] = pt;
		//! Remaninig temporal slices
		for (int ht = 1; ht < sPt; ht++)
		{
			//! Integrate optical flow to new patch center
			cx += (mc_patches ? fflow(cx0, cy0, pt + ht - 1, 0) : 0.f);
			cy += (mc_patches ? fflow(cx0, cy0, pt + ht - 1, 1) : 0.f);
			cx = std::max(float(sPx/2), std::min(float(sz.width  - sPx/2 + sPx), cx));
			cy = std::max(float(sPx/2), std::min(float(sz.height - sPx/2 + sPx), cy));

			cx0 = round(cx); cy0 = round(cy);
			index[k++] = cx0 - sPx/2;
			index[k++] = cy0 - sPx/2;
			index[k++] = pt;
		}
	}

	//! Stack similar patches into 3D groups
	const unsigned w   = sz.width;
	const unsigned wh  = sz.wh;
	const unsigned whc = sz.whc;

	for (unsigned c  = 0; c < sz.channels; c++)
	for (unsigned ht = 0, k = 0; ht < sPt; ht++)
	for (unsigned hy = 0;        hy < sPx; hy++)
	for (unsigned hx = 0;        hx < sPx; hx++)
	for (unsigned n  = 0; n < nSimP; n++, k++)
	{
		int offset = 3*sPt*n + 3*ht;
		groupNoisy[k] = imNoisy(index[offset + 0] + hx,
		                        index[offset + 1] + hy,
		                        index[offset + 2], c);
		if (!step1)
			groupBasic[k] = imBasic(index[offset + 0] + hx,
			                        index[offset + 1] + hy,
			                        index[offset + 2], c);
	}

	return nSimP;
}
 #else //MC_PATCHES
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
unsigned estimateSimilarPatches(
	Video<float> const& imNoisy
,	Video<float> const& imBasic
,	Video<float> const& fflow
,	Video<float> const& bflow
,	std::vector<float> &groupNoisy
,	std::vector<float> &groupBasic
,	std::vector<unsigned> &index
,	const unsigned pidx
,	const nlbParams &params
){
	bool step1 = params.isFirstStep;
	const int sPx   = params.sizePatch;
	const int sPt   = params.sizePatchTime;
	const VideoSize sz = imNoisy.sz;
	const int chnls = sz.channels;
	const int dist_chnls = step1 ? 1 : sz.channels;
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
				float cx_f = cx0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,0) : bflow(cx0,cy0,ct0,0)) : 0.f);
				float cy_f = cy0 + (use_flow ? (dir > 0 ? fflow(cx0,cy0,ct0,1) : bflow(cx0,cy0,ct0,1)) : 0.f);
				cx[dt] = std::max(0.f, std::min(sz.width  - 1.f, round(cx_f)));
				cy[dt] = std::max(0.f, std::min(sz.height - 1.f, round(cy_f)));
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

			// Check if basic estimate was provided
			const Video<float> *p_im = step1 ? &imNoisy : &imBasic;

			//! Compute distance between patches in search range
			for (int qy = rangey[0], dy = 0; qy <= rangey[1]; qy++, dy++)
			for (int qx = rangex[0], dx = 0; qx <= rangex[1]; qx++, dx++)
			{
				//! Squared L2 distance
				float dist = 0.f, dif;
				for (int c = 0; c < dist_chnls; c++)
				for (int ht = 0; ht < sPt; ht++)
				for (int hy = 0; hy < sPx; hy++)
				for (int hx = 0; hx < sPx; hx++)
					dist += (dif = (*p_im)(px + hx, py + hy, pt + ht, c)
					             - (*p_im)(qx + hx, qy + hy, qt + ht, c) ) * dif;

				//! Save distance and corresponding patch index
				distance[nsrch++] = std::make_pair(dist, sz.index(qx, qy, qt, 0));
			}
		}

		distance.resize(nsrch);

		//! Keep only the nSimilarPatches best similar patches
		nSimP = std::min(params.nSimilarPatches, (unsigned)distance.size());
		std::partial_sort(distance.begin(), distance.begin() + nSimP,
		                  distance.end(), comparaisonFirst);

		//! Remove patches with a distance larger than tau_match
		{
			const float tau = params.tau * params.tau * sPx * sPx * sPt;
			int n = 1;
			while (2*n <= nSimP && distance[n].first < tau) n *= 2;
			nSimP = n;
		}

//		//! Find largest power of two smaller than numberSimilarPatches
//		int n = 1; while (2*n <= nSimP) n *= 2;
//
//		//! Add more patches if their distance is below a threshold
//		const float tau = params.tau * params.tau * sPx * sPx * sPt;
//		const float threshold = std::max(tau, distance[nSimP - 1].first);
//		while (n < distance.size() && distance[n].first <= threshold) n *= 2;
//		nSimP = n/2;

		//! Store indices of most similar patches
		for (unsigned n = 0; n < nSimP; n++)
			index[n] = distance[n].second;

//		if (nSimP > params.nSimilarPatches)
//			printf("SR2 [%d,%d,%d] ~ nsim = %d ~ nsim ratio = %f\n", px,py,pt, nSimP, (float)nSimP/(float)(sWx*sWy*sWt));
	}
	else // nSimilarPatches == 1
		index[0] = pidx;

	//! Save similar patches into 3D groups
	const unsigned w   = sz.width;
	const unsigned wh  = sz.wh;
	const unsigned whc = sz.whc;
	for (unsigned c = 0, k = 0; c < sz.channels; c++)
	for (unsigned ht = 0; ht < sPt; ht++)
	for (unsigned hy = 0; hy < sPx; hy++)
	for (unsigned hx = 0; hx < sPx; hx++)
	for (unsigned n = 0 ; n < nSimP; n++, k++)
	{
		groupNoisy[k] = imNoisy(c * wh + index[n] + ht * whc + hy * w + hx);
		if (!step1)
			groupBasic[k] = imBasic(c * wh + index[n] + ht * whc + hy * w + hx);
	}

	return nSimP;
}
 #endif // MC_PATCHES
#endif

/**
 * @brief Implementation of computeBayesEstimateStep1 using an 
 * external basis provided by the user in the workspace i_mat.covEigVecs.
 * This version considers different non-linear thresholding operators.
 *
 * See computeBayesEstimateStep1 for information about the arguments.
 **/
float computeBayesEstimateStep1_vbm3d(
	std::vector<float> &io_group
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	std::vector<float> &aggreWeights
){
	//! Parameters initialization
	const float beta_sigma = p_params.beta * p_params.sigma;
	const float beta_sigma2 = beta_sigma * beta_sigma;
	const float betaM_sigma = p_params.betaMean * p_params.sigma;
	const float betaM_sigma2 = betaM_sigma * betaM_sigma;
#ifndef USE_BETA_FOR_VARIANCE
	const float sigma2 = p_params.sigma * p_params.sigma;
#else
	const float sigma2 = beta_sigma2;
#endif
	const unsigned sPC = p_params.sizePatch * p_params.sizePatch
	                   * p_params.sizePatchTime;

	//! Compute variance to determine the aggregation weight
	float non_zero_coeffs = 0;

	for (unsigned c = 0; c < p_imSize.channels; c++)
	{
		std::vector<float> group_c(io_group.begin() + sPC*p_nSimP * c   ,
		                           io_group.begin() + sPC*p_nSimP *(c+1));

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
		              group_c,
		              i_mat.patch_basis,
		              p_nSimP, sPC, sPC,
		              false, false);

		//! Compute Haar transform inplace
		{
			float *z = i_mat.groupTranspose.data();

			const float isqrt2 = 1./sqrt(2.);
			for (int s = 2, h = 1; s <= p_nSimP; s *= 2, h *= 2)
			for (int i = 0       ; i <  sPC    ; i++)
			for (int n = 0       ; n <  p_nSimP; n += s)
			{
				const float z0 = z[n     + i*p_nSimP]*isqrt2;
				const float z1 = z[n + h + i*p_nSimP]*isqrt2;
				z[n     + i*p_nSimP] = (z0 + z1);
				z[n + h + i*p_nSimP] = (z0 - z1);
			}
		}

		//! Thresholding
		{
			//! Impose order invariance
			if (p_params.orderInvariance)
			{
				for (int i = 0; i < sPC; i++)
				{
					float *z = i_mat.groupTranspose.data() + i*p_nSimP;

					float var_i = 0.; // variance for component i
					for (int n = 1; n < p_nSimP; ++n) var_i += z[n] * z[n];
					var_i /= (float)(p_nSimP - 1);

//					printf("i = %d\n",i);
//					for (int n = 0; n < p_nSimP; ++n) printf("z[n] = % 6.2f z[n] * z[n] = % 6.2f\n", z[n], z[n]*z[n]);
//					printf("var_i = % 6.2f b*sigma2 = % 6.2f\n", var_i, beta_sigma2);

					// Wiener filtering instead of thresholding
					float var_0 = std::max(z[0]*z[0] - sigma2, 0.f);
					float w_0 = var_0 / (var_0 + betaM_sigma2);  
					z[0] *= w_0;

					var_i = std::max(var_i - sigma2, 0.f);
					float w_i = var_i/(var_i + beta_sigma2);
					for (int n = 1; n < p_nSimP; ++n) z[n] *= w_i;

					non_zero_coeffs += w_0 * w_0 + (float)(p_nSimP - 1) * w_i * w_i;

//					// threshold dc component
//					if (z[0] * z[0] > betaM_sigma2 || (i == 0)) non_zero_coeffs++;
//					else z[0] = 0.f;
//
//					// threshold remaining components
//					if (var_i > beta_sigma2) non_zero_coeffs += (p_nSimP - 1);
//					else for (int n = 1; n < p_nSimP; ++n) z[n] = 0.f;
				}

//				while(1) int a = 1;
			}
			else
			{
				float *z = i_mat.groupTranspose.data();
				for (unsigned k = 0; k < sPC; ++k)
				for (unsigned i = 0; i < p_nSimP; ++i, ++z)
					if (*z * *z > beta_sigma2 || (k == 0 && i == 0)) non_zero_coeffs += 1;
					else *z = 0.f;
			}
		}

		//! Invert Haar transform inplace
		{
			float *z = i_mat.groupTranspose.data();
			const float isqrt2 = 1./sqrt(2.);
			for (int s = p_nSimP, h = s/2; s >= 2; s /= 2, h /= 2)
			for (int i = 0       ; i <  sPC    ; i++)
			for (int n = 0       ; n <  p_nSimP; n += s)
			{
				const float z0 = z[n     + i*p_nSimP]*isqrt2;
				const float z1 = z[n + h + i*p_nSimP]*isqrt2;
				z[n     + i*p_nSimP] = (z0 + z1);
				z[n + h + i*p_nSimP] = (z0 - z1);
			}
		}

		//! hX' = Z'*U'
		productMatrix(group_c,
		              i_mat.groupTranspose,
		              i_mat.patch_basis_inv,
		              p_nSimP, sPC, sPC,
		              false, true);

		//! Copy channel back into vector
		std::copy(group_c.begin(), group_c.end(),
		          io_group.begin() + sPC*p_nSimP*c);
	}

	// aggregation weights
	const float variance = 1/sigma2/(float)non_zero_coeffs;
	for (unsigned n = 0; n < p_nSimP; n++)
		aggreWeights[n] = variance;

	// return percentage of captured variance
	return variance;
}

/**
 * @brief Compute the Bayes estimation assuming a low rank covariance matrix.
 *
 * @param io_group: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *    - groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 *    - baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 *    - covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *    - covMatTmp: allocated memory. Used to process the Bayes estimate;
 *    - tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_params: see processStep1 for more explanation.
 * @param p_nSimP: number of similar patches.
 * @param aggreWeights: output aggregation weights.
 *
 * @return none.
 **/
float computeBayesEstimateStep1(
	std::vector<float> &io_group
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

	if (p_nSimP > 1)
	{
		for (unsigned j = 0; j < sPC; j++)
		{
			float v = io_group[j * p_nSimP];
			for (unsigned i = 1; i < p_nSimP; i++)
				if (v != io_group[j * p_nSimP + i])
					goto not_equal;
		}

		//! All patches are equal ~ do nothing
		return 0.f;
	}

not_equal:
	//! Not all patches are equal ~ denoise
	return
		computeBayesEstimateStep1_vbm3d(io_group, i_mat,
			io_nInverseFailed, p_size, p_params, p_nSimP, aggreWeights);
}

/**
 * @brief Implementation of computeBayesEstimateStep2 computing the
 * principal directions of the a priori covariance matrix. This functions
 * computes the eigenvectors/values of the data covariance matrix using LAPACK.
 *
 * See computeBayesEstimateStep2 for information about the arguments.
 **/
float computeBayesEstimateStep2_vbm3d(
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
	const float betaM_sigma = p_params.betaMean * p_params.sigma;
	const float betaM_sigma2 = betaM_sigma * betaM_sigma;
#ifndef USE_BETA_FOR_VARIANCE
	const float sigma2 = p_params.sigma * p_params.sigma;
#else
	const float sigma2 = beta_sigma2;
#endif
	const unsigned sPC  = p_params.sizePatch * p_params.sizePatch
	                    * p_params.sizePatchTime;

	//! Compute variance to determine the aggregation weight
	float variance = 0.;

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

		//! Project noisy patches over basis: Z' = X'*U (to compute variances)
		productMatrix(i_mat.groupTransposeNoisy,
		              groupNoisy_c,
		              i_mat.patch_basis,
		              p_nSimP, sPC, sPC,
		              false, false);

		//! Project basic patches over basis: Z' = X'*U (to compute variances)
		productMatrix(i_mat.groupTranspose,
		              groupBasic_c,
		              i_mat.patch_basis,
		              p_nSimP, sPC, sPC,
		              false, false);

		//! Compute Haar transform inplace
		{
			float *y = i_mat.groupTranspose.data();
			float *z = i_mat.groupTransposeNoisy.data();

			const float isqrt2 = 1./sqrt(2.);
			for (int s = 2, h = 1; s <= p_nSimP; s *= 2, h *= 2)
 			for (int i = 0       ; i <  sPC    ; i++)
			for (int n = 0       ; n <  p_nSimP; n += s)
			{
				const float z0 = z[n     + i*p_nSimP]*isqrt2;
				const float z1 = z[n + h + i*p_nSimP]*isqrt2;
				z[n     + i*p_nSimP] = (z0 + z1);
				z[n + h + i*p_nSimP] = (z0 - z1);

				const float y0 = y[n     + i*p_nSimP]*isqrt2;
				const float y1 = y[n + h + i*p_nSimP]*isqrt2;
				y[n     + i*p_nSimP] = (y0 + y1);
				y[n + h + i*p_nSimP] = (y0 - y1);
			}
		}

		//! Wiener filter
		{
			if (p_params.orderInvariance)
			{
				for (int i = 0; i < sPC; i++)
				{
					float *y = i_mat.groupTranspose.data()      + i*p_nSimP;
					float *z = i_mat.groupTransposeNoisy.data() + i*p_nSimP;

					// dc component
					float var_0 = y[0]*y[0];
					float w_0 = var_0 / (var_0 + betaM_sigma2);  
					z[0] *= w_0;

					// rest of components
					float var_i = 0.; // variance for component i
					for (int n = 1; n < p_nSimP; ++n) var_i += y[n] * y[n];
					var_i /= (float)(p_nSimP - 1);
					float w_i = var_i/(var_i + beta_sigma2);
					for (int n = 1; n < p_nSimP; ++n) z[n] *= w_i;

					variance += w_0 * w_0 + (float)(p_nSimP - 1) * w_i * w_i;
				}
			}
			else
			{
				const float *y = i_mat.groupTranspose.data();
				float *z = i_mat.groupTransposeNoisy.data();

				for (int i = 0; i < sPC; ++i)
				for (int n = 0; n < p_nSimP; ++n, ++y, ++z)
				{
					const float w = *y**y / ( *y**y + beta_sigma2);
					*z *= w;
					variance += w*w;
				}
			}
		}

		//! Invert Haar transform inplace
		{
			float *z = i_mat.groupTransposeNoisy.data();

			const float isqrt2 = 1./sqrt(2.);
			for (int s = p_nSimP, h = s/2; s >= 2; s /= 2, h /= 2)
			for (int i = 0       ; i <  sPC    ; i++)
			for (int n = 0       ; n <  p_nSimP; n += s)
			{
				const float z0 = z[n     + i*p_nSimP]*isqrt2;
				const float z1 = z[n + h + i*p_nSimP]*isqrt2;
				z[n     + i*p_nSimP] = (z0 + z1);
				z[n + h + i*p_nSimP] = (z0 - z1);
			}
		}

		//! invert DCT transform
		productMatrix(groupNoisy_c,
		              i_mat.groupTransposeNoisy,
		              i_mat.patch_basis_inv,
		              p_nSimP, sPC, sPC,
		              false, true);

		//! Copy channel back into vector
		std::copy(groupNoisy_c.begin(), groupNoisy_c.end(),
		          io_groupNoisy.begin() + sPC*p_nSimP*c);
	}

	// aggregation weights
	variance = 1./sigma2/variance;
	for (unsigned n = 0; n < p_nSimP; n++)
		aggreWeights[n] = variance;


	// return percentage of captured variance
	return variance;
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

	if (p_nSimP > 1)
	{
		for (unsigned j = 0; j < sPC; j++)
		{
			float v = io_groupNoisy[j * p_nSimP];
			for (unsigned i = 1; i < p_nSimP; i++)
				if (v != io_groupNoisy[j * p_nSimP + i])
					goto not_equal;
		}

		//! All patches are equal ~ do nothing
		return 0.f;
	}

not_equal:
	//! Not all patches are equal ~ denoise
	return
		computeBayesEstimateStep2_vbm3d(io_groupNoisy, i_groupBasic, i_mat,
			io_nInverseFailed, p_size, p_params, p_nSimP, aggreWeights);
}

#ifdef MC_PATCHES
/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param im: update the image with estimate values;
 * @param weight: update corresponding weight, used later in the weighted aggregation;
 * @param mask: update values of mask: set to true the index of an used patch;
 * @param group: contains estimated values of all similar patches in the 3D group;
 * @param aggreWeights: input aggregation weights.
 * @param index: contains index of all similar patches contained in group;
 * @param imSize: size of im;
 * @param params: see processStep1 for more explanation.
 * @param nSimP: number of similar patches.
 *
 * @return masked: number of processable pixels that were flaged non-processable.
 **/
int computeAggregationStep1(
	Video<float> &im
,	Video<float> &weight
,	Video<char>  &mask
,	std::vector<float> const& group
,	std::vector<float> const& aggreWeights
,	std::vector<float> const& aggreWindow
,	std::vector<unsigned> const& index
,	const nlbParams &params
,	const unsigned nSimP
){
	//! Parameters initializations
	const unsigned chnls = im.sz.channels;
	const unsigned sPx   = params.sizePatch;
	const unsigned sPt   = params.sizePatchTime;

	const unsigned w   = im.sz.width;
	const unsigned h   = im.sz.height;
	const unsigned wh  = im.sz.wh;
	const unsigned whc = im.sz.whc;
	const unsigned sPC = sPx*sPx*sPt;

	//! Aggregate estimates
	int masked = 0;
	for (unsigned n = 0; n < nSimP; n++)
	{
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			unsigned offset = n*sPt*3;
			for (unsigned ht = 0, i = 0; ht < sPt; ht++)
			{
				unsigned px = index[offset++];
				unsigned py = index[offset++];
				unsigned pt = index[offset++];
				for (unsigned hy = 0; hy < sPx; hy++)
				for (unsigned hx = 0; hx < sPx; hx++, i++, k++)
				{
					im(px + hx, py + hy, pt, c) +=
						aggreWeights[n] * aggreWindow[i] * group[k * nSimP + n];

					if (c == 0) weight(px + hx, py + hy, pt) += aggreWeights[n] * aggreWindow[i];
				}
			}
		}

		//! Use Paste Trick
		if (params.doPasteBoost)
		{
			unsigned offset = n*sPt*3;
			unsigned px = index[n*sPt*3 + 0];
			unsigned py = index[n*sPt*3 + 1];
			unsigned pt = index[n*sPt*3 + 2];

			if (mask(px,py,pt)) masked++;
			mask(px,py,pt) = false;

			if ((px >     2*sPx) && mask(px - 1, py, pt)) masked++;
			if ((px < w - 2*sPx) && mask(px + 1, py, pt)) masked++;
			if ((py >     2*sPx) && mask(px, py - 1, pt)) masked++;
			if ((py < h - 2*sPx) && mask(px, py + 1, pt)) masked++;

			if (px >     2*sPx) mask(px - 1, py, pt) = false;
			if (px < w - 2*sPx) mask(px + 1, py, pt) = false;
			if (py >     2*sPx) mask(px, py - 1, pt) = false;
			if (py < h - 2*sPx) mask(px, py + 1, pt) = false;
		}
	}

	if (!params.doPasteBoost)
	{
		masked++;
		mask(index[0], index[1], index[2]) = false;
	}

	return masked;
}
#else //MC_PATCHES
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
,	std::vector<float> const& i_group
,	std::vector<float> const& aggreWeights
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

	//! Aggregate estimates
	int masked = 0;
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			const unsigned ij = ind  + c * wh;
			for (unsigned pt = 0, i = 0; pt < sPt; pt++)
			for (unsigned py = 0       ; py < sPx; py++)
			for (unsigned px = 0       ; px < sPx; px++, k++, i++)
				io_im(ij + pt * whc + py * w + px) +=
					aggreWeights[n] * aggreWindow[i] * i_group[k * p_nSimP + n];
		}

		const unsigned ind1 = (ind / whc) * wh + ind % wh;
		for (unsigned pt = 0, i = 0; pt < sPt; pt++)
		for (unsigned py = 0       ; py < sPx; py++)
		for (unsigned px = 0       ; px < sPx; px++, i++)
			io_weight(ind1 + pt * wh + py * w + px) += aggreWeights[n] * aggreWindow[i];

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

	if (!p_params.doPasteBoost)
	{
		masked++;
		io_mask(i_index[0]) = false;
	}

	return masked;
}
#endif //MC_PATCHES

#ifdef MC_PATCHES
/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param im: update the image with estimate values;
 * @param weight: update corresponding weight, used later in the weighted aggregation;
 * @param mask: update values of mask: set to true the index of an used patch;
 * @param group: contains estimated values of all similar patches in the 3D group;
 * @param aggreWeights: input aggregation weights.
 * @param index: contains index of all similar patches contained in group;
 * @param imSize: size of im;
 * @param params: see processStep2 for more explanation;
 * @param nSimP: number of similar patches.
 *
 * @return masked: number of processable pixels that were flaged non-processable.
 **/
int computeAggregationStep2(
	Video<float> &im
,	Video<float> &weight
,	Video<char>  &mask
,	std::vector<float> const& group
,	std::vector<float> const& aggreWeights
,	std::vector<float> const& aggreWindow
,	Video<float> &variance
,	std::vector<unsigned> const& index
,	const nlbParams &params
,	const unsigned nSimP
){
	//! Parameters initializations
	const unsigned chnls = im.sz.channels;
	const unsigned sPx   = params.sizePatch;
	const unsigned sPt   = params.sizePatchTime;

	const unsigned w   = im.sz.width;
	const unsigned h   = im.sz.height;
	const unsigned wh  = im.sz.wh;
	const unsigned whc = im.sz.whc;
	const unsigned sPC = sPx*sPx*sPt;

	int masked = 0;

	//! Aggregate estimates
	for (unsigned n = 0; n < nSimP; n++)
	{
		for (unsigned c = 0, k = 0; c < chnls; c++)
		for (unsigned ht = 0, i = 0; ht < sPt; ht++)
		{
			unsigned px = index[sPt*3*n + 3*ht + 0];
			unsigned py = index[sPt*3*n + 3*ht + 1];
			unsigned pt = index[sPt*3*n + 3*ht + 2];
			for (unsigned hy = 0; hy < sPx; hy++)
			for (unsigned hx = 0; hx < sPx; hx++, i++, k++)
			{
				im(px + hx, py + hy, pt, c) +=
					aggreWeights[n] * aggreWindow[i] * group[k * nSimP + n];

				if (c == 0) weight(px + hx, py + hy, pt) += aggreWeights[n] * aggreWindow[i];
			}

			//! Use Paste Trick
			if (c == 0 && ht == 0 && params.doPasteBoost)
			{
				if (mask(px,py,pt)) masked++;
				mask(px,py,pt) = false;

				if ((px >     2*sPx) && mask(px - 1, py, pt)) masked++;
				if ((px < w - 2*sPx) && mask(px + 1, py, pt)) masked++;
				if ((py >     2*sPx) && mask(px, py - 1, pt)) masked++;
				if ((py < h - 2*sPx) && mask(px, py + 1, pt)) masked++;

				if (px >     2*sPx) mask(px - 1, py, pt) = false;
				if (px < w - 2*sPx) mask(px + 1, py, pt) = false;
				if (py >     2*sPx) mask(px, py - 1, pt) = false;
				if (py < h - 2*sPx) mask(px, py + 1, pt) = false;
			}
		}
	}

	if (!params.doPasteBoost)
	{
		masked++;
		mask(index[0], index[1], index[2]) = false;
	}

	return masked;
}
#else //MC_PATCHES
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

	//! Aggregate estimates
	int masked = 0;
	for (unsigned n = 0; n < p_nSimP; n++)
	{
		const unsigned ind = i_index[n];
		for (unsigned c = 0, k = 0; c < chnls; c++)
		{
			const unsigned ij = ind + c * wh;
			for (unsigned pt = 0, i = 0; pt < sPt; pt++)
			for (unsigned py = 0       ; py < sPx; py++)
			for (unsigned px = 0       ; px < sPx; px++, k++, i++)
				io_im(ij + pt * whc + py * w + px) +=
					aggreWeights[n] * aggreWindow[i] * i_group[k * p_nSimP + n];
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
		}
	}

	if (!p_params.doPasteBoost)
	{
		masked++;
		io_mask(i_index[0]) = false;
	}

	return masked;
}
#endif //MC_PATCHES

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

} // namespace
