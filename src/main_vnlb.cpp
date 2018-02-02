/*
 * Original work: Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Modified work: Copyright (c) 2014, Pablo Arias <pariasm@gmail.com>
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

#include "utils/utilities.h"
#include "core/nldct.h"
#include "utils/cmd_option.h"

using namespace std;

/**
 * @file   main.cpp
 * @brief  Main executable file
 *
 *
 *
 * @author MARC LEBRUN  <marc.lebrun.ik@gmail.com>
 * @author PABLO ARIAS  <pariasm@gmail.com>
 **/

enum Mode { BSIC_DENO, BSIC_ONLY, DENO_ONLY, NISY_ONLY };

int main(int argc, char **argv)
{
	clo_usage("Video NL-Bayes video denoising");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	//! Paths to input/output sequences
	using std::string;
	const string  input_path = clo_option("-i"    , ""              , "< input sequence");
	const string  inbsc_path = clo_option("-b"    , ""              , "< input basic sequence");
	const string  noisy_path = clo_option("-nisy" , ""              , "> noisy sequence");
	const string  final_path = clo_option("-deno" , "deno_%03d.png" , "> denoised sequence");
	const string  basic_path = clo_option("-bsic" , "bsic_%03d.png" , "> basic denoised sequence");
	const string   diff_path = clo_option("-diff" , "diff_%03d.png" , "> difference sequence");

	const unsigned firstFrame = clo_option("-f", 0, "first frame");
	const unsigned lastFrame  = clo_option("-l", 0, "last frame");
	const unsigned frameStep  = clo_option("-s", 1, "frame step");

	//! Paths to optical flow
	const string  fflow_path = clo_option("-fof", "", "< input forward  optical flow");
	const string  bflow_path = clo_option("-bof", "", "< input backward optical flow");

	//! General parameters
	const float sigma = clo_option("-sigma", 0.f, "Add noise of standard deviation sigma");
	const bool has_noise = (bool) clo_option("-has-noise"   , false, "> input image already has noise");
	const bool do_bias   = (bool) clo_option("-compute-bias", false, "> compute bias outputs");
	const bool verbose   = (bool) clo_option("-verbose"     , true , "> verbose output");
	const unsigned print_prms = (unsigned) clo_option("-print-prms", 0, "> prints parameters for given channels");

	//! Video NLB parameters
	const int time_search1  = clo_option("-wt1", 0  , "> Search window temporal radius, step 1");
	const int time_search2  = clo_option("-wt2", 0  , "> Search window temporal radius, step 2");
	const int space_search1 = clo_option("-wx1",-1  , "> Search window spatial radius, step 1");
	const int space_search2 = clo_option("-wx2",-1  , "> Search window spatial radius, step 2");
	const int patch_sizex1  = clo_option("-px1",-1  , "> Spatial patch size, step 1");
	const int patch_sizex2  = clo_option("-px2",-1  , "> Spatial patch size, step 2");
	const int patch_sizet1  = clo_option("-pt1", 1  , "> Temporal patch size, step 1");
	const int patch_sizet2  = clo_option("-pt2", 1  , "> Temporal patch size, step 2");
	const int num_patches1  = clo_option("-np1",-1  , "> Number of similar patches, step 1");
	const int num_patches2  = clo_option("-np2",-1  , "> Number of similar patches, step 2");
	const float tau1        = clo_option("-t1" ,-1.f, "> Step 1 distance threshold");
	const float tau2        = clo_option("-t2" ,-1.f, "> Step 2 distance threshold");
	const float beta1       = clo_option("-b1" ,-1.f, "> Noise correction factor beta, step 1");
	const float beta2       = clo_option("-b2" ,-1.f, "> Noise correction factor beta, step 2");
#ifdef VBM3D_SEARCH
	const int space_search_f1= clo_option("-wxf1",-1  , "> Search window for predictive search, step 1");
	const int space_search_f2= clo_option("-wxf2",-1  , "> Search window for predictive search, step 2");
	const int num_patches_f1 = clo_option("-npf1",-1  , "> Number of similar patches per frame, step 1");
	const int num_patches_f2 = clo_option("-npf2",-1  , "> Number of similar patches per frame, step 2");
	const float dsub1        = clo_option("-dsub1" ,-1.f, "> Distance bias, step 1");
	const float dsub2        = clo_option("-dsub2" ,-1.f, "> Distance bias, step 2");
#endif
	const float beta_mean1  = clo_option("-bm1",-1.f, "> Noise correction factor beta for mean, step 1");
	const float beta_mean2  = clo_option("-bm2",-1.f, "> Noise correction factor beta for mean, step 2");
	const bool no_paste1  = (bool) clo_option("-no-paste1", false , "> disable paste trick, step 1");
	const bool no_paste2  = (bool) clo_option("-no-paste2", false , "> disable paste trick, step 2");
	const bool no_step1   = (bool) clo_option("-no-step1" , false , "> disable patch skipping, step 1");
	const bool no_step2   = (bool) clo_option("-no-step2" , false , "> disable patch skipping, step 2");
	const bool agg_win1   = (bool) clo_option("-agg-win1" , false , "> aggregation window, step 1");
	const bool agg_win2   = (bool) clo_option("-agg-win2" , false , "> aggregation window, step 2");


	//! Check inputs
	if (input_path == "")
	{
		fprintf(stderr, "%s: no input sequence.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if (patch_sizex1 == 0 && patch_sizex2 > 0 && inbsc_path == "")
	{
		fprintf(stderr, "%s: if px1 = 0 and px2 > 0, a basic sequence path must be given.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if ((patch_sizex1 < 0 && patch_sizex1 != -1) ||
	    (patch_sizex2 < 0 && patch_sizex2 != -1) ||
	    (patch_sizet1 < 0 || patch_sizet2 <   0) )
	{
		fprintf(stderr, "%s: px1, px2, pt1 and pt2 cannot be negative.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if ((num_patches1 < 0 && num_patches1 != -1) ||
	    (num_patches2 < 0 && num_patches2 != -1) )
	{
		fprintf(stderr, "%s: np1 and np2 cannot be negative.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if ((space_search1 < 0 && space_search1 != -1) ||
	    (space_search2 < 0 && space_search2 != -1) ||
	    ( time_search1 < 0 ||  time_search2 <   0) )
	{
		fprintf(stderr, "%s: wx1, wx2, wt1 and wt2 cannot be negative.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	if (patch_sizex1 > 0 && inbsc_path != "")
		fprintf(stderr, "\x1b[33;1mWarning:\x1b[0m a basic sequence path ignored since px1 > 0.\n");

	if ((fflow_path != "" && bflow_path == "") || 
	    (fflow_path == "" && bflow_path != ""))
	{
		fprintf(stderr, "Only one oflow path provided.\nTry `%s --help' for more information.\n",
				argv[0]);
		return EXIT_FAILURE;
	}


	//! Determine mode
	Mode mode;
	if ((patch_sizex1 != 0) && (patch_sizex2 != 0)) mode = BSIC_DENO;
	if ((patch_sizex1 != 0) && (patch_sizex2 == 0)) mode = BSIC_ONLY;
	if ((patch_sizex1 == 0) && (patch_sizex2 != 0)) mode = DENO_ONLY;
	if ((patch_sizex1 == 0) && (patch_sizex2 == 0)) mode = NISY_ONLY;

	bool use_oflow = (fflow_path != "");


	//! Only print parameters
	if (print_prms)
	{
		VideoSize tmp;
		tmp.channels = print_prms;

		//! Compute denoising default parameters
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, 1, sigma, tmp, verbose, time_search1, time_search1, patch_sizet1);
		VideoNLB::initializeNlbParameters(prms2, 2, sigma, tmp, verbose, time_search2, time_search2, patch_sizet2);

		//! Override with command line parameters
		if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
		if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
		if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, tmp, (unsigned)patch_sizex1);;
		if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, tmp, (unsigned)patch_sizex2);;
		if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
		if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
		if (beta1         >= 0) prms1.beta = beta1;
		if (beta2         >= 0) prms2.beta = beta2;

#ifdef VBM3D_SEARCH
		if (space_search_f1 >= 0) prms1.sizeSearchWindowPred = space_search_f1;
		if (space_search_f2 >= 0) prms2.sizeSearchWindowPred = space_search_f2;
		if (num_patches_f1  >= 0) prms1.nSimilarPatchesPred = num_patches_f1;
		if (num_patches_f2  >= 0) prms2.nSimilarPatchesPred = num_patches_f2;
		if (dsub1           >= 0) prms1.dsub = dsub1;
		if (dsub2           >= 0) prms2.dsub = dsub2;
		if (tau1            >= 0) prms1.tau = tau1;
		if (tau2            >= 0) prms2.tau = tau2;
#endif
		if (agg_win1) prms1.agg_window = true;
		if (agg_win2) prms2.agg_window = true;

		if (no_paste1) prms1.doPasteBoost = false;
		if (no_paste2) prms2.doPasteBoost = false;
		if (no_step1)  prms1.offSet = 1;
		if (no_step2)  prms2.offSet = 1;
 
		VideoNLB::printNlbParameters(prms1);
		VideoNLB::printNlbParameters(prms2);

		return EXIT_SUCCESS;
	}


	//! Declarations
	Video<float> original, noisy, basic, final, diff;
	Video<float> fflow, bflow;

	//! Load input videos
	                       original.loadVideo(input_path, firstFrame, lastFrame, frameStep);
	if (mode == DENO_ONLY) basic   .loadVideo(inbsc_path, firstFrame, lastFrame, frameStep);
	if (use_oflow)         fflow   .loadVideo(fflow_path, firstFrame, lastFrame, frameStep);
	if (use_oflow)         bflow   .loadVideo(bflow_path, firstFrame, lastFrame, frameStep);

	//! Add noise
	if (has_noise)
	{
		if (verbose) printf("Input video has noise of sigma = %f\n", sigma);
		noisy = original;
	}
	else if (sigma)
		VideoUtils::addNoise(original, noisy, sigma, verbose);

	//! Save noisy video
	if (noisy_path != "")
	{
		if (verbose) printf("Saving noisy video\n");
		noisy.saveVideo(noisy_path, firstFrame, frameStep);
	}

	if (mode == NISY_ONLY)
	{
		if (!has_noise && sigma)
		{
			float noisy_psnr = -1, noisy_rmse = -1;
			VideoUtils::computePSNR(original, noisy, noisy_psnr, noisy_rmse);
			writingMeasures("measures.txt", sigma, noisy_psnr, noisy_rmse, 0, true, "_noisy");
			return EXIT_SUCCESS;
		}
		else
		{
			fprintf(stderr, "Noisy video not saved since is equal to original input.\n");
			return EXIT_FAILURE;
		}
	}

	//! Denoising
	if (verbose) printf("Running Video NL-Bayes on the noisy video\n");

	//! Compute denoising default parameters
	VideoNLB::nlbParams prms1, prms2;
	VideoNLB::initializeNlbParameters(prms1, 1, sigma, noisy.sz, verbose, time_search1, time_search1, patch_sizet1);
	VideoNLB::initializeNlbParameters(prms2, 2, sigma, noisy.sz, verbose, time_search2, time_search2, patch_sizet2);

	//! Override with command line parameters
	if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
	if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
	if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, noisy.sz, (unsigned)patch_sizex1);;
	if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, noisy.sz, (unsigned)patch_sizex2);;
	if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
	if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
	if (beta1         >= 0) prms1.beta = beta1;
	if (beta2         >= 0) prms2.beta = beta2;

	if (no_paste1) prms1.doPasteBoost = false;
	if (no_paste2) prms2.doPasteBoost = false;
	if (no_step1)  prms1.offSet = 1;
	if (no_step2)  prms2.offSet = 1;

	//! Percentage or processed groups of patches over total number of pixels
	std::vector<float> groupsRatio;

	//! Run denoising algorithm
#ifndef DEBUG_COMPUTE_GROUP_ERROR
	if (use_oflow)
		groupsRatio = VideoNLB::runNlBayes(noisy, fflow, bflow, basic, final, prms1, prms2);
	else
		groupsRatio = VideoNLB::runNlBayes(noisy, basic, final, prms1, prms2);
#else
	if (use_oflow)
		groupsRatio = VideoNLB::runNlBayes(noisy, fflow, bflow, basic, final, prms1, prms2, original);
	else
		groupsRatio = VideoNLB::runNlBayes(noisy, basic, final, prms1, prms2, original);
#endif

	//! Compute PSNR and RMSE
	float final_psnr = -1, final_rmse = -1, basic_psnr = -1, basic_rmse = -1;
	if (prms2.sizePatch) VideoUtils::computePSNR(original, final, final_psnr, final_rmse);
	                     VideoUtils::computePSNR(original, basic, basic_psnr, basic_rmse);

	if (verbose)
	{
	    printf("basic PSNR =\t%f\tRMSE =\t%f\n", basic_psnr, basic_rmse);
	    printf("final PSNR =\t%f\tRMSE =\t%f\n", final_psnr, final_rmse);
	}
	else if (mode != NISY_ONLY)
		printf("%f\n", mode == BSIC_ONLY ? basic_psnr : final_psnr);


	//! Write measures
	writingMeasures("measures.txt", sigma, basic_psnr, basic_rmse, groupsRatio[0], true  , "_basic");
	writingMeasures("measures.txt", sigma, final_psnr, final_rmse, groupsRatio[1], false , "_final");

	//! Compute Difference
	if (prms2.sizePatch) VideoUtils::computeDiff(original, final, diff, sigma);

	//! Save output sequences
	if (verbose) printf("\nSaving output sequences\n");

	if (prms2.sizePatch) final.saveVideo(final_path, firstFrame, frameStep);
	if (prms2.sizePatch) diff .saveVideo( diff_path, firstFrame, frameStep);
	                     basic.saveVideo(basic_path, firstFrame, frameStep);

	if (verbose) printf("Done\n");
	return EXIT_SUCCESS;
}
