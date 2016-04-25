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

#include "Utilities/Utilities.h"
#include "NlBayes/VideoNLBayes.h"
#include "Utilities/cmd_option.h"

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
	const string  noisy_path = clo_option("-nisy" , "nisy_%03d.png" , "> noisy sequence");
	const string  final_path = clo_option("-deno" , "deno_%03d.png" , "> denoised sequence");
	const string  basic_path = clo_option("-bsic" , "bsic_%03d.png" , "> basic denoised sequence");
	const string   diff_path = clo_option("-diff" , "diff_%03d.png" , "> difference sequence");
	// TODO: these should be determined automatically from the other outputs.
	const string   bias_path = clo_option("-bdeno", "bdeno_%03d.png", "> bias sequence");
	const string bbasic_path = clo_option("-bbsic", "bbsic_%03d.png", "> bias basic sequence");
	const string  bdiff_path = clo_option("-bdiff", "bdiff_%03d.png", "> bias difference sequence");

	const unsigned firstFrame = clo_option("-f", 0, "first frame");
	const unsigned lastFrame  = clo_option("-l", 0, "last frame");
	const unsigned frameStep  = clo_option("-s", 1, "frame step");

	//! General parameters
	const float sigma = clo_option("-sigma", 0.f, "Add noise of standard deviation sigma");
	const bool has_noise = (bool) clo_option("-has-noise"   , false, "> input image already has noise");
	const bool do_bias   = (bool) clo_option("-compute-bias", false, "> compute bias outputs");
	const bool verbose   = (bool) clo_option("-verbose"     , true , "> verbose output");
	const unsigned print_prms = (unsigned) clo_option("-print-prms", 0, "> prints parameters for given channels");

	//! Video NLB parameters
	const bool flat_area1 = (bool) clo_option("-flat-area1", false , "> use flat area trick, step 1");
	const bool flat_area2 = (bool) clo_option("-flat-area2", false , "> use flat area trick, step 2");
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
	const int rank1         = clo_option("-r1" , 4  , "> Rank or covariance matrix, step 1");
	const int rank2         = clo_option("-r2" , 4  , "> Rank or covariance matrix, step 2");
	const float beta1       = clo_option("-b1" ,-1.f, "> Noise correction factor beta, step 1");
	const float beta2       = clo_option("-b2" ,-1.f, "> Noise correction factor beta, step 2");
	const float beta_mean1  = clo_option("-bm1",-1.f, "> Noise correction factor beta for mean, step 1");
	const float beta_mean2  = clo_option("-bm2",-1.f, "> Noise correction factor beta for mean, step 2");
	const float tau2        = clo_option("-t2" ,-1.f, "> Step 2 distance threshold");

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
	    (num_patches2 < 0 && num_patches2 != -1) ||
	    (rank1 < 0 || rank2 <   0) )
	{
		fprintf(stderr, "%s: np1, np2, r1 and r2 cannot be negative.\nTry `%s --help' for more information.\n",
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


	//! Determine mode
	Mode mode;
	if ((patch_sizex1 != 0) && (patch_sizex2 != 0)) mode = BSIC_DENO;
	if ((patch_sizex1 != 0) && (patch_sizex2 == 0)) mode = BSIC_ONLY;
	if ((patch_sizex1 == 0) && (patch_sizex2 != 0)) mode = DENO_ONLY;
	if ((patch_sizex1 == 0) && (patch_sizex2 == 0)) mode = NISY_ONLY;


	//! Only print parameters
	if (print_prms)
	{
		VideoSize tmp;
		tmp.channels = print_prms;

		//! Compute denoising default parameters
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, 1, sigma, tmp, flat_area1, verbose, time_search1, time_search1, patch_sizet1);
		VideoNLB::initializeNlbParameters(prms2, 2, sigma, tmp, flat_area2, verbose, time_search2, time_search2, patch_sizet2);

		//! Override with command line parameters
		if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
		if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
		if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, tmp, (unsigned)patch_sizex1);;
		if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, tmp, (unsigned)patch_sizex2);;
		if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
		if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
		if (beta1         >= 0) prms1.beta = beta1;
		if (beta2         >= 0) prms2.beta = beta2;
		if (beta_mean1    >= 0) prms1.betaMean = beta_mean1;
		if (beta_mean2    >= 0) prms2.betaMean = beta_mean2;
		if (tau2          >= 0) VideoNLB::setTau(prms2, tmp, tau2);

		prms1.rank = rank1;
		prms2.rank = rank2;
 
		VideoNLB::printNlbParameters(prms1);
		VideoNLB::printNlbParameters(prms2);

		return EXIT_SUCCESS;
	}


	//! Declarations
	Video<float> original, noisy, basic, final, diff;

	//! Load input videos
	                       original.loadVideo(input_path, firstFrame, lastFrame, frameStep);
	if (mode == DENO_ONLY) basic   .loadVideo(inbsc_path, firstFrame, lastFrame, frameStep);

	//! Add noise
	if (has_noise)
	{
		if (verbose) printf("Input video has noise of sigma = %f\n", sigma);
		noisy = original;
	}
	else if (sigma)
		VideoUtils::addNoise(original, noisy, sigma, verbose);

	//! Save noisy video
	if (mode == NISY_ONLY)
	{
		if (!has_noise && sigma)
		{
			float noisy_psnr = -1, noisy_rmse = -1;
			VideoUtils::computePSNR(original, noisy, noisy_psnr, noisy_rmse);
			writingMeasures("measures.txt", sigma, noisy_psnr, noisy_rmse, 0, true, "_noisy");

			if (verbose) printf("Saving noisy video\n");
			noisy.saveVideo(noisy_path, firstFrame, frameStep);

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
	VideoNLB::initializeNlbParameters(prms1, 1, sigma, noisy.sz, flat_area1, verbose, time_search1, time_search1, patch_sizet1);
	VideoNLB::initializeNlbParameters(prms2, 2, sigma, noisy.sz, flat_area2, verbose, time_search2, time_search2, patch_sizet2);

	//! Override with command line parameters
	if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
	if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
	if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, noisy.sz, (unsigned)patch_sizex1);;
	if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, noisy.sz, (unsigned)patch_sizex2);;
	if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
	if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
	if (beta1         >= 0) prms1.beta = beta1;
	if (beta2         >= 0) prms2.beta = beta2;
	if (beta_mean1    >= 0) prms1.betaMean = beta_mean1;
	if (beta_mean2    >= 0) prms2.betaMean = beta_mean2;
	if (tau2          >= 0) VideoNLB::setTau(prms2, noisy.sz, tau2);

	prms1.rank = rank1;
	prms2.rank = rank2;

	//! Percentage or processed groups of patches over total number of pixels
	std::vector<float> groupsRatio;

	//! Run denoising algorithm
	groupsRatio = VideoNLB::runNlBayes(noisy, basic, final, prms1, prms2);

	//! Compute PSNR and RMSE
	float final_psnr = -1, final_rmse = -1, basic_psnr = -1, basic_rmse = -1;
	if (prms2.sizePatch) VideoUtils::computePSNR(original, final, final_psnr, final_rmse);
	                     VideoUtils::computePSNR(original, basic, basic_psnr, basic_rmse);

	if (verbose)
	{
	    printf("basic PSNR =\t%f\tRMSE =\t%f\n", basic_psnr, basic_rmse);
	    printf("final PSNR =\t%f\tRMSE =\t%f\n", final_psnr, final_rmse);
	}

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

#if 0
	//! Computing bias sequence
	Video<float> bias, bias_basic, bias_diff;
	if (do_bias) {
		if (verbose) cout << "Applying NL-Bayes to the original image :" << endl;

		runNlBayes(im, imBasicBias, imBias, imSize, flat_area1, flat_area2, sigma, verbose);

		if (verbose) cout << endl;

		float psnrBias, psnrBiasBasic, rmseBias, rmseBiasBasic;
		computePsnr(im, imBasicBias, psnrBiasBasic, rmseBiasBasic, "imBiasBasic", verbose);
		computePsnr(im, imBias     , psnrBias     , rmseBias     , "imBiasFinal", verbose);

		//! writing measures
		writingMeasures("measures.txt", sigma, psnrBiasBasic, rmseBiasBasic, false, "_bias_basic");
		writingMeasures("measures.txt", sigma, psnrBias     , rmseBias     , false, "_bias      ");

		
		computeDiff(im, imBias, imDiffBias, sigma, 0.f, 255.f, verbose);

		saveImage(argv[7], imBias     , imSize, 0.f, 255.f);
		saveImage(argv[8], imBasicBias, imSize, 0.f, 255.f);
		saveImage(argv[9], imDiffBias , imSize, 0.f, 255.f);
	}


#endif

	if (verbose) printf("Done\n");
	return EXIT_SUCCESS;
}
