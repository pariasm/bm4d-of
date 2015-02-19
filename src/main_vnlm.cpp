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
#include "NlBayes/VideoNLMeans.h"
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

int main(int argc, char **argv)
{
	clo_usage("Video NL-Means video denoising");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	//! Paths to input/output sequences
	using std::string;
	const string  input_path = clo_option("-i"    , ""                   , "< input sequence");
	const string  noisy_path = clo_option("-nisy" , "nisy_%03d.png" , "> noisy sequence");
	const string  final_path = clo_option("-deno" , "deno_%03d.png" , "> denoised sequence");
	const string  basic_path = clo_option("-bsic" , "bsic_%03d.png" , "> basic denoised sequence");
	const string   diff_path = clo_option("-diff" , "diff_%03d.png" , "> difference sequence");
	// TODO: this should be determined automatically from the other outputs.
	const string   bias_path = clo_option("-bdeno", "bdeno_%03d.png", "> bias sequence");
	const string bbasic_path = clo_option("-bbsic", "bbsic_%03d.png", "> bias basic sequence");
	const string  bdiff_path = clo_option("-bdiff", "bdiff_%03d.png", "> bias difference sequence");

	const unsigned firstFrame = clo_option("-f", 0, "first frame");
	const unsigned lastFrame  = clo_option("-l", 0, "last frame");
	const unsigned frameStep  = clo_option("-s", 1, "frame step");

	//! General parameters
	const float sigma = clo_option("-sigma", 0, "Add noise of standard deviation sigma");
	const bool do_bias  = (bool) clo_option("-compute-bias", false, "> compute bias outputs");
	const bool verbose  = (bool) clo_option("-verbose"     , true , "> verbose output");

	//! Video NLB parameters
	const bool flat_area1 = (bool) clo_option("-flat-area1", true , "> use flat area trick (step 1)");
	const bool flat_area2 = (bool) clo_option("-flat-area2", true , "> use flat area trick (step 2)");
	const unsigned time_search_fwd = clo_option("-Wt-fwd", 2, "> Search window, forward  time range");
	const unsigned time_search_bwd = clo_option("-Wt-bwd", 2, "> Search window, backward time range");

	//! Check inputs
	if (input_path == "")
	{
		fprintf(stderr, "%s: no input sequence.\nTry `%s --help' for more information.\n",
				argv[0], argv[0]);
		return EXIT_FAILURE;
	}

	//! Declarations
	Video<float> original, noisy, basic, final, diff;

	//! Load original video
	original.loadVideo(input_path, firstFrame, lastFrame, frameStep);

	//! Add noise
	if (sigma)
	{
		VideoUtils::addNoise(original, noisy, sigma, verbose);

		//! Save noisy video
		if (verbose) printf("Saving noisy video\n");
		noisy.saveVideo(noisy_path, firstFrame, frameStep);
	}

	//! Denoising
	if (verbose) printf("Running Video NL-Bayes on the noisy video\n");
	VideoNLM::nlbParams prms1, prms2;
	VideoNLM::initializeNlbParameters(prms1, 1, sigma, noisy.sz, flat_area1, verbose, time_search_fwd, time_search_bwd);
	VideoNLM::initializeNlbParameters(prms2, 2, sigma, noisy.sz, flat_area2, verbose, time_search_fwd, time_search_bwd);
	VideoNLM::runNlBayes(noisy, basic, final, prms1, prms2);

	//! Compute PSNR and RMSE
	float final_psnr, final_rmse, basic_psnr, basic_rmse;
	VideoUtils::computePSNR(original, basic, basic_psnr, basic_rmse);
	VideoUtils::computePSNR(original, final, final_psnr, final_rmse);

    if (verbose)
	 {
        printf("basic PSNR =\t%f\tRMSE =\t%f\n", basic_psnr, basic_rmse);
        printf("final PSNR =\t%f\tRMSE =\t%f\n", final_psnr, final_rmse);
    }

	//! Write measures
	writingMeasures("measures.txt", sigma, basic_psnr, basic_rmse, true  , "_basic");
	writingMeasures("measures.txt", sigma, final_psnr, final_rmse, false , "_final");

	//! Compute Difference
	VideoUtils::computeDiff(original, final, diff, sigma);

	//! Save output sequences
	if (verbose) printf("Saving output sequences\n");

	final.saveVideo(final_path, firstFrame, frameStep);
	basic.saveVideo(basic_path, firstFrame, frameStep);
	diff .saveVideo( diff_path, firstFrame, frameStep);

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
