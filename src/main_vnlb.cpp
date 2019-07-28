/* Modified work: Copyright (c) 2019, Pablo Arias <pariasm@gmail.com>
 * Original work: Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 *
 * This program is free software: you can use, modify and/or redistribute it
 * under the terms of the GNU Affero General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version. You should have received a copy of this license
 * along this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
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

enum Mode { BSIC_DENO, BSIC_ONLY, DENO_ONLY };

int main(int argc, char **argv)
{
	clo_usage("BM4D-OF video denoising");
	clo_help(" NOTE: Input (<) and output (>) sequences are specified by their paths in printf format.\n");

	// Paths to input/output sequences
	using std::string;
	const string  input_path = clo_option("-i"    , ""              , "< input sequence");
	const string  inbsc_path = clo_option("-b"    , ""              , "< input basic sequence");
	const string  final_path = clo_option("-deno" , "deno_%03d.png" , "> denoised sequence");
	const string  basic_path = clo_option("-bsic" , "bsic_%03d.png" , "> basic denoised sequence");

	const unsigned first_frame = clo_option("-f", 0, "first frame");
	const unsigned last_frame  = clo_option("-l", 0, "last frame");
	const unsigned frame_step  = clo_option("-s", 1, "frame step");

	// Paths to optical flow
	const string  fflow_path = clo_option("-fof", "", "< input forward  optical flow");
	const string  bflow_path = clo_option("-bof", "", "< input backward optical flow");

	// General parameters
	const float sigma = clo_option("-sigma", 0.f, "< standard deviation of the noise");
	const bool verbose   = (bool) clo_option("-verbose"  , true , "< verbose output");
	const bool order_inv = (bool) clo_option("-order-inv", false, "< invariance to patch order in stack");
	const unsigned print_prms = (unsigned) clo_option("-print-prms", 0, "< prints parameters for given channels");

	// BM4D-OF parameters
	const int time_search1  = clo_option("-wt1", 0  , "< search window temporal radius, step 1");
	const int time_search2  = clo_option("-wt2", 0  , "< search window temporal radius, step 2");
	const int space_search1 = clo_option("-wx1",-1  , "< search window spatial radius, step 1");
	const int space_search2 = clo_option("-wx2",-1  , "< search window spatial radius, step 2");
	const int patch_sizex1  = clo_option("-px1",-1  , "< spatial patch size, step 1");
	const int patch_sizex2  = clo_option("-px2",-1  , "< spatial patch size, step 2");
	const int patch_sizet1  = clo_option("-pt1", 1  , "< temporal patch size, step 1");
	const int patch_sizet2  = clo_option("-pt2", 1  , "< temporal patch size, step 2");
	const int num_patches1  = clo_option("-np1",-1  , "< number of similar patches, step 1");
	const int num_patches2  = clo_option("-np2",-1  , "< number of similar patches, step 2");
	const float tau1        = clo_option("-t1" ,-1.f, "< step 1 distance threshold");
	const float tau2        = clo_option("-t2" ,-1.f, "< step 2 distance threshold");
	const float beta1       = clo_option("-b1" ,-1.f, "< noise correction factor beta, step 1");
	const float beta2       = clo_option("-b2" ,-1.f, "< noise correction factor beta, step 2");
#ifdef VBM3D_SEARCH
	const int space_search_f1 = clo_option("-wxf1" ,-1  , "< search window for predictive search, step 1");
	const int space_search_f2 = clo_option("-wxf2" ,-1  , "< search window for predictive search, step 2");
	const int num_patches_f1  = clo_option("-npf1" ,-1  , "< number of similar patches per frame, step 1");
	const int num_patches_f2  = clo_option("-npf2" ,-1  , "< number of similar patches per frame, step 2");
	const float dsub1         = clo_option("-dsub1",-1.f, "< distance bias, step 1");
	const float dsub2         = clo_option("-dsub2",-1.f, "< distance bias, step 2");
#endif
	const float beta_mean1  = clo_option("-bm1",-1.f, "< noise correction factor beta for mean, step 1");
	const float beta_mean2  = clo_option("-bm2",-1.f, "< noise correction factor beta for mean, step 2");
	const bool no_paste1  = (bool) clo_option("-no-paste1", false , "< disable paste trick, step 1");
	const bool no_paste2  = (bool) clo_option("-no-paste2", false , "< disable paste trick, step 2");
	const bool no_step1   = (bool) clo_option("-no-step1" , false , "< disable patch skipping, step 1");
	const bool no_step2   = (bool) clo_option("-no-step2" , false , "< disable patch skipping, step 2");
	const bool agg_win1   = (bool) clo_option("-agg-win1" , false , "< aggregation window, step 1");
	const bool agg_win2   = (bool) clo_option("-agg-win2" , false , "< aggregation window, step 2");


	// Check inputs
	{
		if (input_path == "")
			return fprintf(stderr, "%s: no input sequence.\nTry `%s --help' for more information.\n",
					argv[0], argv[0]), EXIT_FAILURE;

		if (patch_sizex1 == 0 && patch_sizex2 > 0 && inbsc_path == "")
			return fprintf(stderr, "%s: if px1 = 0 and px2 > 0, a basic sequence path must be given.\n"
					"Try `%s --help' for more information.\n", argv[0], argv[0]), EXIT_FAILURE;

		if ((patch_sizex1 < 0 && patch_sizex1 != -1) ||
		    (patch_sizex2 < 0 && patch_sizex2 != -1) ||
		    (patch_sizet1 < 0 || patch_sizet2 <   0) )
			return fprintf(stderr, "%s: px1, px2, pt1 and pt2 cannot be negative.\n"
					"Try `%s --help' for more information.\n", argv[0], argv[0]), EXIT_FAILURE;

		if ((num_patches1 < 0 && num_patches1 != -1) ||
		    (num_patches2 < 0 && num_patches2 != -1) )
			return fprintf(stderr, "%s: np1 and np2 cannot be negative.\n"
					"Try `%s --help' for more information.\n", argv[0], argv[0]), EXIT_FAILURE;

		if ((space_search1 < 0 && space_search1 != -1) ||
		    (space_search2 < 0 && space_search2 != -1) ||
		    ( time_search1 < 0 ||  time_search2 <   0) )
			return fprintf(stderr, "%s: wx1, wx2, wt1 and wt2 cannot be negative.\n"
					"Try `%s --help' for more information.\n", argv[0], argv[0]), EXIT_FAILURE;

		if (patch_sizex1 > 0 && inbsc_path != "")
			fprintf(stderr, "\x1b[33;1mWarning:\x1b[0m a basic sequence path ignored since px1 > 0.\n");

		if ((fflow_path != "" && bflow_path == "") ||
		    (fflow_path == "" && bflow_path != ""))
			return fprintf(stderr, "Only one oflow path provided.\n"
					"Try `%s --help' for more information.\n", argv[0]), EXIT_FAILURE;

		if ((patch_sizex1 == 0) && (patch_sizex2 == 0))
			return fprintf(stderr, "Given patch sizes for both steps are 0. Nothing to do.\n"),
				EXIT_FAILURE;
	}

	// Determine mode
	Mode mode;
	if ((patch_sizex1 != 0) && (patch_sizex2 != 0)) mode = BSIC_DENO;
	if ((patch_sizex1 != 0) && (patch_sizex2 == 0)) mode = BSIC_ONLY;
	if ((patch_sizex1 == 0) && (patch_sizex2 != 0)) mode = DENO_ONLY;

	bool use_oflow = (fflow_path != "") && (last_frame - first_frame > frame_step) ;

	// Only print parameters
	if (print_prms)
	{
		VideoSize tmp;
		tmp.channels = print_prms;

		// Compute denoising default parameters
		VideoNLB::nlbParams prms1, prms2;
		VideoNLB::initializeNlbParameters(prms1, 1, sigma, tmp, verbose);
		VideoNLB::initializeNlbParameters(prms2, 2, sigma, tmp, verbose);

		// Override with command line parameters
		if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
		if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
		if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, tmp, (unsigned)patch_sizex1);;
		if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, tmp, (unsigned)patch_sizex2);;
		if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
		if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
		if (beta1         >= 0) prms1.beta = beta1;
		if (beta2         >= 0) prms2.beta = beta2;
		if (beta1         >= 0) prms1.betaMean = beta1;
		if (beta2         >= 0) prms2.betaMean = beta2;
		if (beta_mean1    >= 0) prms1.betaMean = beta_mean1;
		if (beta_mean2    >= 0) prms2.betaMean = beta_mean2;

		if (patch_sizet1  >= 0) prms1.sizePatchTime = patch_sizet1;
		if (patch_sizet2  >= 0) prms2.sizePatchTime = patch_sizet2;

#ifdef VBM3D_SEARCH
		if (space_search_f1 >= 0) prms1.sizeSearchWindowPred = space_search_f1;
		if (space_search_f2 >= 0) prms2.sizeSearchWindowPred = space_search_f2;
		if (num_patches_f1  >= 0) prms1.nSimilarPatchesPred = num_patches_f1;
		if (num_patches_f2  >= 0) prms2.nSimilarPatchesPred = num_patches_f2;
		if (dsub1           >= 0) prms1.dsub = dsub1;
		if (dsub2           >= 0) prms2.dsub = dsub2;
#endif
		if (tau1            >= 0) prms1.tau = tau1;
		if (tau2            >= 0) prms2.tau = tau2;
		if (agg_win1) prms1.agg_window = true;
		if (agg_win2) prms2.agg_window = true;

		if (no_paste1) prms1.doPasteBoost = false;
		if (no_paste2) prms2.doPasteBoost = false;
		if (no_step1)  prms1.offSet = 1;
		if (no_step2)  prms2.offSet = 1;

		if (order_inv) prms1.orderInvariance = true;
		if (order_inv) prms2.orderInvariance = true;
 
		VideoNLB::printNlbParameters(prms1);
		VideoNLB::printNlbParameters(prms2);

		return EXIT_SUCCESS;
	}

	// Declarations
	Video<float> noisy, basic, final;
	Video<float> fflow, bflow;

	// Load input videos
	                       noisy.loadVideo(input_path, first_frame, last_frame, frame_step);
	if (mode == DENO_ONLY) basic.loadVideo(inbsc_path, first_frame, last_frame, frame_step);
	if (use_oflow)         fflow.loadVideo(fflow_path, first_frame, last_frame, frame_step);
	if (use_oflow)         bflow.loadVideo(bflow_path, first_frame, last_frame, frame_step);

	// Denoising
	if (verbose) printf("Running BM4D-OF on the noisy video\n");

	// Compute denoising default parameters
	VideoNLB::nlbParams prms1, prms2;
	VideoNLB::initializeNlbParameters(prms1, 1, sigma, noisy.sz, verbose);
	VideoNLB::initializeNlbParameters(prms2, 2, sigma, noisy.sz, verbose);

	// Override with command line parameters
	if (space_search1 >= 0) VideoNLB::setSizeSearchWindow(prms1, (unsigned)space_search1);
	if (space_search2 >= 0) VideoNLB::setSizeSearchWindow(prms2, (unsigned)space_search2);
	if (patch_sizex1  >= 0) VideoNLB::setSizePatch(prms1, noisy.sz, (unsigned)patch_sizex1);;
	if (patch_sizex2  >= 0) VideoNLB::setSizePatch(prms2, noisy.sz, (unsigned)patch_sizex2);;
	if (num_patches1  >= 0) VideoNLB::setNSimilarPatches(prms1, (unsigned)num_patches1);
	if (num_patches2  >= 0) VideoNLB::setNSimilarPatches(prms2, (unsigned)num_patches2);
	if (beta1         >= 0) prms1.beta = beta1;
	if (beta2         >= 0) prms2.beta = beta2;
	if (beta1         >= 0) prms1.betaMean = beta1;
	if (beta2         >= 0) prms2.betaMean = beta2;
	if (beta_mean1    >= 0) prms1.betaMean = beta_mean1;
	if (beta_mean2    >= 0) prms2.betaMean = beta_mean2;

	if (patch_sizet1  >= 0) prms1.sizePatchTime = patch_sizet1;
	if (patch_sizet2  >= 0) prms2.sizePatchTime = patch_sizet2;

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

	if (order_inv) prms1.orderInvariance = true;
	if (order_inv) prms2.orderInvariance = true;

	// Percentage or processed groups of patches over total number of pixels
	std::vector<float> groupsRatio;

	// Run denoising algorithm
	groupsRatio = VideoNLB::runNlBayes(noisy, fflow, bflow, basic, final, prms1, prms2);

	if (verbose)
		printf("Done. Processed %5.2f%% of possible patch groups in 1st step, and\n"
		       "%5.2f%% in 2nd step.\n", groupsRatio[0], groupsRatio[1]);

	// Save output sequences
	if (verbose) printf("\nSaving output sequences\n");
	if (prms2.sizePatch) final.saveVideo(final_path, first_frame, frame_step);
	                     basic.saveVideo(basic_path, first_frame, frame_step);

	return EXIT_SUCCESS;
}
