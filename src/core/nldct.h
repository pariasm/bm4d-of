/*
 * Original work Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * Modified work Copyright (c) 2014, Pablo Arias <parias@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef VIDEO_NL_BAYES_H_INCLUDED
#define VIDEO_NL_BAYES_H_INCLUDED

#include "../utils/lib_videot.hpp"
#include "../utils/utilities.h"


//#define VBM3D_SEARCH             		// Use VBM3D predictive search
#define VBM3D_HAAR_TRANSFORM     		// Use VBM3D 3D transform
#define MC_PATCHES               		// Motion compensated 3D patches

namespace VideoNLB
{

enum ColorSpace {RGB, YUV};
enum Transform {dct, bior1_5};

/**
 * @brief Structures of parameters dedicated to NL-Bayes process
 *
 * @param sigma: value of noise;
 * @param sizePatch: size of patches (sizePatch x sizePatch);
 * @param nSimilarPatches: minimum number of similar patches wanted;
 * @param sizeSearchWindow: size of the search window around the reference patch;
 * @param boundary: must be > sizeSearchWindow. Boundary kept around sub-images when the image is
 *        subdivided for parallelization;
 * @param offSet: step between two similar patches;
 * @param useHomogeneousArea: if true, use the homogeneous area trick;
 * @param gamma: threshold to detect homogeneous area;
 * @param beta: parameter used to estimate the covariance matrix;
 * @param tau: parameter used to determine similar patches;
 * @param isFirstStep: true if the first step of the algorithm is wanted;
 * @param doPasteBoost: if true, patches near denoised similar patches will not be used as reference
 *        patches;
 * @param verbose: if true, print informations.
 **/
struct nlbParams
{
	float sigma;
	unsigned sizePatch;        // depends on sigma
	unsigned sizePatchTime;    // user given
	unsigned nSimilarPatches;  // depends on sigma, sizeSearchTimeRange (1 channel) or sizePatch (3 channels)
	unsigned sizeSearchWindow; // depends on nSimilarPatches
	unsigned sizeSearchTimeRangeFwd; // how many forward  frames in search cube
	unsigned sizeSearchTimeRangeBwd; // how many backward frames in search cube
	unsigned boundary;         // depends on sizeSearchWindow
	unsigned offSet;           // depends on sizePatch
	unsigned offSetTime;       // depends on sizePatchTime
	float beta;                // noise correction factor
	float betaMean;            // noise correction factor
	float tau;                 // depends on sizePatch
	bool orderInvariance;
	bool isFirstStep;
	bool doPasteBoost;
	bool verbose;
	ColorSpace colorSpace;

	Transform transform;           // VBM3D: patch transform
	unsigned sizeSearchWindowPred; // VBM3D: search window for predictive search
	unsigned nSimilarPatchesPred;  // VBM3D: number of sim patches for pred. search
	float dsub;                    // VBM3D: patch distance bias towards zero motion
	bool agg_window;               // VBM3D: aggregation window (for the moment only Kaiser)

};

/**
 * @brief Structure containing matrices used in the Bayesian estimation.
 *
 * @param groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 * @param baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 * @param covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 * @param covMatTmp: allocated memory. Used to process the Bayes estimate;
 * @param tmpMat: allocated memory. Used to process the Bayes estimate.
 **/
struct matWorkspace
{
	std::vector<float> groupTranspose;
	std::vector<float> groupTransposeNoisy;
	std::vector<float> baricenter;
	std::vector<float> baricenterNoisy;
	std::vector<float> covMat;
	std::vector<float> covMatTmp;
	std::vector<float> tmpMat;

	// the following have been added for computing a low rank
	// approximation of the covariance matrix C

	// used if low-rank approximation is done via eigendecomp. of C
	std::vector<float> covEigVecs;
	std::vector<float> covEigVals;

	// used if low-rank approximation of C is done via SVD of data matrix X
	std::vector<float> svd_U;      // left  sing. vecs of X ~ eigen vecs of C
	std::vector<float> svd_V;      // right sing. vecs of X
	std::vector<float> svd_VT;     // right sing. vecs of X (transposed for LAPACKE)
	std::vector<float> svd_S;      // sing. values of X ~ sqrt of eigen vecs of C
	std::vector<float> svd_work;
	std::vector<int  > svd_iwork;

	// double matrix to compute svd using idd
	std::vector<double> svd_dU;    // left sing. vecs of X ~ eigen vecs of C
	std::vector<double> svd_dV;    // right sing. vecs of X
	std::vector<double> svd_dS;    // sing. values of X ~ sqrt of eigen vecs of C
	std::vector<double> svd_ddata; //
	std::vector<double> svd_dwork;

	// store a fixed patch basis (eg. DCT)
	std::vector<float> patch_basis;
	std::vector<float> patch_basis_inv; // for non-orthogonal transforms

	// store a sepparable fixed patch basis (eg. DCT)
	std::vector<float> patch_basis_x;
	std::vector<float> patch_basis_y;
	std::vector<float> patch_basis_t;

	// aggregation window to minimize blocking artifacts
	std::vector<float> agg_window;
};

/**
 * @brief Initialize Parameters of the NL-Bayes algorithm.
 *
 * @param params   : will contain the nlbParams for the first step of the algorithm;
 * @param step     : select first or second step;
 * @param sigma    : standard deviation of the noise;
 * @param size     : size of the video;
 * @param verbose  : if true, print some informations.
 *
 * @return none.
 **/
void initializeNlbParameters(
	nlbParams &params
,	const unsigned step
,	const float sigma
,	const VideoSize &size
,	const bool verbose
);

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
void setSizeSearchWindow(
	nlbParams& prms
,	unsigned sizeSearchWindow
);

/**
 * @brief Sets size of the patch. It sets the pixel offset as half the patch
 * size (this is BM3D speed-up).
 *
 * @param prms      : nlbParams for first or second step of the algorithm;
 * @param sizePatch : size of the patch;
 *
 * @return none.
 **/
void setSizePatch(
	nlbParams& prms
,	const VideoSize &p_size
,	unsigned sizePatch
);

/**
 * @brief Sets number of similar patches, ensuring that the number of similar
 * patches is not larger that the total number of available patches.
 *
 * @param prms            : nlbParams for first or second step of the algorithm;
 * @param nSimilarPatches : number of similar patches;
 *
 * @return none.
 **/
void setNSimilarPatches(
	nlbParams& prms
,	unsigned nSimilarPatches
);

/**
 * @brief Sets the distance threshold relative to the patch size.
 *
 * @param prms : nlbParams for first or second step of the algorithm;
 * @param tau  : distance threshold;
 *
 * @return none.
 **/
void setTau(
	nlbParams& prms
,	const VideoSize &p_size
,	float tau
);

/**
 * @brief Display parameters of the NL-Bayes algorithm.
 *
 * @param i_params : nlbParams for first or second step of the algorithm;
 *
 * @return none.
 **/
void printNlbParameters(
	const nlbParams &i_params
);

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * @param i_imNoisy: contains the noisy video;
 * @param i_fflow  : forward optical flow;
 * @param i_bflow  : backward optical flow;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_params1 : parameters for first step
 * @param p_params1 : parameters for second step
 *
 * @return Percentage of processed groups over number of pixels.
 **/
std::vector<float> runNlBayes(
	Video<float> const& i_imNoisy
,  Video<float> const& i_fflow
,  Video<float> const& i_bflow
,	Video<float> &o_imBasic
,	Video<float> &o_imFinal
,	const nlbParams p_params1
,	const nlbParams p_params2
);

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy video;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_params: parameters of the method, contains:
 *    - sigma: standard deviation of the noise;
 *    - sizePatch: size of patches (sizePatch x sizePatch);
 *    - nSimilarPatches: number of similar patches;
 *    - sizeSearchWindow: size of the neighbourhood searching window;
 *    - useHomogeneousArea: if true, the trick of using homogeneous area will be used;
 *    - gamma: parameter used to determine if we are in an homogeneous area;
 *    - maxAvoid: parameter used to stop the paste trick;
 *    - beta: parameter used during the estimate of the denoised patch;
 *    - coefBaricenter: parameter to determine if the covariance matrix inversion is correct;
 *    - isFirstStep: true if it's the first step of the algorithm which is needed;
 *    - verbose: if true, print some informations, do nothing otherwise.
 *
 * @return Percentage of processed groups over number of pixels.
 **/
unsigned processNlBayes(
	Video<float> const& i_imNoisy
,	Video<float> const& i_fflow
,	Video<float> const& i_bflow
,	Video<float> &io_imBasic
,	Video<float> &o_imFinal
,	nlbParams const& p_params
,	VideoUtils::CropPosition p_crop = VideoUtils::CropPosition()
);

/**
 * @brief Keep from all near patches the similar ones to the reference patch for the second step.
 *
 * @param i_imNoisy: contains the original noisy video;
 * @param i_imBasic: contains the basic estimation;
 * @param o_groupNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_groupBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatches(
	Video<float> const& i_imNoisy
,	Video<float> const& i_imBasic
,	Video<float> const& i_fflow
,	Video<float> const& i_bflow
,	std::vector<float> &o_groupNoisy
,	std::vector<float> &o_groupBasic
,	std::vector<unsigned> &o_index
,	const unsigned p_ij
,	const nlbParams &p_params
);

/**
 * @brief Compute the Bayes estimation assuming a low rank covariance matrix.
 *
 * @param io_group: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *       - groupTranspose: allocated memory. Used to contain the transpose of io_groupNoisy;
 *       - baricenter: allocated memory. Used to contain the baricenter of io_groupBasic;
 *       - covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *       - covMatTmp: allocated memory. Used to process the Bayes estimate;
 *       - tmpMat: allocated memory. Used to process the Bayes estimate;
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
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	std::vector<float> &aggreWeights
);

/**
 * @brief Compute the Bayes estimation assuming a low rank covariance matrix.
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
 * @param aggreWeights: output aggregation weights.
 *
 * @return none.
 **/
float computeBayesEstimateStep2(
	std::vector<float> &io_groupNoisy
,	std::vector<float>  &i_groupBasic
,	matWorkspace &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams const& p_params
,	const unsigned p_nSimP
,	std::vector<float> &aggreWeights
);

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param im: update the video with estimate values;
 * @param weight: update corresponding weight, used later in the weighted aggregation;
 * @param mask: update values of mask: set to true the index of an used patch;
 * @param group: contains estimated values of all similar patches in the 3D group;
 * @param aggreWeights: input aggregation weights.
 * @param index: contains index of all similar patches contained in i_group;
 * @param params: see processStep2 for more explanation;
 * @param nSimP: number of similar patches.
 *
 * @return masked: number of processable pixels that were flaged non-processable.
 *
 **/
int computeAggregation(
	Video<float> &im
,	Video<float> &weight
,	Video<char>  &mask
,	std::vector<float> const& group
,	std::vector<float> const& aggreWeights
,	std::vector<float> const& aggreWindow
,	std::vector<unsigned> const& index
,	const nlbParams& params
,	const unsigned nSimP
);

/**
 * @brief Compute the final weighted aggregation.
 *
 * i_imReference: video of reference, when the weight if null;
 * io_imResult: will contain the final video;
 * i_weight: associated weight for each estimate of pixels.
 *
 * @return none.
 **/
void computeWeightedAggregation(
	Video<float> const& i_im
,	Video<float> &io_im
,	Video<float> const& i_weight
);

} // namespace

#endif // VIDEO_NL_BAYES_H_INCLUDED
