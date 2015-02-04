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

#include "../Utilities/LibVideo.h"
#include "../Utilities/Utilities.h"

namespace VideoNLB
{

/**
 * @brief Structures of parameters dedicated to NL-Bayes process
 *
 * @param sigma: value of noise;
 * @param sizePatch: size of patches (sizePatch x sizePatch);
 * @param nSimilarPatches: minimum number of similar patches wanted;
 * @param sizeSearchWindow: size of the search window around the reference patch;
 * @param boundary: must be > sizeSearchWindow. Boundary kept around sub-images when the image is
 *      subdivided for parallelization;
 * @param offSet: step between two similar patches;
 * @param useHomogeneousArea: if true, use the homogeneous area trick;
 * @param gamma: threshold to detect homogeneous area;
 * @param beta: parameter used to estimate the covariance matrix;
 * @param tau: parameter used to determine similar patches;
 * @param isFirstStep: true if the first step of the algorithm is wanted;
 * @param doPasteBoost: if true, patches near denoised similar patches will not be used as reference
 *		patches;
 * @param verbose: if true, print informations.
 **/
struct nlbParams
{
	float sigma;
	unsigned sizePatch;
	unsigned nSimilarPatches;
	unsigned sizeSearchWindow;
	unsigned sizeSearchTimeRangeFwd; //< VIDEO how many forward  frames in search cube
	unsigned sizeSearchTimeRangeBwd; //< VIDEO how many backward frames in search cube
	unsigned boundary;
	unsigned offSet; //< ASK MARC any experiment quantifying implact of this parameter
	bool useHomogeneousArea;
	float gamma;
	float beta;
	float tau;
	bool isFirstStep;
	bool doPasteBoost;
	bool verbose;
};

/**
 * @brief Structure containing usefull matrices for the Bayes estimations.
 *
 * @param group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 * @param baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 * @param covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 * @param covMatTmp: allocated memory. Used to process the Bayes estimate;
 * @param tmpMat: allocated memory. Used to process the Bayes estimate.
 **/
struct matParams
{
	std::vector<float> group3dTranspose;
	std::vector<float> baricenter;
	std::vector<float> covMat;
	std::vector<float> covMatTmp;
	std::vector<float> tmpMat;
};

/**
 * @brief Initialize Parameters of the NL-Bayes algorithm.
 *
 * @param o_paramStep1 : will contain the nlbParams for the first step of the algorithm;
 * @param o_paramStep2 : will contain the nlbParams for the second step of the algorithm;
 * @param p_sigma : standard deviation of the noise;
 * @param p_imSize: size of the image;
 * @param p_useArea1 : if true, use the homogeneous area trick for the first step;
 * @param p_useArea2 : if true, use the homogeneous area trick for the second step;
 * @param p_verbose : if true, print some informations.
 *
 * @return none.
 **/
void initializeNlbParameters(
	nlbParams &o_paramStep1
,	nlbParams &o_paramStep2
,	const float p_sigma
,	const VideoSize &p_imSize
,	const bool p_useArea1
,	const bool p_useArea2
,	const bool p_verbose
);

/**
 * @brief Display parameters of the NL-Bayes algorithm.
 *
 * @param i_paramStep1 : nlbParams for the first step of the algorithm;
 * @param i_paramStep2 : nlbParams for the second step of the algorithm;
 *
 * @return none.
 **/
void printNlbParameters(
	nlbParams &i_paramStep1
,	nlbParams &i_paramStep2
);

/**
 * @brief Main function to process the whole NL-Bayes algorithm.
 *
 * @param i_imNoisy: contains the noisy video;
 * @param o_imBasic: will contain the basic estimate image after the first step;
 * @param o_imFinal: will contain the final denoised image after the second step;
 * @param p_useArea1 : if true, use the homogeneous area trick for the first step;
 * @param p_useArea2 : if true, use the homogeneous area trick for the second step;
 * @param p_sigma : standard deviation of the noise;
 * @param p_verbose : if true, print some informations.
 *
 * @return EXIT_FAILURE if something wrong happens during the whole process.
 **/
int runNlBayes(
	Video_f32 const& i_imNoisy
,	Video_f32 &o_imBasic
,	Video_f32 &o_imFinal
,	const bool p_useArea1
,	const bool p_useArea2
,	const float p_sigma
,	const bool p_verbose
);

/**
 * @brief Generic step of the NL-Bayes denoising (could be the first or the second).
 *
 * @param i_imNoisy: contains the noisy video;
 * @param io_imBasic: will contain the denoised image after the first step (basic estimation);
 * @param o_imFinal: will contain the denoised image after the second step;
 * @param p_params: parameters of the method, contains:
 *			- sigma: standard deviation of the noise;
 *			- sizePatch: size of patches (sizePatch x sizePatch);
 *			- nSimilarPatches: number of similar patches;
 *			- sizeSearchWindow: size of the neighbourhood searching window;
 *			- useHomogeneousArea: if true, the trick of using homogeneous area will be used;
 *			- gamma: parameter used to determine if we are in an homogeneous area;
 *			- maxAvoid: parameter used to stop the paste trick;
 *			- beta: parameter used during the estimate of the denoised patch;
 *			- coefBaricenter: parameter to determine if the covariance matrix inversion is correct;
 *			- isFirstStep: true if it's the first step of the algorithm which is needed;
 *			- verbose: if true, print some informations, do nothing otherwise.
 *
 * @return none.
 **/
void processNlBayes(
	Video_f32 const& i_imNoisy
,	Video_f32 &io_imBasic
,	Video_f32 &o_imFinal
,	nlbParams &p_params
);

/**
 * @brief Estimate the best similar patches to a reference one.
 *
 * @param i_im: contains the noisy video on which distances are processed;
 * @param o_group3d: will contain values of similar patches;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
void estimateSimilarPatchesStep1(
	Video_f32 const& i_im
,	std::vector<std::vector<float> > &o_group3d
,	std::vector<unsigned> &o_index
,	const unsigned p_ij
,	const nlbParams &p_params
);

/**
 * @brief Keep from all near patches the similar ones to the reference patch for the second step.
 *
 * @param i_imNoisy: contains the original noisy video;
 * @param i_imBasic: contains the basic estimation;
 * @param o_group3dNoisy: will contain similar patches for all channels of i_imNoisy;
 * @param o_group3dBasic: will contain similar patches for all channels of i_imBasic;
 * @param o_index: will contain index of similar patches;
 * @param p_ij: index of the reference patch;
 * @param p_params: see processStep2 for more explanations.
 *
 * @return number of similar patches kept.
 **/
unsigned estimateSimilarPatchesStep2(
	Video_f32 const& i_imNoisy
,	Video_f32 const& i_imBasic
,	std::vector<float> &o_group3dNoisy
,	std::vector<float> &o_group3dBasic
,	std::vector<unsigned> &o_index
,	const unsigned p_ij
,	const nlbParams &p_params
);

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group3d: contains for each channels values of similar patches. If an homogeneous area
 *			is detected, will contain the average of all pixels in similar patches;
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_imSize: size of the video.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep1(
	std::vector<std::vector<float> > &io_group3d
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
);

/**
 * @brief Detect if we are in an homogeneous area. In this case, compute the mean.
 *
 * @param io_group3dNoisy: contains values of similar patches for the noisy video;
 * @param io_group3dBasic: contains values of similar patches for the basic video. If an homogeneous
 *		area is detected, will contain the average of all pixels in similar patches;
 * @param p_sP2: size of each patch (sP x sP);
 * @param p_nSimP: number of similar patches;
 * @param p_threshold: threshold below which an area is declared homogeneous;
 * @param p_imSize: size of the video.
 *
 * @return 1 if an homogeneous area is detected, 0 otherwise.
 **/
int computeHomogeneousAreaStep2(
	std::vector<float> const& i_group3dNoisy
,	std::vector<float> &io_group3dBasic
,	const unsigned p_sP
,	const unsigned p_nSimP
,	const float p_threshold
,	const VideoSize &p_imSize
);

/**
 * @brief Compute the Bayes estimation.
 *
 * @param io_group3d: contains all similar patches. Will contain estimates for all similar patches;
 * @param i_mat: contains :
 *		- group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *		- baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 *		- covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *		- covMatTmp: allocated memory. Used to process the Bayes estimate;
 *		- tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
void computeBayesEstimateStep1(
	std::vector<std::vector<float> > &io_group3d
,	matParams &i_mat
,	unsigned &io_nInverseFailed
,	nlbParams &p_params
);

/**
 * @brief Compute the Bayes estimation.
 *
 * @param i_group3dNoisy: contains all similar patches in the noisy video;
 * @param io_group3dBasic: contains all similar patches in the basic video. Will contain estimates
 *			for all similar patches;
 * @param i_mat: contains :
 *		- group3dTranspose: allocated memory. Used to contain the transpose of io_group3dNoisy;
 *		- baricenter: allocated memory. Used to contain the baricenter of io_group3dBasic;
 *		- covMat: allocated memory. Used to contain the covariance matrix of the 3D group;
 *		- covMatTmp: allocated memory. Used to process the Bayes estimate;
 *		- tmpMat: allocated memory. Used to process the Bayes estimate;
 * @param io_nInverseFailed: update the number of failed matrix inversion;
 * @param p_imSize: size of the video;
 * @param p_params: see processStep2 for more explanations;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeBayesEstimateStep2(
	std::vector<float> &i_group3dNoisy
,	std::vector<float> &io_group3dBasic
,	matParams &i_mat
,	unsigned &io_nInverseFailed
,	const VideoSize &p_imSize
,	nlbParams p_params
,	const unsigned p_nSimP
);

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the video with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_params: see processStep1 for more explanation.
 *
 * @return none.
 **/
void computeAggregationStep1(
	Video_f32 &io_im
,	Video_f32 &io_weight
,	std::vector<bool> &io_mask
,	std::vector<std::vector<float> > const& i_group3d
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
);

/**
 * @brief Aggregate estimates of all similar patches contained in the 3D group.
 *
 * @param io_im: update the video with estimate values;
 * @param io_weight: update corresponding weight, used later in the weighted aggregation;
 * @param io_mask: update values of mask: set to true the index of an used patch;
 * @param i_group3d: contains estimated values of all similar patches in the 3D group;
 * @param i_index: contains index of all similar patches contained in i_group3d;
 * @param p_params: see processStep2 for more explanation;
 * @param p_nSimP: number of similar patches.
 *
 * @return none.
 **/
void computeAggregationStep2(
	Video_f32 &io_im
,	Video_f32 &io_weight
,	std::vector<bool> &io_mask
,	std::vector<float> const& i_group3d
,	std::vector<unsigned> const& i_index
,	const nlbParams &p_params
,	const unsigned p_nSimP
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
	Video_f32 const& i_im
,	Video_f32 &io_im
,	Video_f32 const& i_weight
);

} // namespace

#endif // VIDEO_NL_BAYES_H_INCLUDED