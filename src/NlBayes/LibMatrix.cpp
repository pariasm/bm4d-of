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
 * @file LibMatrix.cpp
 * @brief Tools for matrix manipulation, based on ccmath functions
 *        by Daniel A. Atkinson.
 *
 * @author Marc Lebrun <marc.lebrun.ik@gmail.com>
 **/

#include "LibMatrix.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include <cblas.h>
#include <lapacke.h>

using namespace std;

/**
 * @brief Invert (in place) a symmetric real matrix, V -> Inv(V).
 *			The input matrix V is symmetric (V[i,j] = V[j,i]).
 *
 * @param io_mat : array containing a symmetric input matrix. This is converted to the inverse
			matrix;
 * @param p_N : dimension of the system (dim(v)=n*n)
 *
 * @return	EXIT_SUCCESS -> normal exit
 *			EXIT_FAILURE -> input matrix not positive definite
 **/
int inverseMatrix(
	vector<float> &io_mat
,	const unsigned p_N
){
	//! Initializations
	float z;
	unsigned p, q, r, s, t, j, k;

	for (j = 0, p = 0; j < p_N ; j++, p += p_N + 1) {
		for (q = j * p_N; q < p ; q++)
			io_mat[p] -= io_mat[q] * io_mat[q];

		if (io_mat[p] <= 0.f) return EXIT_FAILURE;

		io_mat[p] = sqrtf(io_mat[p]);

		for (k = j + 1, q = p + p_N; k < p_N ; k++, q += p_N) {
			for (r = j * p_N, s = k * p_N, z = 0.f; r < p; r++, s++) {
				z += io_mat[r] * io_mat[s];
			}

			io_mat[q] -= z;
			io_mat[q] /= io_mat[p];
		}
	}

	transposeMatrix(io_mat, p_N);

	for (j = 0, p = 0; j < p_N; j++, p += p_N + 1) {
		io_mat[p] = 1.f / io_mat[p];

		for (q = j, t = 0; q < p; t += p_N + 1, q += p_N) {
			for (s = q, r = t, z = 0.f; s < p; s += p_N, r++) {
				z -= io_mat[s] * io_mat[r];
			}

			io_mat[q] = z * io_mat[p];
		}
	}

	for (j = 0, p = 0; j < p_N; j++, p += p_N + 1) {
		for (q = j, t = p - j; q <= p; q += p_N, t++) {
			for (k = j, r = p, s = q, z = 0.f; k < p_N; k++, r++, s++)
				z += io_mat[r] * io_mat[s];

			io_mat[t] = io_mat[q] = z;
		}
	}

	return EXIT_SUCCESS;
}

/**
 * @brief Transpose a real square matrix in place io_mat -> io_mat~.
 *
 * @param io_mat : pointer to an array of p_N by p_N input matrix io_mat. This is overloaded by
			the transpose of io_mat;
 * @param p_N : dimension (dim(io_mat) = p_N * p_N).
 *
 * @return none.
 **/
void transposeMatrix(
	vector<float> &io_mat
,	const unsigned p_N
){
	for (unsigned i = 0; i < p_N - 1; i++) {
		unsigned p = i * (p_N + 1) + 1;
		unsigned q = i * (p_N + 1) + p_N;

		for (unsigned j = 0; j < p_N - 1 - i; j++, p++, q += p_N) {
			const float s = io_mat[p];
			io_mat[p] = io_mat[q];
			io_mat[q] = s;
		}
	}
}

/**
 * @brief Compute the covariance matrix.
 *
 * @param i_patches: set of patches of size (nb x N);
 * @param o_covMat: will contain patches' * patches;
 * @param p_N : size of a patch;
 * @param p_nb: number of similar patches in the set of patches.
 *
 * @return none.
 **/
void covarianceMatrix(
	vector<float> const& i_patches
,	vector<float> &o_covMat
,	const unsigned p_nb
,	const unsigned p_N
){
	const float coefNorm = 1.f / (float) (p_nb);

	for (unsigned i = 0; i < p_N; i++)
	for (unsigned j = 0; j < i + 1; j++) {
		float val = 0.f;
		for (unsigned k = 0; k < p_nb; k++)
			val += i_patches[i * p_nb + k] * i_patches[j * p_nb + k];

		o_covMat[i * p_N + j] = val * coefNorm;
		o_covMat[i + j * p_N] = val * coefNorm;
	}
}

/**
 * @brief Multiply two matrix A * B. (all matrices stored in row order).
 *
 * @param o_mat = array containing n by l product matrix at exit;
 * @param i_A = input array containing n by m matrix;
 * @param i_B = input array containing m by l matrix;
 * @param p_n, p_m, p_l = dimension parameters of arrays.
 *
 * @return  none.
 **/
void productMatrix(
	vector<float> &o_mat
,	vector<float> const& i_A
,	vector<float> const& i_B
,	const unsigned p_n
,	const unsigned p_m
,	const unsigned p_l
){
	vector<float> q0(p_m, 0.f);

	for (unsigned i = 0; i < p_l; i++) {
		for (unsigned k = 0; k < p_m; k++) {
			q0[k] = i_B[i + k * p_l];
		}

		for (unsigned j = 0; j < p_n; j++) {
			float z = 0.f;

			for (unsigned k = 0; k < p_m; k++) {
				z += i_A[j * p_m + k] * q0[k];
			}

			o_mat[i + j * p_l] = z;
		}
	}
}

/**
 * @brief Multiply two matrix A * B. It uses BLAS SGEMM.
 *
 * @param o_AB = array containing nxl product matrix at exit;
 * @param i_A = input array containing nxm (or mxn if p_transA) matrix;
 * @param i_B = input array containing mxl (or lxm if p_transB) matrix;
 * @param p_n, p_l, p_m = dimension parameters of arrays.
 * @param p_transA = true for transposing A.
 * @param p_transA = true for transposing B.
 * @param p_colMajor = true for if matrices should be read by columns.
 *
 * @return  none.
 **/
void productMatrix(
	vector<float> &o_AB
,	vector<float> const& i_A
,	vector<float> const& i_B
,	const unsigned p_n
,	const unsigned p_l
,	const unsigned p_m
,	const bool p_transA
,	const bool p_transB
,	const bool p_colMajor
,	unsigned lda
,	unsigned ldb
){
	if (!lda) lda = p_colMajor ? (p_transA ? p_m : p_n) : (p_transA ? p_n : p_m);
	if (!ldb) ldb = p_colMajor ? (p_transB ? p_l : p_m) : (p_transB ? p_m : p_l);

	cblas_sgemm(p_colMajor ? CblasColMajor : CblasRowMajor,  // matrix storage mode
	            p_transA   ? CblasTrans    : CblasNoTrans,   // op(A)
	            p_transB   ? CblasTrans    : CblasNoTrans,   // op(B)
	            p_n,                                         // rows(op(A)) [= rows(AB)   ]
	            p_l,                                         // cols(op(B)) [= cols(AB)   ]
	            p_m,                                         // cols(op(A)) [= rows(op(B))]
	            1.f,                                         // alpha
	            i_A.data(), lda,                             // A, lda
	            i_B.data(), ldb,                             // B, ldb
	            0.f,                                         // beta
	            o_AB.data(), p_colMajor ? p_n : p_l          // AB, ldab
	            );
}

/**
 * @brief Compute a specified number of eigenvectors and eigenvalues of a
 * symmetric matrix.
 *
 * NOTES:
 * - matrices are stored in column-major ordering
 * - columns of input matrices are contiguous in memory
 * - only the upper triangular triangular part of o_mat is used
 * - the upper triangular part of o_mat is destroyed
 * - the output o_U contains the eigenvectors as columns, and is
 *   stored ini column-major ordering (i.e. it returns the eigenvectors
 *   as rows in row-major ordering)
 *
 * @param i_mat: contains input matrix;
 * @param p_n  : size of the matrix;
 * @param p_r  : number of eigenvectors and eigenvalues.
 * @param o_S  : vector with the r eigenvalues
 * @param o_U  : matrix with the r eigenvectors
 *
 * @return none.
 **/
int matrixEigs(
	vector<float> &i_mat
,	const unsigned p_n
,	const unsigned p_r
,	vector<float> &o_S
,	vector<float> &o_U
){
	// set parameters for LAPACKE SSYEV function
	// SSYEVX: Single SYmmetric EigenValues and eigenvectors eXpert
	lapack_int m;          //< total values of eigenvalues found

	o_S.resize(p_n);       //< array of dimension n. The first m entries
	                       //< contain the eigenvalues found in ascending order.

	o_U.resize(p_n*p_r);   //< the first m columns contain the output eigenvectors

	lapack_int lda = p_n;  //< leading dimension for input matrix
	lapack_int ldu = p_n;  //< leading dimension for output eigenvectors matrix

	lapack_int ifail[p_n]; //< if jobz == 'V' and info > 0, then ifail contains
	                       //< the indices of the eigenvectors that didn't converge
	
	lapack_int info;       //< info =  0 : successful exit
	                       //< info = -i : ith argument is wrong
	                       //< info =  i : i eigenvectors failed to converge

	info = LAPACKE_ssyevx(LAPACK_COL_MAJOR,
			'V',                // compute eigenvalues and eigenvectors
			'I',                // range 'I': eigenvals/vecs between IL-th and IU-th
			'U',                // use upper triangular part of A
			p_n,                // order of matrix A
			i_mat.data(),       // matrix A
			lda,                // stride of matrix
			-1, -1,             // not used (used only when range is 'V'
			p_n - p_r + 1, p_n, // IL and IU indices for eigenvals/vecs range
			0,                  // abstol for stopping criterion
			&m,                 // total values of eigenvectors found
			o_S.data(),         // eigenvalues output
			o_U.data(),         // eigenvectors matrix output
			ldu,                // eigenvectors matrix stride
			ifail               // eigenvectors that did not converge
			);

	return(info);
}

/**
 * @brief Compute a complete SVD.
 *
 * NOTES:
 * - matrices are stored in column-major ordering
 * - columns of input matrices are contiguous in memory
 * - the output o_U contains the left singular vectors as columns, and is
 *   stored in column-major ordering (or as rows in row-major ordering)
 * - the output o_VT contains the right singular vectors as rows, stored 
 *   in column-major ordering (or as columns in row-major ordering)
 * - if the workspace vectors are empty, they are resized internally
 *
 * @param i_mat: m x n input matrix;
 * @param p_n  : cols of the matrix;
 * @param p_m  : rows of the matrix;
 * @param o_S  : vector with min(m,n) singular values
 * @param o_U  : m x min(m,n) matrix with min(m,n) left singular values
 * @param o_VT : min(m,n) x n matrix with min(m,n) right singular values (transposed)
 * @param i_work  : LAPACK's workspace
 * @param i_iwork : LAPACK's integer workspace
 *
 * @return none.
 **/
int matrixSVD(
	vector<float> &i_mat
,	const unsigned p_n
,	const unsigned p_m
,	vector<float> &o_S
,	vector<float> &o_U
,	vector<float> &o_VT
,	vector<float> &i_work
,	vector<int> &i_iwork
){
	// set parameters for LAPACKE SGESDD function
	// SGESDD: Single GEneric Singular value Decomposition using Divide-and-conquer

	o_S .resize(std::min(p_m, p_n));      // singular values in descending order.
	o_U .resize(p_m*std::min(p_m, p_n));  // columns contain left  sigular vectors
	o_VT.resize(p_n*std::min(p_m, p_n));  // rows    contain right sigular vectors
	i_iwork.resize(8*std::min(p_m, p_n)); // integer workspace for SGESDD

	lapack_int info;       // info =  0 : successful exit
	                       // info = -i : ith argument is wrong
	                       // info =  i : failed to converge

	if (!i_work.size())
	{
		float workspace_size;
		// query workspace
		info = LAPACKE_sgesdd_work(LAPACK_COL_MAJOR,
			'S',                   // economic SVD
			p_m,                   // rows of matrix A
			p_n,                   // cols of matrix A
			i_mat.data(),          // matrix A
			p_m,                   // stride of matrix A
			o_S.data(),            // singular values output
			o_U.data(),            // matrix U output
			p_m,                   // stride of matrix U
			o_VT.data(),           // matrix VT output
			std::min(p_m, p_n),    // stride of VT
			&workspace_size,       // size of workspace (output)
			-1,                    // because we are quering
			i_iwork.data()         // integer workspace
			);

		// allocate workspace
		i_work.resize((int)workspace_size);
	}

	info = LAPACKE_sgesdd_work(LAPACK_COL_MAJOR,
		'S',                   // economic SVD
		p_m,                   // rows of matrix A
		p_n,                   // cols of matrix A
		i_mat.data(),          // matrix A
		p_m,                   // stride of matrix A
		o_S.data(),            // singular values output
		o_U.data(),            // matrix U output
		p_m,                   // stride of matrix U
		o_VT.data(),           // matrix VT output
		std::min(p_m, p_n),    // stride of VT
		i_work.data(),         // workspace
		i_work.size(),         // workspace size
		i_iwork.data()         // integer workspace
		);

	return(info);
}

void printMatrix(
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
