/*
 * Copyright (c) 2013, Marc Lebrun <marc.lebrun.ik@gmail.com>
 * All rights reserved.
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the GNU General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later
 * version. You should have received a copy of this license along
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef LIB_MATRIX_H_INCLUDED
#define LIB_MATRIX_H_INCLUDED

#include <vector>
#include <string>

/**
 * @brief Invert (in place) a symmetric real matrix, V -> Inv(V).
 *			The input matrix V is symmetric (V[i,j] = V[j,i]).
 *
 * @param io_mat = array containing a symmetric input matrix. This is converted to the inverse
			matrix;
 * @param p_N = dimension of the system (dim(v)=n*n)
 *
 * @return	EXIT_SUCCESS -> normal exit
 *			EXIT_FAILURE -> input matrix not positive definite
 **/
int inverseMatrix(
	std::vector<float> &io_mat
,	const unsigned p_N
);

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
	std::vector<float> &io_mat
,	const unsigned p_N
);

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
	std::vector<float> const& i_patches
,	std::vector<float> &o_covMat
,	const unsigned p_nb
,	const unsigned p_N
);

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
	std::vector<float> &o_mat
,	std::vector<float> const& i_A
,	std::vector<float> const& i_B
,	const unsigned p_n
,	const unsigned p_m
,	const unsigned p_l
);

/**
 * @brief Multiply two matrix A * B. It uses BLAS SGEMM.
 *
 * @param o_AB = array containing n by l product matrix at exit;
 * @param i_A = input array containing n by m matrix;
 * @param i_B = input array containing m by l matrix;
 * @param p_n, p_m, p_l = dimension parameters of arrays.
 * @param p_transA = true for transposing A.
 * @param p_transA = true for transposing B.
 * @param p_colMajor = true for if matrices should be read by columns.
 *
 * @return  none.
 **/
void productMatrix(
	std::vector<float> &o_AB
,	std::vector<float> const& i_A
,	std::vector<float> const& i_B
,	const unsigned p_n
,	const unsigned p_m
,	const unsigned p_l
,	const bool p_transA
,	const bool p_transB
,	const bool p_colMajor = true
,	unsigned lda = 0
,	unsigned ldb = 0
);

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
	std::vector<float> &i_mat
,	const unsigned p_n
,	const unsigned p_r
,	std::vector<float> &o_S
,	std::vector<float> &o_U
);

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
	std::vector<float> &i_mat
,	const unsigned p_n
,	const unsigned p_m
,	std::vector<float> &o_S
,	std::vector<float> &o_U
,	std::vector<float> &o_VT
,	std::vector<float> &i_work
,	std::vector<int> &i_iwork
);

/**
 * @brief Compute approximated low-rank SVD.
 *
 * NOTES:
 * - matrices are stored in column-major ordering
 * - columns of input matrices are contiguous in memory
 * - the output o_U contains the left singular vectors as columns, and is
 *   stored in column-major ordering (or as rows in row-major ordering)
 * - the output o_V contains the right singular vectors as columns, stored
 *   in column-major ordering (or as rows in row-major ordering)
 * - if the workspace vectors are empty, they are resized internally
 *
 * @param i_mat: m x n input matrix;
 * @param p_m  : rows of the matrix;
 * @param p_n  : cols of the matrix;
 * @param p_k  : rank;
 * @param o_S  : vector with min(m,n,k) singular values
 * @param o_U  : m x min(m,n,k) matrix with min(m,n,k) left singular values
 * @param o_V  : n x min(m,n,k) matrix with min(m,n,k) right singular values
 * @param i_work : id_dist's workspace
 *
 * @return none.
 **/
int matrixLRSVD(
	std::vector<double> &i_mat
,	int p_n
,	int p_m
,	int p_k
,	std::vector<double> &o_S
,	std::vector<double> &o_U
,	std::vector<double> &o_V
,	std::vector<double> &i_work
);

void printMatrix(
	std::vector<float> &matrix
,	unsigned rows
,	unsigned cols
,	std::string filename
);

#endif // LIB_MATRIX_H_INCLUDED
