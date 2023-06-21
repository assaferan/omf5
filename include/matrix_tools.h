/*****************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : matrix_tools.h
 *
 * Description: Functions for handling matrices.
 *
 *****************************************************************
 */

#ifndef __MATRIX_TOOLS_H__
#define __MATRIX_TOOLS_H__

// External packages dependencies

#include <antic/nf.h>
#include <antic/nf_elem.h>

#include <carat/matrix.h>

#include <flint/fmpz_mat.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>

// Internal dependencies

#include "isometry.h"
#include "square_matrix.h"

/***************************************************************************
 *
 * Function: init_sym_matrix
 *
 * Description: create a symmetric matrix from an array of integers,
 *              that describe the matrix according to a given input type.
 *              Possible input types are:
 *              "A" - the array is the upper triangular part of the matrix
 *                    in column major order (1st column, 2nd column,...)
 *              "GG" - the array is the upper triangular part of the matrix
 *                    in row major order (1st row, 2nd row, ....), with the
 *                    items on the diagonal being halved
 *              If a different input type is specified, returns NULL.
 *
 * Arguments:
 *     + coeff_vec (const int*) - the array of integers
 *     + inp_type (const char*) - the input type ("A" or "GG")
 *
 * Returns:
 *     + (matrix_TYP*) - the symmetric matrix
 *
 ***************************************************************************
 */

/* function to initialize a symmetric matrix from a vector of coefficients */
matrix_TYP* init_sym_matrix(const int* coeff_vec, const char* inp_type);

/***************************************************************************
 *
 * Function: print_mat
 *
 * Description: prints out a matrix nicely as a 2-d array to stdout.
 *
 * Arguments:
 *     + Q (const matrix_TYP*) - a matrix
 *
 ***************************************************************************
 */

/* function to print matrices nicely*/
void print_mat(const matrix_TYP* Q);

/***************************************************************************
 *
 * Function: print_mat_dense
 *
 * Description: prints out a matrix as an array of arrays, for
 *              loading in python or magma.
 *
 * Arguments:
 *     + Q (const matrix_TYP*) - a matrix
 *
 ***************************************************************************
 */

void print_mat_dense(const matrix_TYP* Q);

/***************************************************************************
 *
 * Function: swap
 *
 * Description: Modifies a matrix,
 *              by swapping two entries in the matrix.
 *              Assumes the first and the second entry are
 *              at different indices. 
 *
 * Arguments:
 *     + Q (int**) - a matrix, given as an array of arrays of integers
 *     + row1 (int) - row of first entry
 *     + col1 (int) - column of first entry
 *     + row2 (int) - row of second entry
 *     + col2 (int) - column of second entry
 *
 * Returns:
 *     + (int) - Returns 0. If (row1,col1) = (row2,col2) (illegal action),
 *               returns -1.
 *
 ***************************************************************************
 */

/* swapping, assumes they do not point to the same thing!!! */
int swap(int** Q, int row1, int col1, int row2, int col2);

/***************************************************************************
 *
 * Function: resymmetrize
 *
 * Description: Given a matrix of 64-bit integers,
 *              reflect its upper triangular part to the lower triangular
 *              part to create a symmetric matrix.
 *              (Useful in order to act only on the upper half, and then
 *               reflect, to save operations)
 *
 * Arguments:
 *     + Q (Z64**) - a matrix, given as an array of arrays of 64-bit
 *                   signed integers
 *
 ***************************************************************************
 */

/* resymmetrizing a matrix after working with upper triangular part */
void resymmetrize(Z64 **Q);

/***************************************************************************
 *
 * Function: automorphism_group
 *
 * Description: Returns the automorphism group of a lattice.
 *              This wraps some functionality of CARAT into a
 *              single function.
 *
 * Arguments:
 *     + Q (matrix_TYP*) - the gram matrix of a lattice L
 *
 * Returns:
 *     + (bravais_TYP*) - the automorphism group of L, Aut(L),
 *                        as a CARAT type
 *
 ***************************************************************************
 */

bravais_TYP* automorphism_group(matrix_TYP* Q);

// This function is obsolete for testing isometry, the new funciton is in isometry.h
// matrix_TYP* is_isometric(matrix_TYP* Q1, matrix_TYP* Q2);

/***************************************************************************
 *
 * Function: minkowski_reduce
 *
 * Description: Returns a Minkowski reduction of a symmetric p.d. matrix,
 *              i.e. a reduced representative, lying on the boundary of
 *              the Minkowski domain.
 *              This wraps functionality from CARAT into a single
 *              function.
 *
 * Arguments:
 *     + Q (matrix_TYP*) - a symmetric p.d. matrix
 *
 * Returns:
 *     + (matrix_TYP*) - a reduced symmetric p.d. matrix isometric to Q
 *
 ***************************************************************************
 */

matrix_TYP* minkowski_reduce(matrix_TYP* Q);

/***************************************************************************
 *
 * Function: greedy
 *
 * Description: Returns a greedy reduction of a symmetric p.d. matrix,
 *              according to the algorithm of [Nguyen, Stelhe].
 *              This does not always yield a Minkowski-reduced lattice,
 *              but it is good enough to test for isometries.
 *
 * Arguments:
 *     + gram (square_matrix_t) - a symmetric p.d. matrix
 *     + s (isometry_t) - a basis in which the lattice has this gram matrix
 *     + dim (int) - dimension of the lattice (up to 5), specifies dimension
 *                   of a submatrix to consider, i.e. we only reduce
 *                   the submatrix gram[:dim,:dim]
 *
 * Returns:
 *     + gram (square_matrix_t) - the reduced gram matrix
 *     + s (isometry_t) - a basis in which the lattice has this gram matrix
 *
 ***************************************************************************
 */

void greedy(square_matrix_t gram, isometry_t s, int dim);

/****************************************************************************
 *
 * Function: get_eigenvector_on_subspace
 *
 * Description: Computes an eigenvector of a rational matrix T on a linear
 *              invariant subspace W, indicated by basis vectors.
 *              Assumes the subspace W is invariant under T and irreducible. 
 *
 * Arguments:
 *     + T (const fmpq_mat_t) - a square matrix of rationals
 *     + basis_W (const fmpq_mat_t) - a matrix whose rows are the basis
 *                                    vectors of the subspace W
 *
 * Returns:
 *     + (bool) - true if W is irreducible, false otherwise
 *     + evec (nf_elem_t*) - an array of number field elements,
 *                           which is an eigenvector of T on W
 *     + nf (nf_t) - the number field over which evec is defined
 *
 ****************************************************************************
 */

bool get_eigenvector_on_subspace(nf_elem_t* evec, nf_t nf, const fmpq_mat_t T, const fmpq_mat_t basis_W);

void fmpq_mat_init_set_matrix_TYP(fmpq_mat_t M, const matrix_TYP* mat);
void fmpz_mat_init_set_matrix_TYP(fmpz_mat_t M, const matrix_TYP* mat);
void matrix_TYP_init_set_fmpz_mat(matrix_TYP** mat, const fmpz_mat_t M);
void nmod_mat_init_set_fmpz_mat(nmod_mat_t dest, const fmpz_mat_t mat, mp_limb_t n);
void fq_nmod_mat_init_set_fmpz_mat(fq_nmod_mat_t dest, const fmpz_mat_t mat, const fq_nmod_ctx_t F);

void restrict_mat(fmpq_mat_t res_T, const fmpq_mat_t T, const fmpq_mat_t basis_W);

void kernel_on(fmpq_mat_t ker, const fmpq_mat_t A, const fmpq_mat_t B);

void fmpq_mat_kernel(fmpq_mat_t ker, const fmpq_mat_t mat);
void fmpq_mat_left_kernel(fmpq_mat_t ker, const fmpq_mat_t mat);

#endif // __MATRIX_TOOLS_H__
