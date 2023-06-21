/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : isometry.h
 *
 * Description: Data structure for isometries.
 *
 **********************************************************
 */

#ifndef __ISOMETRY_H__
#define __ISOMETRY_H__

// Required packages dependencies

#include <carat/matrix.h>
#include <flint/fmpz_mat.h>

// Self dependencies

#include "square_matrix.h"
#include "typedefs.h"

/***********************************************************************
 *
 * Type: isometry_t
 *
 * Description: this struct holds data for an isometry
 *              between lattices of rank 5,
 *              which can be in general a Q-isometry.
 *
 * Fields:
 *     + s (square_matrix_t) - a 5x5 matrix representing the isometry
 *     + s_inv (square_matrix_t) - a 5x5 matrix represting the inverse
 *                                 isometry
 *     + denom (int) - denominator of s (so the isometry is really
 *                     the rational matrix (s / denom))
 *     + inv_denom (int) - denominator of _invs (so the inverse
 *                         isometry is really
 *                         the rational matrix (s_inv / inv_denom))
 *
 ***********************************************************************
 */

// This structure keeps all the isometry information (denominator and inverse) in one place

typedef struct
{
  square_matrix_t s;
  square_matrix_t s_inv;
  int denom;
  int inv_denom;
  
} isometry_struct;

typedef isometry_struct isometry_t[1];

/**********************************************************************
 *
 * Function: isometry_init
 *
 * Description: Initializes an isometry.
 *
 * Returns:
 *     + isom (isometry_t) - the identity isometry.
 *
 **********************************************************************
 */

void isometry_init(isometry_t isom);

/**********************************************************************
 *
 * Function: isometry_init_set
 *
 * Description: Initializes an isometry, and sets it to a given one.
 *
 * Arguments:
 *     + src (isometry_t) - the given isometry.
 *
 * Returns:
 *     + dest (isometry_t) - the cloned isometry.
 *
 **********************************************************************
 */

void isometry_init_set(isometry_t dest, const isometry_t src);

/**********************************************************************
 *
 * Function: isometry_init_set_square_matrix
 *
 * Description: Initializes an isometry, and sets it to a given
 *              rational matrix.
 *
 * Arguments:
 *     + s (square_matrix_t) - an integral 5x5 matrix
 *     + denom (int) - the denominator
 *
 * Returns:
 *     + isom (isometry_t) - the isometry s / denom.
 *
 **********************************************************************
 */

void isometry_init_set_square_matrix(isometry_t isom, const square_matrix_t s, int denom);

/**********************************************************************
 *
 * Function: isometry_init_set_fmpz_mat
 *
 * Description: Initializes an isometry, and sets it to a given
 *              rational matrix.
 *
 * Arguments:
 *     + s (fmpz_mat_t) - an integral 5x5 matrix
 *     + denom (int) - the denominator
 *
 * Returns:
 *     + isom (isometry_t) - the isometry s / denom.
 *
 **********************************************************************
 */

void isometry_init_set_fmpz_mat(isometry_t, const fmpz_mat_t s, int denom);

/**********************************************************************
 *
 * Function: isometry_mul
 *
 * Description: Computes the composition (product) of two isometries.
 *
 * Arguments:
 *     + sL (const isometry_t) - the left operand
 *     + sR (const isometry_t) - the right operand
 *
 * Returns:
 *     + prod (isometry_t) - the isometry sL * sR
 *
 **********************************************************************
 */

void isometry_mul(isometry_t prod, const isometry_t sL, const isometry_t sR);

/**********************************************************************
 *
 * Function: isometry_muleq_right
 *
 * Description: Multiply an isometry in place from the right.
 *              (sL := sL*sR;)
 *
 * Arguments:
 *     + sL (isometry_t) - the left operand
 *     + sR (const isometry_t) - the right operand
 *
 * Returns:
 *     + sL (isometry_t) - the isometry sL * sR
 *
 **********************************************************************
 */

void isometry_muleq_right(isometry_t sL, const isometry_t sR);

/**********************************************************************
 *
 * Function: isometry_muleq_right
 *
 * Description: Multiply an isometry in place from the left.
 *              (sR := sL*sR;)
 *
 * Arguments:
 *     + sR (isometry_t) - the right operand
 *     + sL (const isometry_t) - the left operand
 *
 * Returns:
 *     + sR (isometry_t) - the isometry sL * sR
 *
 **********************************************************************
 */

void isometry_muleq_left(isometry_t sR, const isometry_t sL);

/************************************************************************
 *
 * Function: isometry_mul_mat_left
 *
 * Description: Multiply an isometry from the left by a rational matrix.
 *              (computes (sL / sL_denom) * sR)
 *
 * Arguments:
 *     + sL (const square_matrix_t) - the numerator of the left operand
 *     + sL_denom (int) - the denominator of the left operand
 *     + sR (const isometry_t) - the right operand
 *
 * Returns:
 *     + prod (isometry_t) - the isometry (sL / sL_denom) * sR
 *
 ************************************************************************
 */

void isometry_mul_mat_left(isometry_t prod, const square_matrix_t sL, int sL_denom, const isometry_t sR);

/************************************************************************
 *
 * Function: isometry_mul_mat_right
 *
 * Description: Multiply an isometry from the right by a rational matrix.
 *              (computes sL * (sR / sR_denom))
 *
 * Arguments:
 *     + sL (const isometry_t) - the left operand
 *     + sR (const square_matrix_t) - the numerator of the right operand
 *     + sR_denom (int) - the denominator of the right operand
 *
 * Returns:
 *     + prod (isometry_t) - the isometry sL * (sR / sR_denom)
 *
 ************************************************************************
 */

void isometry_mul_mat_right(isometry_t prod, const isometry_t sL, const square_matrix_t sR, int sR_denom);

/************************************************************************
 *
 * Function: isometry_mul_mateq_left
 *
 * Description: Multiply an isometry from the left by a rational matrix,
 *              in place (sR = (sL / sL_denom) * sR)
 *
 * Arguments:
 *     + sR (isometry_t) - the right operand
 *     + sL (const square_matrix_t) - the numerator of the left operand
 *     + sL_denom (int) - the denominator of the left operand
 *
 * Returns:
 *     + sR (isometry_t) - the isometry (sL / sL_denom) * sR
 *
 ************************************************************************
 */

void isometry_mul_mateq_left(isometry_t sR, const square_matrix_t sL, int sL_denom);

/************************************************************************
 *
 * Function: isometry_muleq_mat_right
 *
 * Description: Multiply an isometry from the right by a rational matrix,
 *              in place (sL = sL * (sR / sR_denom))
 *
 * Arguments:
 *     + sL (isometry_t) - the left operand
 *     + sR (const square_matrix_t) - the numerator of the right operand
 *     + sR_denom (int) - the denominator of the right operand
 *
 * Returns:
 *     + sL (isometry_t) - the isometry sL * (sR / sR_denom)
 *
 ************************************************************************
 */

void isometry_mul_mateq_right(isometry_t sL, const square_matrix_t sR, int sR_denom);

/************************************************************************
 *
 * Function: isometry_inv
 *
 * Description: Computes the inverse of an isometry.
 *
 * Arguments:
 *     + isom (const isometry_t) - a given isometry
 *
 * Returns:
 *     + inv (isometry_t) - the isometry isom^(-1)
 *
 ************************************************************************
 */

void isometry_inv(isometry_t inv, const isometry_t isom);

/************************************************************************
 *
 * Function: isometry_transform_gram
 *
 * Description: Given the gram matrix of an integral lattice L,
 *              computes the gram matrix of the image of the lattice
 *              under an isometry g.
 *              Assumes that the lattice g(L) is also integral.
 *
 * Arguments:
 *     + g (const isometry_t) - the given isometry
 *     + Q (const square_matrix_t) - the gram matrix of the lattice L
 *
 * Returns:
 *     + gtQg (square_matrix_t) - the gram matrix of g(L), namely
 *                                the matrix g^t * Q * g.
 *
 ************************************************************************
 */

void isometry_transform_gram(square_matrix_t gtQg, const isometry_t g, const square_matrix_t Q);

/************************************************************************
 *
 * Function: isometry_is_isom
 *
 * Description: Given the gram matrix of two integral lattices L1, L2,
 *              and an isometry s, checks whether s(L1) = L2.
 *
 * Arguments:
 *     + s (isometry_t) - an isometry
 *     + q1 (const square_matrix_t) - the gram matrix of the lattice L1
 *     + q2 (const square_matrix_t) - the gram matrix of the lattice L2
 *
 * Returns:
 *     + (bool) - true if s(L1) = L2, i.e s^t * q1 * s = q2
 *
 ************************************************************************
 */

bool isometry_is_isom(isometry_t s, const square_matrix_t q1, const square_matrix_t q2);

/************************************************************************
 *
 * Function: isometry_clear
 *
 * Description: Clears the memory allocated for an isometry.
 *
 * Arguments:
 *     + isom (isometry_t) - an isometry
 *
 ************************************************************************
 */

void isometry_clear(isometry_t isom);

/************************************************************************
 *
 * Function: isometry_is_equal
 *
 * Description: Checks whether two isometries are equal.
 *
 * Arguments:
 *     + isom1 (isometry_t) - an isometry
 *     + isom2 (isometry_t) - another isometry
 *
 * Returns:
 *     + (bool) - true if isom1 == isom2, false otherwise
 *
 ************************************************************************
 */

bool isometry_is_equal(const isometry_t isom1, const isometry_t isom2);

// This function was for checking if a rational matrix is an isometry between
// two lattices, where the gram matrices are given as matrix_TYP* (carat type).
// However, it is no onger needed.
// bool is_isometry(matrix_TYP* s, matrix_TYP* q1, matrix_TYP* q2, int denom);

/************************************************************************
 *
 * Function: is_isometric
 *
 * Description: Checks whether two integral lattices are isometric.
 *              If so, finds an isometry between them.
 *
 * Arguments:
 *     + q1 (const square_matrix_t) - gram matrix of a lattice L1
 *     + q2 (const square_matrix_t) - gram matrix of a lattice L2
 *
 * Returns:
 *     + (bool) - true if L1 is isometric to L2
 *     + s (isometry_t) - an isometry s such that s(L1) = L2,
 *                        i.e. s^t * q1 * s = q2
 *                        If L1 is not isometric to L2, s is unchanged
 *
 ************************************************************************
 */

bool is_isometric(isometry_t s, const square_matrix_t q1, const square_matrix_t q2);

// What follows are basic matrix operations.
// We have special functions for isometries, since we want to update
// the inverse isometry at the same time.
// This saves computation of the inverse matrices eventually.

/*************************************************************************
 *
 * Function: isometry_swap_vecs
 *
 * Description: Swaps two basis vectors in an isometry.
 *
 * Arguments:
 *     + isom (isometry_t) - an isometry
 *     + col1, col2 (int) - the indices of the two basis vectors to swap
 *
 * Returns:
 *     + isom (isometry_t) - the resulting isometry
 *                           (with col1, col2 swapped)
 *
 *************************************************************************
 */

// basis vectors are columns (?check?)
void isometry_swap_vecs(isometry_t isom, int col1, int col2);

/*************************************************************************
 *
 * Function: isometry_add_vec
 *
 * Description: Adds a multiple of one basis vector to another.
 *
 * Arguments:
 *     + isom (isometry_t) - an isometry
 *     + target_col (int) - the index of basis vector (column)
 *                          to be updated
 *     + c (Z64) - the scalar to multiply by
 *     + src_col (int) - the index of the basis vector (column)
 *                       whose multiple we are adding
 *
 * Returns:
 *     + isom (isometry_t) - the resulting isometry
 *                           (isom[target_col] += c * isom[src_col])
 *
 *************************************************************************
 */

void isometry_add_vec(isometry_t isom, int target_col, Z64 c, int src_col);

/*************************************************************************
 *
 * Function: isometry_replace_vec
 *
 * Description: Replaces a basis vector by a given vector.
 *              This operation is useful when constructing neighbors,
 *              since we want to set a basis vector to be an isotropic
 *              vector.
 *
 * Arguments:
 *     + isom (isometry_t) - an isometry
 *     + col (int) - the index of basis vector (column)
 *                   to be updated
 *     + x (const vector_t) - the new column (as a vector)
 *
 * Returns:
 *     + isom (isometry_t) - the resulting isometry
 *                           (isom[col] = x)
 *
 *************************************************************************
 */

void isometry_replace_vec(isometry_t isom, int col, const vector_t x);

/*************************************************************************
 *
 * Function: isometry_vec_scalar_div
 *
 * Description: Divides a basis vector by a scalar.
 *
 * Arguments:
 *     + isom (isometry_t) - an isometry
 *     + col (int) - the index of basis vector (column)
 *                   to be updated
 *     + scalar (int) - the scalar to divide by
 *
 * Returns:
 *     + isom (isometry_t) - the resulting isometry
 *                           (isom[col] /= scalar)
 *
 *************************************************************************
 */

void isometry_vec_scalar_div(isometry_t isom, int col, int scalar);

/*************************************************************************
 *
 * Function: isometry_vec_scalar_mul
 *
 * Description: Multiplies a basis vector by a scalar.
 *
 * Arguments:
 *     + isom (isometry_t) - an isometry
 *     + col (int) - the index of basis vector (column)
 *                   to be updated
 *     + scalar (int) - the scalar to multiply by
 *
 * Returns:
 *     + isom (isometry_t) - the resulting isometry
 *                           (isom[col] *= scalar)
 *
 *************************************************************************
 */

void isometry_vec_scalar_mul(isometry_t isom, int col, int scalar);

#endif // __ISOMETRY_H__
