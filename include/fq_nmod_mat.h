/***********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : fq_nmod_mat.h
 *
 * Description: functions for handling the type fq_nmod_mat
 *              from FLINT, matrices over finite fields.
 *
 ***********************************************************
 */

#ifndef __FQ_NMOD_MAT_H__
#define __FQ_NMOD_MAT_H__

// Required packages dependencies

#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>

// Self dependencies

#include "typedefs.h"

/*******************************************************************
 *
 * Function: fq_nmod_mat_is_zero_row
 *
 * Description: Check if a certain row is all zeros.
 *
 * Arguments:
 *     + mat (const fq_nmod_mat_t) - the matrix
 *     + row (slong) - the index of the row
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + (bool) - true if the row is a zero row
 *
 ******************************************************************
 */

bool fq_nmod_mat_is_zero_row(const fq_nmod_mat_t mat, slong row, const fq_nmod_ctx_t F);

// These basic functions are no longer needed as they have analogues implemented in FLINT
// void fq_nmod_mat_swap_rows(fq_nmod_mat_t mat, slong row1, slong row2, const fq_nmod_ctx_t F);
// void fq_nmod_mat_swap_cols(fq_nmod_mat_t mat, slong col1, slong col2, const fq_nmod_ctx_t F);

/*******************************************************************
 *
 * Function: fq_nmod_mat_add_row
 *
 * Description: add a multiple of a row to another row,
 *              (R_i <- R_i + scalar * R_j).
 *
 * Arguments:
 *     + mat (fq_nmod_mat_t) - the matrix
 *     + dst_row (slong) - the index of the destination row (i)
 *     + src_row (slong) - the index of the source row (j)
 *     + scalar (const fq_nmod_t) - the scalar we multiply by
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + mat (fq_nmod_mat_t) - the matrix
 *
 ******************************************************************
 */

void fq_nmod_mat_add_row(fq_nmod_mat_t mat, slong dst_row, slong src_row, const fq_nmod_t scalar, const fq_nmod_ctx_t F);

/*******************************************************************
 *
 * Function: fq_nmod_mat_add_col
 *
 * Description: add a multiple of a column to another column,
 *              (C_i <- C_i + scalar * C_j).
 *
 * Arguments:
 *     + mat (fq_nmod_mat_t) - the matrix
 *     + dst_col (slong) - the index of the destination column (i)
 *     + src_col (slong) - the index of the source column (j)
 *     + scalar (const fq_nmod_t) - the scalar we multiply by
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + mat (fq_nmod_mat_t) - the matrix
 *
 ******************************************************************
 */

void fq_nmod_mat_add_col(fq_nmod_mat_t mat, slong dst_col, slong src_col, const fq_nmod_t scalar, const fq_nmod_ctx_t F);

/*******************************************************************
 *
 * Function: fq_nmod_mat_mul_row
 *
 * Description: multiply a row by a scalar,
 *              (R_i <- scalar * R_i).
 *
 * Arguments:
 *     + mat (fq_nmod_mat_t) - the matrix
 *     + row (slong) - the index of the row (i)
 *     + scalar (const fq_nmod_t) - the scalar we multiply by
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + mat (fq_nmod_mat_t) - the matrix
 *
 ******************************************************************
 */

void fq_nmod_mat_mul_row(fq_nmod_mat_t mat, slong row, const fq_nmod_t scalar, const fq_nmod_ctx_t F);

/*******************************************************************
 *
 * Function: fq_nmod_mat_mul_col
 *
 * Description: multiply a column by a scalar,
 *              (C_i <- scalar * C_i).
 *
 * Arguments:
 *     + mat (fq_nmod_mat_t) - the matrix
 *     + col (slong) - the index of the column (i)
 *     + scalar (const fq_nmod_t) - the scalar we multiply by
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + mat (fq_nmod_mat_t) - the matrix
 *
 ******************************************************************
 */

void fq_nmod_mat_mul_col(fq_nmod_mat_t mat, slong col, const fq_nmod_t scalar, const fq_nmod_ctx_t F);

/*******************************************************************
 *
 * Function: fq_nmod_mat_init_set_fmpz_mat
 *
 * Description: Reduce an integer matrix modulo a prime p.
 *              Assumes F = F_p is a prime field.
 *
 * Arguments:
 *     + mat (const fmpz_mat_t) - the integral matrix to be reduced
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + dest (fq_nmod_mat_t) - the reduced matrix
 *
 ******************************************************************
 */

void fq_nmod_mat_init_set_fmpz_mat(fq_nmod_mat_t dest, const fmpz_mat_t mat, const fq_nmod_ctx_t F);

/*******************************************************************
 *
 * Function: fq_nmod_mat_transpose
 *
 * Description: Transpose a matrix.
 *
 * Arguments:
 *     + mat (const fq_nmod_mat_t) - the matrix to be transposed
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + mat_t (fq_nmod_mat_t) - the transposed matrix
 *
 ******************************************************************
 */

void fq_nmod_mat_transpose(fq_nmod_mat_t mat_t, const fq_nmod_mat_t mat, const fq_nmod_ctx_t F);

/*******************************************************************
 *
 * Function: fq_nmod_mat_rref_trans
 *
 * Description: Bring to a row reduce echelon form, together with
 *              the transformation matrix.
 *
 * Arguments:
 *     + mat (fq_nmod_mat_t) - the matrix to be echelonized, M
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + mat (fq_nmod_mat_t) - the echelonized matrix, E
 *     + trans (fq_nmod_mat_t) - the transformation marix T, such
 *                               that T*E = M
 *
 ******************************************************************
 */

void fq_nmod_mat_rref_trans(fq_nmod_mat_t mat, fq_nmod_mat_t trans, const fq_nmod_ctx_t F);

/*******************************************************************
 *
 * Function: fq_nmod_mat_kernel
 *
 * Description: Find the kernel of a matrix.
 *
 * Arguments:
 *     + mat (fq_nmod_mat_t) - the matrix
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + ker (fq_nmod_mat_t) - a matrix whose rows form a basis for
 *                             the null space of the matrix mat.
 *                             its size is initialized inside the
 *                             function.
 *
 ******************************************************************
 */

void fq_nmod_mat_kernel(fq_nmod_mat_t ker, const fq_nmod_mat_t mat, const fq_nmod_ctx_t F);

#endif // __FQ_NMOD_MAT_H__
