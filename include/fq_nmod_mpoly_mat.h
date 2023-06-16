/*************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : fq_nmod_mpoly_mat.h
 *
 * Description: Matrices of multivariable polynomials over
 *              finite fields.
 *
 *************************************************************
 */

#ifndef __FQ_NMOD_MPOLY_MAT_H__
#define __FQ_NMOD_MPOLY_MAT_H__

#ifdef FQ_NMOD_MPOLY_MAT_INLINES_C
#define FQ_NMOD_MPOLY_MAT_INLINE FLINT_DLL
#else
#define FQ_NMOD_MPOLY_MAT_INLINE static __inline__
#endif

// Required packages dependencies

#include <flint/flint.h>
#include <flint/fq_nmod_mpoly.h>

/****************************************************************
 *
 * Type: fq_nmod_mpoly_mat_t
 *
 * Description: matrix whose entries are multivariable polys. 
 *
 * Fields:
 *     + entries (fq_nmod_mpoly_struct*) - a list of the 
 *                     entries, multivariable polynomials
 *     + r (slong) - number of rows
 *     + c (slong) - number of columns
 *     + rows (fq_nmod_mpoly_struct**) - a list of pointers to
 *                                       the rows of the matrix 
 *
 ****************************************************************
 */

typedef struct
{
  fq_nmod_mpoly_struct* entries;
  slong r;
  slong c;
  fq_nmod_mpoly_struct** rows;
} fq_nmod_mpoly_mat_struct;

typedef fq_nmod_mpoly_mat_struct fq_nmod_mpoly_mat_t[1];

/******************************************************************
 *
 * Function: fq_nmod_mpoly_mat_init
 *
 * Description: Allocates memory for a matrix.
 *
 * Arguments:
 *     + rows (slong) - number of rows
 *     + cols (slong) - number of columns
 *     + R (const fq_nmod_mpoly_ctx_t) - the relevant ring of
 *                                       multivariable polynomials
 *
 * Returns:
 *     + mat (fq_nmod_mpoly_mat_t) - the empty matrix
 *
 ******************************************************************
 */

void fq_nmod_mpoly_mat_init(fq_nmod_mpoly_mat_t mat, slong rows, slong cols, const fq_nmod_mpoly_ctx_t R);

/******************************************************************
 *
 * Function: fq_nmod_mpoly_mat_clear
 *
 * Description: Deallocates memory for a matrix.
 *
 * Arguments:
 *     + mat (fq_nmod_mpoly_mat_t) - the matrix
 *     + R (const fq_nmod_mpoly_ctx_t) - the ring of multivariable
 *                                       polynomials
 *
 ******************************************************************
 */

void fq_nmod_mpoly_mat_clear(fq_nmod_mpoly_mat_t mat, const fq_nmod_mpoly_ctx_t R);

/******************************************************************
 *
 * Function: fq_nmod_mpoly_mat_print
 *
 * Description: Prints a matrix.
 *
 * Arguments:
 *     + mat (const fq_nmod_mpoly_mat_t) - the matrix
 *     + var_names (const char**) - a list of strings for
 *                                  the variable names
 *     + R (const fq_nmod_mpoly_ctx_t) - the ring of multivariable
 *                                       polynomials
 *
 ******************************************************************
 */

void fq_nmod_mpoly_mat_print(const fq_nmod_mpoly_mat_t mat, const char** var_names, const fq_nmod_mpoly_ctx_t R);

/******************************************************************
 *
 * Function: fq_nmod_mpoly_mat_entry
 *
 * Description: Returns an entry in the matrix.
 *
 * Arguments:
 *     + mat (const fq_nmod_mpoly_mat_t) - the matrix
 *     + i (slong) - the index of the row
 *     + j (slong) - the index of the column
 *
 * Returns:
 *     + (fq_nmod_mpoly_struct*) - the polynomial mat[i,j] 
 *
 ******************************************************************
 */

FQ_NMOD_MPOLY_MAT_INLINE
fq_nmod_mpoly_struct * fq_nmod_mpoly_mat_entry(const fq_nmod_mpoly_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

/******************************************************************
 *
 * Function: fq_nmod_mpoly_mat_nrows
 *
 * Description: Returns the number of rows in the matrix.
 *
 * Arguments:
 *     + mat (const fq_nmod_mpoly_mat_t) - the matrix
 *
 * Returns:
 *     + (slong) - the number of rows in mat
 *
 ******************************************************************
 */

FQ_NMOD_MPOLY_MAT_INLINE
slong fq_nmod_mpoly_mat_nrows(const fq_nmod_mpoly_mat_t mat)
{
   return mat->r;
}

/******************************************************************
 *
 * Function: fq_nmod_mpoly_mat_ncols
 *
 * Description: Returns the number of columns in the matrix.
 *
 * Arguments:
 *     + mat (const fq_nmod_mpoly_mat_t) - the matrix
 *
 * Returns:
 *     + (slong) - the number of columns in mat
 *
 ******************************************************************
 */

FQ_NMOD_MPOLY_MAT_INLINE
slong fq_nmod_mpoly_mat_ncols(const fq_nmod_mpoly_mat_t mat)
{
   return mat->c;
}

#endif // __FQ_NMOD_MPOLY_MAT_H__
