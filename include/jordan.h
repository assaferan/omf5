/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : jordan.h
 *
 * Description: Data structure for Jordan decomposition
 *              of a quadratic form over Z_p.
 *
 **********************************************************
 */

#ifndef __JORDAN_H__
#define __JORDAN_H__

// Required packages dependencies

#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>

// self dependencies

#include "square_matrix.h"

/***************************************************************************
 *
 * Type: jordan_data_t
 *
 * Description: this struct holds data for a Jordan decomposition
 *              of a quadratic form over the ring of p-adic numbers Z_p.
 *
 * Fields:
 *     + matrices (fmpq_mat_t*) - an array of rational matrices,
 *                                representing the basis vectors
 *                                of each of the summands in the
 *                                Jordan decomposition
 *     + grams (fmpq_mat_t*) - an array of rational matrices, which are
 *                             the gram matrices of each of the subspaces,
 *                             scaled to be an atomic gram matrix
 *     + exponents (ulong*) - the exponents of the scaling factors
 *     + num_blocks (ulong) - the number of blocks in the decomposition
 *
 ***************************************************************************
 */

typedef struct {
  fmpq_mat_t* matrices;
  fmpq_mat_t* grams;
  ulong* exponents;
  ulong num_blocks;
} jordan_data;

typedef jordan_data jordan_data_t[1];

/**********************************************************************
 *
 * Function: jordan_data_init
 *
 * Description: Allocates memory for a jordan_data_t,
 *              with enough room for a given number of blocks.
 *
 * Arguments:
 *     + n (ulong) - maximal number of blocks in the decomposition
 *
 * Returns:
 *     + jordan (jordan_data_t) - an empty Jordan decomposition,
 *                                with enough room to store n blocks
 *
 **********************************************************************
 */

void jordan_data_init(jordan_data_t jordan, ulong n);

/**********************************************************************
 *
 * Function: jordan_data_clear
 *
 * Description: Clears memory allocated for a jordan_data_t.
 *
 * Arguments:
 *     + jordan (jordan_data_t) - the jordan_data_t to be cleared
 *
 **********************************************************************
 */

void jordan_data_clear(jordan_data_t jordan);

/**********************************************************************
 *
 * Function: inner_product
 *
 * Description: Computes the inner product between two rows of a
 *              matrix with rational entries, with respect to
 *              a given bilinear form with integral coefficients.
 *
 * Arguments:
 *     + G (const fmpz_mat_t) - the bilinear form (gram matrix)
 *     + S (const fmpq_mat_t) - the rational matrix whose rows we
 *                              wish to compute inner product for
 *     + idx1, idx2 (ulong) - indices of the rows whose inner product
 *                            we wish to compute
 *
 * Returns:
 *     + res (fmpq_t) - the inner product of S[idx1] and S[idx2]
 *                      with respect to G, i.e. S[idx1]*G*S[idx2]^T
 *
 **********************************************************************
 */

void inner_product(fmpq_t res, const fmpz_mat_t G,
		   const fmpq_mat_t S, ulong idx1, ulong idx2);

/**********************************************************************
 *
 * Function: jordan_decomposition
 *
 * Description: Computes the Jordan decomposition of a quadratic form
 *              with integral gram matrix over Z_p.
 *
 * Arguments:
 *     + q (const square_matrix_t) - the gram matrix of the quadratic
 *                                   form
 *     + p (const fmpz_t) - the prime p
 *
 * Returns:
 *     + jordan (jordan_data_t) - the Jordan decomposition
 *
 **********************************************************************
 */

void jordan_decomposition(jordan_data_t jordan, const square_matrix_t q, const fmpz_t p);

#endif // __JORDAN_H__
