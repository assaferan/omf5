/*****************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : mass.h
 *
 * Description: Functions for computing mass formula
 *
 *****************************************************************
 */

#ifndef __MASS_H__
#define __MASS_H__

// Required packages dependencies

#include <carat/matrix.h>

#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

// Self dependencies

#include "jordan.h"
#include "square_matrix.h"

/***************************************************************************
 *
 * Function: local_factor
 *
 * Description: Given a gram matrix of a lattice in Q \cap Z_p, compute the
 *              factor it contributes to the mass formula.             
 *
 * Arguments:
 *     + q (const fmpq_mat_t) - the gram matrix
 *     + p (const fmpz_t) - the prime p
 *
 * Returns:
 *     + f (fmpq_t) - the local factor
 *
 ***************************************************************************
 */

void local_factor(fmpq_t f, const fmpq_mat_t g, const fmpz_t p);

/***************************************************************************
 *
 * Function: diagonal_join
 *
 * Description: Given a Jordan decomposition of a lattice,
 *              returns a rational matrix which is the direct sum of
 *              all gram matrices in the decomposition.
 *
 * Arguments:
 *     + jordan (const jordan_data_t) - the Jordan decomposition
 *
 * Returns:
 *     + joined (fmpq_mat_t) - the direct sum of the matrices
 *
 ***************************************************************************
 */

void diagonal_join(fmpq_mat_t joined, const jordan_data_t jordan);

/***************************************************************************
 *
 * Function: get_mass
 *
 * Description: Given the gram matrix of a lattice, computes the mass
 *              of its genus.
 *
 * Arguments:
 *     + q (const square_matrix_t) - the gram matrix of L
 *
 * Returns:
 *     + mass (fmpq_t) - Mass(L) = sum(aut(L')^(-1)) over L'~L
 *
 ***************************************************************************
 */

/* compute the mass of the genus represented by q */
void get_mass(fmpq_t mass, const square_matrix_t q);

#endif // __MASS_H__
