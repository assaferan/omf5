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
