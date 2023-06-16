/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : fmpq_poly.h
 *
 * Description: functions for handling the type fmpq_poly
 *              from FLINT, polynomials over rationals.
 *
 **********************************************************
 */

#ifndef __FMPQ_POLY_H__
#define __FMPQ_POLY_H__

// Required packages dependencies

#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>

// Self dependencies

#include "typedefs.h"

/*************************************************************
 *
 * Type: fmpq_poly_factor_t
 *
 * Description: polynomial factorization data.
 *
 * Fields:
 *     + p (fmpq_poly_struct*) - a list of factors
 *     + exp (slong*)          - a list of exponents
 *     + num (slong)           - length of both lists
 *
 * Used to hold the factorization data of the polynomial
 * prod p[i]^exp[i]. Based on a similar structure for
 * fmpz_poly.
 *
 ************************************************************
 */

typedef struct
{
  fmpq_poly_struct *p;
  slong *exp;
  slong num;
}
fmpq_poly_factor_struct;

typedef fmpq_poly_factor_struct fmpq_poly_factor_t[1];

/***************************************************************
 *
 * Function: fmpq_poly_factor
 *
 * Description: Factors a polynomial over the rationals.
 *
 * Arguments:
 *     + f (const fmpq_poly_t) - the polynomial to be factored
 *
 * Returns:
 *     + factors (fmpq_poly_factor_t) - the factorization data
 *
 **************************************************************
 */

void fmpq_poly_factor(fmpq_poly_factor_t factors, const fmpq_poly_t f);

/*******************************************************************
 *
 * Function: fmpq_poly_factor_clear
 *
 * Description: Clears the memory allocated for factorization data.
 *
 * Arguments:
 *     + factors (fmpq_poly_factor_t) - the factorization data
 *
 ******************************************************************
 */

void fmpq_poly_factor_clear(fmpq_poly_factor_t factors);

/*******************************************************************
 *
 * Function: fmpq_poly_is_irreducible
 *
 * Description: Checks irreduciblity of a polynomial.
 *
 * Arguments:
 *     + f (const fmpq_poly_t) - the polynomial
 *
 * Returns:
 *     + bool - true if f is irreducible
 *
 ******************************************************************
 */

bool fmpq_poly_is_irreducible(const fmpq_poly_t f);

/*******************************************************************
 *
 * Function: fmpq_poly_evaluate_fmpq_mat
 *
 * Description: Evaluate a polynomial at a matrix.
 *
 * Arguments:
 *     + f (const fmpq_poly_t) - the polynomial
 *     + A (const fmpa_mat_t) - the matrix
 *
 * Returns:
 *     + res (fmpq_mat_t) - the matrix f(A)
 *
 ******************************************************************
 */

void fmpq_poly_evaluate_fmpq_mat(fmpq_mat_t res, const fmpq_poly_t f, const fmpq_mat_t A);

#endif // __FMPQ_POLY_H__
