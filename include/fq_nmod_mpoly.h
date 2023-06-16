/*************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : fq_nmod_mpoly.h
 *
 * Description: functions for handling the type fq_nmod_mpoly
 *              from FLINT, multivariable polynomials
 *              over finite fields.
 *
 *************************************************************
 */

#ifndef __FQ_NMOD_MPOLY_H__
#define __FQ_NMOD_MPOLY_H__

// Required packages dependencies

#include <flint/fq_nmod_mpoly.h>

/*******************************************************************
 *
 * Function: fq_nmod_mpoly_quadratic_part
 *
 * Description: Get the degree 2 homogenous component of
 *              a polynomial f.
 *
 * Arguments:
 *     + f (const fq_nmod_mpoly_t) - the polynomial
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + quad (fq_nmod_mpoly_t) - the quadratic part of f
 *
 ******************************************************************
 */

void fq_nmod_mpoly_quadratic_part(fq_nmod_mpoly_t quad, const fq_nmod_mpoly_t f, const fq_nmod_mpoly_ctx_t ctx);

/*******************************************************************
 *
 * Function: fq_nmod_mpoly_linear_part
 *
 * Description: Get the degree 1 homogenous component of
 *              a polynomial f.
 *
 * Arguments:
 *     + f (const fq_nmod_mpoly_t) - the polynomial
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + lin (fq_nmod_mat_t) - a row vector representing the linear
 *                             component of f. 
 *
 ******************************************************************
 */

void fq_nmod_mpoly_linear_part(fq_nmod_mat_t lin, const fq_nmod_mpoly_t f, const fq_nmod_mpoly_ctx_t ctx);

#endif // __FQ_NMOD_MPOLY_H__
