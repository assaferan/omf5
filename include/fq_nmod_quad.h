/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : fq_nmod_quad.h
 *
 * Description: functions for handling quadratic forms
 *              over finite fields.
 *
 **********************************************************
 */

#ifndef __FQ_NMOD_QUAD_H__
#define __FQ_NMOD_QUAD_H__

// Required packages dependencies

#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>
#include <flint/fq_nmod_mpoly.h>

// Self dependencies

#include "typedefs.h"

/*************************************************************************
 *
 * Function: fq_nmod_quad_evaluate
 *
 * Description: Evaluate a quadratic form at a vector.
 *
 * Arguments:
 *     + q (const fq_nmod_mat_t) - the gram matrix of the quadratic form
 *     + v (const fq_nmod_mat_t) - the vector (as a matrix with one row)
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + value (fq_nmod_t) - the value q(v)
 *
 *************************************************************************
 */

void fq_nmod_quad_evaluate(fq_nmod_t value, const fq_nmod_mat_t q, const fq_nmod_mat_t v, const fq_nmod_ctx_t F);

/****************************************************************************
 *
 * Function: fq_nmod_quad_transform
 *
 * Description: Transform a gram matrix under a change of basis (isometry)
 *
 * Arguments:
 *     + gram (const fq_nmod_mat_t) - the gram matrix of the quadratic form
 *     + isom (const fq_nmod_mat_t) - an isometry (its rows are the basis)
 *     + F (const fq_nmod_ctx_t) - the finite field
 *
 * Returns:
 *     + new_gram (fq_nmod_mat_t) - the new gram matrix isom * gram * isom^T
 *
 ****************************************************************************
 */

void fq_nmod_quad_transform(fq_nmod_mat_t new_gram, const fq_nmod_mat_t gram, const fq_nmod_mat_t isom, const fq_nmod_ctx_t F);

/****************************************************************************
 *
 * Function: fq_nmod_quad_isotropic_vector
 *
 * Description: Find an isotropic vector for a quadratic form.
 *
 * Arguments:
 *     + q (const fq_nmod_mat_t) - the gram matrix of the quadratic form
 *     + F (const fq_nmod_ctx_t) - the finite field
 *     + start (slong) - consider only the rows and columns from "start",
 *                       i.e. find an isotropic vector for q[-start:,-start:]
 *     + deterministic (bool) - true to avoid using probabilistic means,
 *                              useful for testing and debugging.
 *
 * Returns:
 *     + vec (fq_nmod_mat_t) - an isotropic vector for q, precisely v 
 *                             such that q[-start:-start:](v[-start:]) = 0
 *     + (bool) - true if succeeded, false if the form is anisotropic
 *
 ****************************************************************************
 */


bool fq_nmod_quad_isotropic_vector(fq_nmod_mat_t vec, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, slong start, bool deterministic);

/******************************************************************************
 *
 * Function: fq_nmod_quad_split_hyperbolic_plane
 *
 * Description: Split a hyperbolic plane from a quadratic form.
 *
 * Arguments:
 *     + v (const fq_nmod_mat_t) - an isotropic vector (as a row matrix)
 *     + q (const fq_nmod_mat_t) - the gram matrix of the quadratic form
 *     + F (const fq_nmod_ctx_t) - the finite field
 *     + start (slong) - consider only the rows and columns from "start",
 *                       i.e. split a hyperbolic plane from q[-start:,-start:]
 *
 * Returns:
 *     + basis (fq_nmod_mat_t) - a basis, two of whose elements form a
 *                               hyperbolic plane containing v
 *     + gram (fq_nmod_mat_t) - the gram matrix with respect to basis 
 *
 ******************************************************************************
 */

void fq_nmod_quad_split_hyperbolic_plane(const fq_nmod_mat_t v, fq_nmod_mat_t gram, fq_nmod_mat_t basis, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, slong start);

/******************************************************************************
 *
 * Function: fq_nmod_quad_hyperbolize
 *
 * Description: Split all hyperbolic planes from a quadratic form.
 *
 * Arguments:
 *     + q (const fq_nmod_mat_t) - the gram matrix of the quadratic form
 *     + F (const fq_nmod_ctx_t) - the finite field
 *     + deterministic (bool) - true to avoid using probabilistic means,
 *                              useful for testing and debugging.
 *     + start (slong) - consider only the rows and columns from "start",
 *                       i.e. split hyperbolic planes from q[-start:,-start:]
 *
 * Returns:
 *     + basis (fq_nmod_mat_t) - a basis, 2h of whose elements form h
 *                               hyperbolic planes, where h is the Witt index
 *     + gram (fq_nmod_mat_t) - the gram matrix with respect to basis 
 *
 ******************************************************************************
 */

void fq_nmod_quad_hyperbolize(fq_nmod_mat_t gram, fq_nmod_mat_t basis, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, bool deterministic, slong start);

/******************************************************************************
 *
 * Function: fq_nmod_quad_decompose
 *
 * Description: Decompose a quadratic form as q = h + a + r
 *
 * Arguments:
 *     + q (const fq_nmod_mat_t) - the gram matrix of the quadratic form
 *     + F (const fq_nmod_ctx_t) - the finite field
 *     + deterministic (bool) - true to avoid using probabilistic means,
 *                              useful for testing and debugging.
 *
 * Returns:
 *     + basis (fq_nmod_mat_t) - a basis, whose first 2h elements form h
 *                               hyperbolic planes, where h is the Witt index,
 *                               its next a elements are a basis for the
 *                               anisotropic subspace, and its last r elements
 *                               form a basis for the radical.
 *     + gram (fq_nmod_mat_t) - the gram matrix with respect to basis 
 *
 ******************************************************************************
 */

void fq_nmod_quad_decompose(fq_nmod_mat_t gram, fq_nmod_mat_t basis, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, bool deterministic);

/******************************************************************************
 *
 * Function: fq_nmod_poly_set_fq_nmod_quad
 *
 * Description: Construct the polynomial representing the quadratic form.
 *
 * Arguments:
 *     + q (const fq_nmod_mat_t) - the gram matrix of the quadratic form
 *     + F (const fq_nmod_ctx_t) - the finite field
 *     + R (const fq_nmod_mpoly_ctx_t) - the ring of polynomials
 *
 * Returns:
 *     + poly (fq_nmod_mpoly_t) - the polynomial q(x),
 *                                where x = x_1,...,x_n are the variables of R
 *
 ******************************************************************************
 */


void fq_nmod_poly_set_fq_nmod_quad(fq_nmod_mpoly_t poly, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, const fq_nmod_mpoly_ctx_t R);

#endif // __FQ_NMOD_QUAD_H__
