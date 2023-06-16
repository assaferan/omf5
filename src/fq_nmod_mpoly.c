/*************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : fq_nmod_mpoly.c
 *
 * Description: functions for handling the type fq_nmod_mpoly
 *              from FLINT, multivariable polynomials
 *              over finite fields.
 *
 *************************************************************
 */

// Self dependencies

#include "fq_nmod_mpoly.h"

// Get the degree 2 homogeneous part of the polynomial
void fq_nmod_mpoly_quadratic_part(fq_nmod_mpoly_t quad, const fq_nmod_mpoly_t f, const fq_nmod_mpoly_ctx_t ctx)
{
  slong i, j;
  ulong deg;
  fq_nmod_t coeff;
  ulong* exp;

  exp = (ulong*)malloc((ctx->minfo->nvars)*sizeof(ulong));
  
  // !! TODO - maube use get_term and get_exp instead
  fq_nmod_init(coeff, ctx->fqctx);
  fq_nmod_mpoly_zero(quad, ctx);

  // we just loop over the monomials
  for (i = 0; i < f->length; i++) {
    fq_nmod_mpoly_get_term_exp_ui(exp, f, i, ctx);
    deg = 0;
    for (j = 0; j < ctx->minfo->nvars; j++)
      deg += exp[j];
    // if the degree of the monomial is 2, we add it (multiplied by its coefficient) to the result
    if (deg == 2) {
      fq_nmod_mpoly_get_coeff_fq_nmod_ui(coeff, f, exp, ctx);
      fq_nmod_mpoly_set_coeff_fq_nmod_ui(quad, coeff, exp, ctx);
    }
  }
  fq_nmod_clear(coeff, ctx->fqctx);
  free(exp);
  return;
}

// Get the degree 1 homogeneous part of the polynomial, as a row vector
void fq_nmod_mpoly_linear_part(fq_nmod_mat_t lin, const fq_nmod_mpoly_t f, const fq_nmod_mpoly_ctx_t ctx)
{
  slong i, nvars;
  ulong* exp;

  nvars = fq_nmod_mpoly_ctx_nvars(ctx);
  
  exp = (ulong*)malloc(nvars * sizeof(ulong));
  for (i = 0; i < nvars; i++)
    exp[i] = 0;

  for (i = 0; i < nvars; i++) {
    // set the exponent to be e_i, so we get the coefficient of x_i
    exp[i] = 1;
    if (i >= 1) exp[i-1] = 0;
    fq_nmod_mpoly_get_coeff_fq_nmod_ui(fq_nmod_mat_entry(lin, 0, i), f, exp, ctx);
  }

  free(exp);
  return;
}
