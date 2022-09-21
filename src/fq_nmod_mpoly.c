#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_quadratic_part(fq_nmod_mpoly_t quad, const fq_nmod_mpoly_t f, const fq_nmod_mpoly_ctx_t ctx)
{
  slong i, j, N;
  ulong deg;
  fq_nmod_t coeff;
  ulong* exp;

  // !! TODO - maube use get_term and get_exp instead
  fq_nmod_init(coeff, ctx->fqctx);
  N = mpoly_words_per_exp(f->bits, ctx->minfo);
  fq_nmod_mpoly_zero(quad, ctx);
  for (i = 0; i < f->length; i++) {
    deg = 0;
    exp = f->exps + i*N;
    for (j = 0; j < N; j++)
      deg += exp[j];
    if (deg == 2) {
      fq_nmod_mpoly_get_coeff_fq_nmod_ui(coeff, f, exp, ctx);
      fq_nmod_mpoly_set_coeff_fq_nmod_ui(quad, coeff, exp, ctx);
    }
  }
  fq_nmod_clear(coeff, ctx->fqctx);
  return;
}

void fq_nmod_mpoly_linear_part(fq_nmod_mat_t lin, const fq_nmod_mpoly_t f, const fq_nmod_mpoly_ctx_t ctx)
{
  slong i, nvars;
  ulong* exp;

  nvars = fq_nmod_mpoly_ctx_nvars(ctx);
  
  exp = (ulong*)malloc(nvars * sizeof(ulong));
  for (i = 0; i < nvars; i++)
    exp[i] = 0;

  for (i = 0; i < nvars; i++) {
    exp[i] = 1;
    if (i >= 1) exp[i-1] = 0;
    fq_nmod_mpoly_get_coeff_fq_nmod_ui(fq_nmod_mat_entry(lin, 0, i), f, exp, ctx);
  }

  free(exp);
  return;
}
