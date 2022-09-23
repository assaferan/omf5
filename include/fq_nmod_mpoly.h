#ifndef __FQ_NMOD_MPOLY_H__
#define __FQ_NMOD_MPOLY_H__

#include "flint/fq_nmod_mpoly.h"

void fq_nmod_mpoly_quadratic_part(fq_nmod_mpoly_t quad, const fq_nmod_mpoly_t f, const fq_nmod_mpoly_ctx_t ctx);
void fq_nmod_mpoly_linear_part(fq_nmod_mat_t lin, const fq_nmod_mpoly_t f, const fq_nmod_mpoly_ctx_t ctx);

#endif // __FQ_NMOD_MPOLY_H__
