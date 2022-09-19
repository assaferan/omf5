#ifndef __FQ_NMOD_QUAD_H__
#define __FQ_NMOD_QUAD_H__

#include "flint/fq_nmod.h"
#include "flint/fq_nmod_mat.h"

void fq_nmod_quad_evaluate(fq_nmod_t value, const fq_nmod_mat_t q, const fq_nmod_mat_t vec, const fq_nmod_ctx_t F);

void fq_nmod_quad_transform(fq_nmod_mat_t new_gram, const fq_nmod_mat_t gram, const fq_nmod_mat_t isom, const fq_nmod_ctx_t F);

bool fq_nmod_quad_isotropic_vector(fq_nmod_mat_t vec, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, slong start, bool deterministic);

void fq_nmod_quad_split_hyperbolic_plane(const fq_nmod_mat_t vec, fq_nmod_mat_t gram, fq_nmod_mat_t basis, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, slong start);

void fq_nmod_quad_hyperbolize(fq_nmod_mat_t gram, fq_nmod_mat_t basis, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, bool deterministic, slong start);

void fq_nmod_quad_decompose(fq_nmod_mat_t gram, fq_nmod_mat_t basis, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, bool deterministic);

void fq_nmod_poly_set_fq_nmod_quad(fq_nmod_poly_t poly, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, const fq_nmod_mpoly_ctx_t R);

#endif // __FQ_NMOD_QUAD_H__
