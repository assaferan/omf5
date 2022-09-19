#ifndef __FQ_NMOD_MAT_H__
#define __FQ_NMOD_MAT_H__

#include "flint/fq_nmod.h"
#include "flint/fq_nmod_mat.h"

bool fq_nmod_mat_is_zero_row(const fq_nmod_mat_t mat, slong row, const fq_nmod_ctx_t F);
void fq_nmod_mat_swap_rows(fq_nmod_mat_t mat, slong row1, slong row2, const fq_nmod_ctx_t F);
void fq_nmod_mat_swap_cols(fq_nmod_mat_t mat, slong col1, slong col2, const fq_nmod_ctx_t F);
void fq_nmod_mat_add_row(fq_nmod_mat_t mat, slong dst_row, slong src_row, const fq_nmod_t scalar, const fq_nmod_ctx_t F);
void fq_nmod_mat_add_col(fq_nmod_mat_t mat, slong dst_col, slong src_col, const fq_nmod_t scalar, const fq_nmod_ctx_t F);
void fq_nmod_mat_mul_row(fq_nmod_mat_t mat, slong row, const fq_nmod_t scalar, const fq_nmod_ctx_t F);
void fq_nmod_mat_mul_col(fq_nmod_mat_t mat, slong col, const fq_nmod_t scalar, const fq_nmod_ctx_t F);

#endif // __FQ_NMOD_MAT_H__
