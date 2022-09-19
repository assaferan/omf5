#ifndef __FQ_NMOD_POLY_MAT_H__
#define __FQ_NMOD_POLY_MAT_H__

#include "flint/fq_nmod_poly.h"

typedef struct
{
  fq_nmod_poly* entries;
  slong r;
  slong c;
  fq_nmod_poly** rows;
} fq_nmod_poly_mat_struct;

typedef fq_nmod_poly_mat_struct fq_nmod_poly_mat_t[1];

void fq_nmod_poly_mat_init(fq_nmod_poly_mat_t mat, slong rows, slong cols, const fq_nmod_ctx_t F);
void fq_nmod_poly_mat_clear(fq_nmod_poly_mat_t mat, const fq_nmod_ctx_t F);

#endif // __FQ_NMOD_POLY_MAT_H__
