#ifndef __FQ_NMOD_MPOLY_MAT_H__
#define __FQ_NMOD_MPOLY_MAT_H__

#ifdef FQ_NMOD_MPOLY_MAT_INLINES_C
#define FQ_NMOD_MPOLY_MAT_INLINE FLINT_DLL
#else
#define FQ_NMOD_MPOLY_MAT_INLINE static __inline__
#endif

#include "flint/flint.h"
#include "flint/fq_nmod_mpoly.h"

typedef struct
{
  fq_nmod_mpoly_struct* entries;
  slong r;
  slong c;
  fq_nmod_mpoly_struct** rows;
} fq_nmod_mpoly_mat_struct;

typedef fq_nmod_mpoly_mat_struct fq_nmod_mpoly_mat_t[1];

void fq_nmod_mpoly_mat_init(fq_nmod_mpoly_mat_t mat, slong rows, slong cols, const fq_nmod_mpoly_ctx_t R);
void fq_nmod_mpoly_mat_clear(fq_nmod_mpoly_mat_t mat, const fq_nmod_mpoly_ctx_t R);
void fq_nmod_mpoly_mat_print(const fq_nmod_mpoly_mat_t mat, const char** var_names, const fq_nmod_mpoly_ctx_t R);

FQ_NMOD_MPOLY_MAT_INLINE
fq_nmod_mpoly_struct * fq_nmod_mpoly_mat_entry(const fq_nmod_mpoly_mat_t mat, slong i, slong j)
{
   return mat->rows[i] + j;
}

FQ_NMOD_MPOLY_MAT_INLINE
slong fq_nmod_mpoly_mat_nrows(const fq_nmod_mpoly_mat_t mat)
{
   return mat->r;
}

FQ_NMOD_MPOLY_MAT_INLINE
slong fq_nmod_mpoly_mat_ncols(const fq_nmod_mpoly_mat_t mat)
{
   return mat->c;
}

#endif // __FQ_NMOD_MPOLY_MAT_H__
