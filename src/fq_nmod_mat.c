#include "flint/fq_nmod_mat.h"

#include "fq_nmod_mat.h"
#include "typedfes.h"

// This file implements additional utilities for performing linear algebra with fq_nmod_mat

bool fq_nmod_mat_is_zero_row(const fq_nmod_mat_t mat, slong row, const fq_nmod_ctx_t F)
{
  slong col;

  for (col = 0; col < fq_nmod_mat_ncols(mat); col++)
    if (!fq_nmod_is_zero(fq_nmod_mat_entry(mat, row, col),F) )
      return false;

  return true;
}

void fq_nmod_mat_swap_rows(fq_nmod_mat_t mat, slong row1, slong row2, const fq_nmod_ctx_t F)
{
  slong col;
  fq_nmod_t tmp;

  fq_nmod_init(tmp, F);

  for (col = 0; col < fq_nmod_mat_ncols(mat); col++) {
    fq_nmod_set(tmp, fq_nmod_mat_entry(mat, row1, col), F);
    fq_nmod_set(fq_nmod_mat_entry(mat, row1, col), fq_nmod_mat_entry(mat, row2, col), F);
    fq_nmod_set(fq_nmod_mat_entry(mat, row2, col), tmp, F);
  }
  
  fq_nmod_clear(tmp);

  return;
}

void fq_nmod_mat_swap_cols(fq_nmod_mat_t mat, slong col1, slong col2, const fq_nmod_ctx_t F)
{
  slong row;
  fq_nmod_t tmp;

  fq_nmod_init(tmp, F);

  for (row = 0; row < fq_nmod_mat_nrows(mat); row++) {
    fq_nmod_set(tmp, fq_nmod_mat_entry(mat, row, col1), F);
    fq_nmod_set(fq_nmod_mat_entry(mat, row, col1), fq_nmod_mat_entry(mat, row, col2), F);
    fq_nmod_set(fq_nmod_mat_entry(mat, row, col2), tmp, F);
  }
  
  fq_nmod_clear(tmp);

  return;
}

void fq_nmod_mat_add_row(fq_nmod_mat_t mat, slong dst_row, slong src_row, const fq_nmod_t scalar, const fq_nmod_ctx_t F)
{
  slong col;
  fq_nmod_t tmp;

  fq_nmod_init(tmp, F);

  for (col = 0; col < fq_nmod_mat_ncols(mat); col++) {
    fq_nmod_mul(tmp, fq_nmod_mat_entry(mat, src_row, col), scalar, F);
    fq_nmod_add(fq_nmod_mat_entry(mat, dst_row, col), fq_nmod_mat_entry(mat, dst_row, col), tmp, F);
  }

  fq_nmod_clear(tmp);
  
  return;
}

void fq_nmod_mat_add_col(fq_nmod_mat_t mat, slong dst_col, slong src_col, const fq_nmod_t scalar, const fq_nmod_ctx_t F)
{
  slong row;
  fq_nmod_t tmp;

  fq_nmod_init(tmp, F);

  for (row = 0; row < fq_nmod_mat_nrows(mat); row++) {
    fq_nmod_mul(tmp, fq_nmod_mat_entry(mat, row, src_col), scalar, F);
    fq_nmod_add(fq_nmod_mat_entry(mat, row, dst_col), fq_nmod_mat_entry(mat, row, dst_col), tmp, F);
  }

  fq_nmod_clear(tmp);
  
  return;
}

void fq_nmod_mat_mul_row(fq_nmod_mat_t mat, slong row, const fq_nmod_t scalar, const fq_nmod_ctx_t F)
{
  slong col;

  for (col = 0; col < fq_nmod_mat_ncols(mat); col++) {
    fq_nmod_mul(fq_nmod_mat_entry(mat,row,col),fq_nmod_mat_entry(mat,row,col),scalar,F);
  }

  return;
}

void fq_nmod_mat_mul_col(fq_nmod_mat_t mat, slong col, const fq_nmod_t scalar, const fq_nmod_ctx_t F)
{
  slong row;

  for (row = 0; row < fq_nmod_mat_nrows(mat); row++) {
    fq_nmod_mul(fq_nmod_mat_entry(mat,row,col),fq_nmod_mat_entry(mat,row,col),scalar,F);
  }

  return;
}
