#include <assert.h>

#include "flint/fq_nmod_mat.h"

#include "fq_nmod_mat.h"
#include "typedefs.h"

// This file implements additional utilities for performing linear algebra with fq_nmod_mat

bool fq_nmod_mat_is_zero_row(const fq_nmod_mat_t mat, slong row, const fq_nmod_ctx_t F)
{
  slong col;

  for (col = 0; col < fq_nmod_mat_ncols(mat,F); col++)
    if (!fq_nmod_is_zero(fq_nmod_mat_entry(mat, row, col),F) )
      return false;

  return true;
}

// these seem to already exist

/* void fq_nmod_mat_swap_rows(fq_nmod_mat_t mat, slong row1, slong row2, const fq_nmod_ctx_t F) */
/* { */
/*   slong col; */
/*   fq_nmod_t tmp; */

/*   fq_nmod_init(tmp, F); */

/*   for (col = 0; col < fq_nmod_mat_ncols(mat); col++) { */
/*     fq_nmod_set(tmp, fq_nmod_mat_entry(mat, row1, col), F); */
/*     fq_nmod_set(fq_nmod_mat_entry(mat, row1, col), fq_nmod_mat_entry(mat, row2, col), F); */
/*     fq_nmod_set(fq_nmod_mat_entry(mat, row2, col), tmp, F); */
/*   } */
  
/*   fq_nmod_clear(tmp); */

/*   return; */
/* } */

/* void fq_nmod_mat_swap_cols(fq_nmod_mat_t mat, slong col1, slong col2, const fq_nmod_ctx_t F) */
/* { */
/*   slong row; */
/*   fq_nmod_t tmp; */

/*   fq_nmod_init(tmp, F); */

/*   for (row = 0; row < fq_nmod_mat_nrows(mat); row++) { */
/*     fq_nmod_set(tmp, fq_nmod_mat_entry(mat, row, col1), F); */
/*     fq_nmod_set(fq_nmod_mat_entry(mat, row, col1), fq_nmod_mat_entry(mat, row, col2), F); */
/*     fq_nmod_set(fq_nmod_mat_entry(mat, row, col2), tmp, F); */
/*   } */
  
/*   fq_nmod_clear(tmp); */

/*   return; */
/* } */

void fq_nmod_mat_add_row(fq_nmod_mat_t mat, slong dst_row, slong src_row, const fq_nmod_t scalar, const fq_nmod_ctx_t F)
{
  slong col;
  fq_nmod_t tmp;

  fq_nmod_init(tmp, F);

  for (col = 0; col < fq_nmod_mat_ncols(mat,F); col++) {
    fq_nmod_mul(tmp, fq_nmod_mat_entry(mat, src_row, col), scalar, F);
    fq_nmod_add(fq_nmod_mat_entry(mat, dst_row, col), fq_nmod_mat_entry(mat, dst_row, col), tmp, F);
  }

  fq_nmod_clear(tmp,F);
  
  return;
}

void fq_nmod_mat_add_col(fq_nmod_mat_t mat, slong dst_col, slong src_col, const fq_nmod_t scalar, const fq_nmod_ctx_t F)
{
  slong row;
  fq_nmod_t tmp;

  fq_nmod_init(tmp, F);

  for (row = 0; row < fq_nmod_mat_nrows(mat,F); row++) {
    fq_nmod_mul(tmp, fq_nmod_mat_entry(mat, row, src_col), scalar, F);
    fq_nmod_add(fq_nmod_mat_entry(mat, row, dst_col), fq_nmod_mat_entry(mat, row, dst_col), tmp, F);
  }

  fq_nmod_clear(tmp,F);
  
  return;
}

void fq_nmod_mat_mul_row(fq_nmod_mat_t mat, slong row, const fq_nmod_t scalar, const fq_nmod_ctx_t F)
{
  slong col;

  for (col = 0; col < fq_nmod_mat_ncols(mat,F); col++) {
    fq_nmod_mul(fq_nmod_mat_entry(mat,row,col),fq_nmod_mat_entry(mat,row,col),scalar,F);
  }

  return;
}

void fq_nmod_mat_mul_col(fq_nmod_mat_t mat, slong col, const fq_nmod_t scalar, const fq_nmod_ctx_t F)
{
  slong row;

  for (row = 0; row < fq_nmod_mat_nrows(mat,F); row++) {
    fq_nmod_mul(fq_nmod_mat_entry(mat,row,col),fq_nmod_mat_entry(mat,row,col),scalar,F);
  }

  return;
}

void fq_nmod_mat_init_set_fmpz_mat(fq_nmod_mat_t dest, const fmpz_mat_t mat, const fq_nmod_ctx_t F)
{
  slong row, col;

  fq_nmod_mat_init(dest, fmpz_mat_nrows(mat), fmpz_mat_ncols(mat), F);
  
  for (row = 0; row < fmpz_mat_nrows(mat); row++)
    for (col = 0; col < fmpz_mat_ncols(mat); col++)
      fq_nmod_set_fmpz(fq_nmod_mat_entry(dest, row, col), fmpz_mat_entry(mat, row, col), F);

  return;
}


void fq_nmod_mat_transpose(fq_nmod_mat_t mat_t, const fq_nmod_mat_t mat, const fq_nmod_ctx_t F)
{
  slong row,col;

  assert(mat_t != mat);
  assert(fq_nmod_mat_nrows(mat_t,F) == fq_nmod_mat_ncols(mat,F));
  assert(fq_nmod_mat_nrows(mat,F) == fq_nmod_mat_ncols(mat_t,F));

  for (row = 0; row < fq_nmod_mat_nrows(mat,F); row++)
    for (col = 0; col < fq_nmod_mat_ncols(mat,F); col++)
      fq_nmod_set(fq_nmod_mat_entry(mat_t,col,row),fq_nmod_mat_entry(mat,row,col),F);

  return;
}

void fq_nmod_mat_rref_trans(fq_nmod_mat_t mat, fq_nmod_mat_t trans, const fq_nmod_ctx_t F)
{
  slong* P;
  slong row, col, nrows, ncols, pivot, prev_row;
  fq_nmod_mat_t LU, L_inv;
  fq_nmod_t scalar;

#ifdef DEBUG
  fq_nmod_mat_t test1, test2;
#endif // DEBUG

  nrows = fq_nmod_mat_nrows(mat, F);
  ncols = fq_nmod_mat_ncols(mat, F);
  P = (slong*)malloc(nrows*sizeof(slong));

#ifdef DEBUG
  fq_nmod_mat_init_set(test1, mat, F);
#endif // DEBUG
  
  fq_nmod_mat_init_set(LU, mat, F);
  
  fq_nmod_mat_lu(P, LU, false, F); // returns L\U in LU such that PA = LU, P permutation matrix

  // initializing L_inv to be 1
  fq_nmod_mat_init(L_inv, nrows, nrows, F);
  fq_nmod_mat_zero(L_inv,F);
  for (row = 0; row < nrows; row++)
    fq_nmod_one(fq_nmod_mat_entry(L_inv, row, row),F); 

  // sets L_inv to L^(-1)
  fq_nmod_mat_solve_tril(L_inv, LU, L_inv, 1, F);

  // setting the transformation to be L^(-1)*P
  assert((fq_nmod_mat_nrows(trans,F) == nrows) && (fq_nmod_mat_ncols(trans,F) == nrows));
  
  for (row = 0; row < nrows; row++) {
    for (col = 0; col < nrows; col++) {
      fq_nmod_set(fq_nmod_mat_entry(trans,row,P[col]), fq_nmod_mat_entry(L_inv,row,col),F);
    }
  }

  // setting mat to be U
  for (row = 0; row < nrows; row++) {
    for (col = 0; col < row; col++) {
      fq_nmod_zero(fq_nmod_mat_entry(mat,row,col), F);
    }
    for (col = row; col < ncols; col++) {
      fq_nmod_set(fq_nmod_mat_entry(mat,row,col), fq_nmod_mat_entry(LU,row,col),F);
    }
  }

#ifdef DEBUG
  fq_nmod_mat_init(test2, nrows, nrows, F);
  fq_nmod_mat_mul(test2, trans, test1, F);

  assert(fq_nmod_mat_equal(test2, mat, F));
  
  fq_nmod_mat_clear(test1, F);
  fq_nmod_mat_clear(test2, F);
#endif // DEBUG

#ifdef DEBUG_LEVEL_FULL
  printf("Before zeroing above the pivots, get: \n");
  fq_nmod_mat_print_pretty(mat,F);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  // zeroing above the pivots
  fq_nmod_init(scalar, F);
  for (row = 0; row < nrows; row++) {
    for (pivot = 0; (pivot < ncols) && (fq_nmod_is_zero(fq_nmod_mat_entry(mat,row,pivot),F)); pivot++);
    assert((pivot == ncols) || fq_nmod_is_one(fq_nmod_mat_entry(mat,row,pivot),F));
    if (pivot < ncols) {
      for (prev_row = 0; prev_row < row; prev_row++) {
	fq_nmod_neg(scalar, fq_nmod_mat_entry(mat,prev_row,pivot), F);
	fq_nmod_mat_add_row(mat,prev_row,row,scalar,F);
	fq_nmod_mat_add_row(trans,prev_row,row,scalar,F);
      }
    }
  }
  fq_nmod_clear(scalar,F);

  fq_nmod_mat_clear(LU, F);
  fq_nmod_mat_clear(L_inv, F);
  free(P);	      
  return;
}

// !! TODO - reuse the code from above
void fq_nmod_mat_kernel(fq_nmod_mat_t ker, const fq_nmod_mat_t mat, const fq_nmod_ctx_t F)
{
  slong* P;
  slong rank, row, col, n;
  fq_nmod_mat_t LU, L_inv;
#ifdef DEBUG_LEVEL_FULL
  fq_nmod_mat_t zero;
#endif // DEBUG_LEVEL_FULL

  n = fq_nmod_mat_nrows(mat, F);
  P = (slong*)malloc(n*sizeof(slong));
  
  fq_nmod_mat_init_set(LU, mat, F);
  
  rank = fq_nmod_mat_lu(P, LU, false, F); // returns L\U in LU such that PA = LU, P permutation matrix

  // initializing L_inv to be 1
  fq_nmod_mat_init(L_inv, n, n, F);
  fq_nmod_mat_zero(L_inv,F);
  for (row = 0; row < n; row++)
    fq_nmod_one(fq_nmod_mat_entry(L_inv, row, row),F); 

  // sets L_inv to L^(-1)
  fq_nmod_mat_solve_tril(L_inv, LU, L_inv, 1, F);

  // setting the kernel to the last n-r rows of L^(-1)*P
  fq_nmod_mat_init(ker, n - rank, n, F);
  fq_nmod_mat_zero(ker,F);
  for (row = rank; row < n; row++) {
    for (col = 0; col < n; col++) {
      fq_nmod_set(fq_nmod_mat_entry(ker,row-rank,P[col]), fq_nmod_mat_entry(L_inv,row,col),F);
    }
  }

#ifdef DEBUG_LEVEL_FULL
  fq_nmod_mat_init(zero, fq_nmod_mat_nrows(ker,F), fq_nmod_mat_ncols(mat,F),F);
  fq_nmod_mat_mul(zero, ker, mat, F);
  printf("ker = \n");
  fq_nmod_mat_print_pretty(ker, F);
  printf("\n");
  printf("mat = \n");
  fq_nmod_mat_print_pretty(mat, F);
  printf("ker * mat = \n");
  fq_nmod_mat_print_pretty(zero, F);
  printf("\n");
  printf("LU = \n");
  fq_nmod_mat_print_pretty(LU, F);
  printf("\n");
  printf("L_inv = \n");
  fq_nmod_mat_print_pretty(L_inv, F);
  printf("\n");
  assert(fq_nmod_mat_is_zero(zero, F));
  fq_nmod_mat_clear(zero,F);
#endif // DEBUG_LEVEL_FULL

  fq_nmod_mat_clear(LU, F);
  fq_nmod_mat_clear(L_inv, F);
  free(P);	      
  return;
}
