#include <assert.h>

#include <carat/matrix.h>

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>

#include "fq_nmod_mat.h"
#include "typedefs.h"

// Here mat/denom is the isometry, and we're computing the spinor norm at p
slong fmpz_mat_spinor_norm(const fmpz_mat_t q, const fmpz_mat_t mat, const fmpz_t denom, slong p_int)
{
  fq_nmod_ctx_t F;
  fq_nmod_t denom_p, ev;
  fq_nmod_mat_t q_p, mat_p, rad, rad_mat;
  fmpz_t tmp, p;
  slong idx, row, col, pivot, norm;
  
#ifdef DEBUG
  fmpz_t disc;
#endif // DEBUG

  fmpz_init_set_si(p, p_int);
  
#ifdef DEBUG
  // verify that p divides the discirminant
  fmpz_init(disc);
  fmpz_mat_det(disc, q);
  // (half) discriminant = det/2 when N = 5 
  fmpz_divexact_si(disc, disc, 2);
  assert(fmpz_divisible(disc, p));
  fmpz_clear(disc);
#endif // DEBUG
  
  fq_nmod_ctx_init(F, p, 1, "1");
  fq_nmod_mat_init_set_fmpz_mat(q_p, q, F);
  if (p_int == 2) {
    fmpz_init(tmp);
    for (idx = 0; idx < N; idx++) {
      fmpz_divexact_si(tmp, fmpz_mat_entry(q, idx, idx), 2); 
      fq_nmod_set_fmpz(fq_nmod_mat_entry(q_p,idx,idx), tmp, F);
    }
    fmpz_clear(tmp);
  }

  fq_nmod_mat_init_set_fmpz_mat(mat_p, mat, F);
  fq_nmod_init(denom_p, F);
  fq_nmod_set_fmpz(denom_p, denom, F);
  fq_nmod_inv(denom_p, denom_p, F);
  // scaling the reduced matrix
  for (row = 0; row < fq_nmod_mat_nrows(mat_p, F); row++)
    for (col = 0; col < fq_nmod_mat_ncols(mat_p, F); col++)
      fq_nmod_mul(fq_nmod_mat_entry(mat_p,row,col), fq_nmod_mat_entry(mat_p,row,col), denom_p, F);
  
  fq_nmod_mat_kernel(rad, q_p, F);

  fq_nmod_mat_init(rad_mat, fq_nmod_mat_nrows(rad, F), fq_nmod_mat_ncols(mat_p, F), F);
  fq_nmod_mat_mul(rad_mat, rad, mat_p, F);

  // for now assume rad is a single vector (we choose our lattices this way)
  // !! TODO !! -  modify to determinant so it will work in the general case
  assert(fq_nmod_mat_nrows(rad, F) == 1);
  
  for (pivot = 0; pivot < fq_nmod_mat_ncols(rad, F);) {
    if (!(fq_nmod_is_zero(fq_nmod_mat_entry(rad, 0, pivot), F))) {
      break;
    }
    pivot++;
  }

  fq_nmod_init(ev, F);
  // !! TODO !! - we should always set the pivot in the kernel to be 1
  // so that no arithmetic will be needed here.
  fq_nmod_inv(ev, fq_nmod_mat_entry(rad, 0, pivot), F);
  fq_nmod_mul(ev, ev, fq_nmod_mat_entry(rad_mat, 0, pivot), F);

  norm = (fq_nmod_is_one(ev, F) ? 1 : -1);

#ifdef DEBUG
  if (norm == -1) {
    fq_nmod_neg(ev, ev, F);
    assert(fq_nmod_is_one(ev, F));
  }
#endif // DEBUG

  fq_nmod_clear(ev, F);
  fq_nmod_mat_clear(rad_mat, F);
  fq_nmod_clear(denom_p, F);
  fq_nmod_mat_clear(rad, F);
  fq_nmod_mat_clear(mat_p,F);
  fq_nmod_mat_clear(q_p,F);
  fq_nmod_ctx_clear(F);
  fmpz_clear(p);

  return norm;
}
