#include <assert.h>

#include <carat/matrix.h>

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>

#include "fq_nmod_mat.h"
#include "matrix_tools.h"
#include "spinor.h"
#include "typedefs.h"

void spinor_init(spinor_t spinor, const fmpz_mat_t q)
{
  slong prime_idx, idx, n;
  fmpz_t tmp, det;
  fq_nmod_mat_t q_p;
  fmpz_factor_t bad_primes;

  n = fmpz_mat_nrows(q);
  assert(n == fmpz_mat_ncols(q));

  fmpz_init(det);
  fmpz_mat_det(det, q);
  fmpz_divexact_si(det, det, 2);
  fmpz_factor_init(bad_primes);
  fmpz_factor(bad_primes, det);
  fmpz_clear(det);
  
  spinor->num_primes = bad_primes->num;
  spinor->twist = (1LL << (bad_primes->num)) - 1;
  spinor->rads = (fq_nmod_mat_t*)malloc((bad_primes->num) * sizeof(fq_nmod_mat_t));
  spinor->fields = (fq_nmod_ctx_t*)malloc((bad_primes->num) * sizeof(fq_nmod_ctx_t));

  // !! TODO - when n is not squarefree take only the ones with odd exponents
  for (prime_idx = 0; prime_idx < bad_primes->num; prime_idx++) {
    fq_nmod_ctx_init(spinor->fields[prime_idx], &(bad_primes->p[prime_idx]), 1, "1");
    fq_nmod_mat_init_set_fmpz_mat(q_p, q, spinor->fields[prime_idx]);
    if (fmpz_equal_si(&(bad_primes->p[prime_idx]),2)) {
      fmpz_init(tmp);
      for (idx = 0; idx < n; idx++) {
	fmpz_divexact_si(tmp, fmpz_mat_entry(q, idx, idx), 2); 
	fq_nmod_set_fmpz(fq_nmod_mat_entry(q_p,idx,idx), tmp, spinor->fields[prime_idx]);
      }
      fmpz_clear(tmp);
    }
    fq_nmod_mat_kernel(spinor->rads[prime_idx], q_p, spinor->fields[prime_idx]);
    fq_nmod_mat_clear(q_p, spinor->fields[prime_idx]);
  }
  fmpz_factor_clear(bad_primes);
  return;
}

void spinor_clear(spinor_t spinor)
{
  slong prime_idx;
  
  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++) {
    fq_nmod_mat_clear(spinor->rads[prime_idx], spinor->fields[prime_idx]);
  }
  free(spinor->rads);
  free(spinor->fields);

  return;
}

W64 spinor_compute_vals(const spinor_t spinor, const fq_nmod_t* a)
{
  W64 val, mask;
  slong i;
#ifdef DEBUG
  fq_nmod_t neg;
#endif // DEBUG

  val = 0;
  mask = 1;
  for (i = 0; i < spinor->num_primes; i++) {
    if (!fq_nmod_is_one(a[i],spinor->fields[i])) {
#ifdef DEBUG
      fq_nmod_init(neg, spinor->fields[i]);
      fq_nmod_neg(neg, a[i], spinor->fields[i]);
      assert(fq_nmod_is_one(neg, spinor->fields[i]));
      fq_nmod_clear(neg, spinor->fields[i]);
#endif // DEBUG
      val ^= mask;
    }
    mask <<= 1;
  }
  return val;
}

W64 spinor_norm(const spinor_t spinor, matrix_TYP* mat, int denom)
{
  fmpz_mat_t mat_fmpz;
  fmpz_t denom_fmpz;
  W64 val;

  fmpz_init_set_si(denom_fmpz, denom);
  fmpz_mat_init_set_matrix_TYP(mat_fmpz, mat);
  
  val = spinor_norm_fmpz_mat(spinor, mat_fmpz, denom_fmpz);

  fmpz_mat_clear(mat_fmpz);
  fmpz_clear(denom_fmpz);
  return val;
}

// Here mat/denom is the isometry, and we're computing the spinor norm at all spinor primes
W64 spinor_norm_fmpz_mat(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom)
{
  fq_nmod_t denom_p;
  fq_nmod_mat_t mat_p, mat_p_t, rad_mat;
  slong idx, row, col, pivot, prime_idx, n;
  fmpz_t det; // for some reason det is not implemented for fq_nmod_mat. we reduce the det mod p instead
  fq_nmod_t det_p;
  fq_nmod_t* evs;
  W64 norm;

  n = fmpz_mat_nrows(mat);
  assert(n == fmpz_mat_ncols(mat));
  
  evs = (fq_nmod_t*)malloc((spinor->num_primes)*sizeof(fq_nmod_t));

  fmpz_init(det);
  fmpz_mat_det(det, mat);
  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++) {
    fq_nmod_mat_init_set_fmpz_mat(mat_p, mat, spinor->fields[prime_idx]);
    fq_nmod_mat_init(mat_p_t, n, n, spinor->fields[prime_idx]);
    fq_nmod_init(denom_p, spinor->fields[prime_idx]);
    fq_nmod_set_fmpz(denom_p, denom, spinor->fields[prime_idx]);
    // at the moment, when mat has scale divisible by p (at the bad primes)
    // we replace it by the trivial matrix. In general, should implement the spinor norm using reflection decomposition
    if (fq_nmod_is_zero(denom_p, spinor->fields[prime_idx])) {
      fq_nmod_mat_one(mat_p, spinor->fields[prime_idx]);
    }
    else {
      fq_nmod_inv(denom_p, denom_p, spinor->fields[prime_idx]);
      // scaling the reduced matrix - should be replaced bu a one-line function
      for (row = 0; row < fq_nmod_mat_nrows(mat_p, spinor->fields[prime_idx]); row++)
	for (col = 0; col < fq_nmod_mat_ncols(mat_p, spinor->fields[prime_idx]); col++)
	  fq_nmod_mul(fq_nmod_mat_entry(mat_p,row,col), fq_nmod_mat_entry(mat_p,row,col),
		      denom_p, spinor->fields[prime_idx]);
    }

    // !! TODO - we could instead multiply the end result by the determinant?
    fq_nmod_init(det_p, spinor->fields[prime_idx]);
    fq_nmod_set_fmpz(det_p, det, spinor->fields[prime_idx]);
    // scaling the determinant
    for (idx = 0; idx < n; idx++)
      fq_nmod_mul(det_p, det_p, denom_p, spinor->fields[prime_idx]);
    // If the matrix has determinant -1, we modify it to have a matrix in SO
    if (!fq_nmod_is_one(det_p, spinor->fields[prime_idx])) {
      fq_nmod_mat_neg(mat_p, mat_p, spinor->fields[prime_idx]);
    }
    fq_nmod_clear(det_p, spinor->fields[prime_idx]);

#ifdef DEBUG_LEVEL_FULL
    printf("rad = \n");
    fq_nmod_mat_print_pretty(spinor->rads[prime_idx], spinor->fields[prime_idx]);
    printf("\n");
    printf("s_p =  \n");
    fq_nmod_mat_print_pretty(mat_p, spinor->fields[prime_idx]);
    printf("\n");
    printf("scale (inverse) =  \n");
    fq_nmod_print_pretty(denom_p, spinor->fields[prime_idx]);
    printf("\n");
#endif // DEBUG_LEVEL_FULL
    
    fq_nmod_mat_init(rad_mat, fq_nmod_mat_nrows(spinor->rads[prime_idx], spinor->fields[prime_idx]),
		     fq_nmod_mat_ncols(mat_p, spinor->fields[prime_idx]), spinor->fields[prime_idx]);
    // we need to transpose, given the direction of our isometries
    fq_nmod_mat_transpose(mat_p_t, mat_p, spinor->fields[prime_idx]);
    fq_nmod_mat_mul(rad_mat, spinor->rads[prime_idx], mat_p_t, spinor->fields[prime_idx]);

    // for now assume rad is a single vector (we choose our lattices this way)
    // !! TODO !! -  modify to determinant so it will work in the general case
    assert(fq_nmod_mat_nrows(spinor->rads[prime_idx], spinor->fields[prime_idx]) == 1);
  
    for (pivot = 0; pivot < fq_nmod_mat_ncols(spinor->rads[prime_idx], spinor->fields[prime_idx]);) {
      if (!(fq_nmod_is_zero(fq_nmod_mat_entry(spinor->rads[prime_idx], 0, pivot), spinor->fields[prime_idx]))) {
	break;
      }
      pivot++;
    }
    
    fq_nmod_init(evs[prime_idx], spinor->fields[prime_idx]);
    // !! TODO !! - we should always set the pivot in the kernel to be 1
    // so that no arithmetic will be needed here.
    // we can also save the kernels upon initialization
    fq_nmod_inv(evs[prime_idx], fq_nmod_mat_entry(spinor->rads[prime_idx], 0, pivot), spinor->fields[prime_idx]);
    fq_nmod_mul(evs[prime_idx], evs[prime_idx], fq_nmod_mat_entry(rad_mat, 0, pivot), spinor->fields[prime_idx]);
    
    //    norm = (fq_nmod_is_one(ev, spinor->fields[prime_idx]) ? 1 : -1);

#ifdef DEBUG
    if (!fq_nmod_is_one(evs[prime_idx], spinor->fields[prime_idx]) ) {
      fq_nmod_init(det_p, spinor->fields[prime_idx]);
      fq_nmod_neg(det_p, evs[prime_idx], spinor->fields[prime_idx]);
      assert(fq_nmod_is_one(det_p, spinor->fields[prime_idx]));
      fq_nmod_clear(det_p, spinor->fields[prime_idx]);
    }
#endif // DEBUG

    fq_nmod_mat_clear(rad_mat, spinor->fields[prime_idx]);
    fq_nmod_clear(denom_p, spinor->fields[prime_idx]);
    fq_nmod_mat_clear(mat_p_t, spinor->fields[prime_idx]);
    fq_nmod_mat_clear(mat_p, spinor->fields[prime_idx]);
  }
  fmpz_clear(det);

  norm = spinor_compute_vals(spinor, evs);

  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++)
    fq_nmod_clear(evs[prime_idx], spinor->fields[prime_idx]);
  
  free(evs);
  
  return norm;
}
