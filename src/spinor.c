#include <assert.h>

#include <carat/matrix.h>

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>

#include "arith.h"
#include "fq_nmod_mat.h"
#include "matrix_tools.h"
#include "spinor.h"
#include "typedefs.h"

void spinor_init(spinor_t spinor, const fmpz_mat_t q)
{
  slong prime_idx /*, idx */;
  fmpz_t /* tmp, */ det;
  fq_nmod_mat_t q_p;
  fmpz_factor_t bad_primes;

  fmpz_init(det);
  fmpz_mat_det(det, q);
  fmpz_divexact_si(det, det, 2);
  fmpz_factor_init(bad_primes);
  fmpz_factor(bad_primes, det);
  fmpz_clear(det);

  fmpz_mat_init_set(spinor->Q,q);
  spinor->num_primes = bad_primes->num;
  spinor->twist = (1LL << (bad_primes->num)) - 1;
  spinor->rads = (fq_nmod_mat_t*)malloc((bad_primes->num) * sizeof(fq_nmod_mat_t));
  spinor->fields = (fq_nmod_ctx_t*)malloc((bad_primes->num) * sizeof(fq_nmod_ctx_t));

  // !! TODO - when n is not squarefree take only the ones with odd exponents
  for (prime_idx = 0; prime_idx < bad_primes->num; prime_idx++) {
    fq_nmod_ctx_init(spinor->fields[prime_idx], &(bad_primes->p[prime_idx]), 1, "1");
    fq_nmod_mat_init_set_fmpz_mat(q_p, q, spinor->fields[prime_idx]);
    // Is this necessary? Check!!
    /*
    if (fmpz_equal_si(&(bad_primes->p[prime_idx]),2)) {
      fmpz_init(tmp);
      for (idx = 0; idx < n; idx++) {
	fmpz_divexact_si(tmp, fmpz_mat_entry(q, idx, idx), 2); 
	fq_nmod_set_fmpz(fq_nmod_mat_entry(q_p,idx,idx), tmp, spinor->fields[prime_idx]);
      }
      fmpz_clear(tmp);
    }
    */
    fq_nmod_mat_kernel(spinor->rads[prime_idx], q_p, spinor->fields[prime_idx]);
    fq_nmod_mat_clear(q_p, spinor->fields[prime_idx]);
  }
  fmpz_factor_clear(bad_primes);
  return;
}

void spinor_clear(spinor_t spinor)
{
  slong prime_idx;

  fmpz_mat_clear(spinor->Q);
  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++) {
    fq_nmod_mat_clear(spinor->rads[prime_idx], spinor->fields[prime_idx]);
    fq_nmod_ctx_clear(spinor->fields[prime_idx]);
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

  if (fmpz_equal_si(fq_nmod_ctx_prime(spinor->fields[0]),2))
    val = spinor_norm_cd_fmpz_mat(spinor, mat_fmpz, denom_fmpz);
  else {
    val = spinor_norm_fmpz_mat(spinor, mat_fmpz, denom_fmpz);
    assert(val==spinor_norm_cd_fmpz_mat(spinor, mat_fmpz, denom_fmpz));
  }

  fmpz_mat_clear(mat_fmpz);
  fmpz_clear(denom_fmpz);
  return val;
}

W64 spinor_norm_isom(const spinor_t spinor, const isometry_t isom)
{
  fmpz_mat_t mat_fmpz;
  fmpz_t denom_fmpz;
  W64 val;

  fmpz_init_set_si(denom_fmpz, isom->denom);
  fmpz_mat_init_set_square_matrix(mat_fmpz, isom->s);

  if (fmpz_equal_si(fq_nmod_ctx_prime(spinor->fields[0]),2))
    val = spinor_norm_cd_fmpz_mat(spinor, mat_fmpz, denom_fmpz);
  else {
    val = spinor_norm_fmpz_mat(spinor, mat_fmpz, denom_fmpz);
    assert(val==spinor_norm_cd_fmpz_mat(spinor, mat_fmpz, denom_fmpz));
  }

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
    // we compute the Cartan-Dieudonne instead
    if (fq_nmod_is_zero(denom_p, spinor->fields[prime_idx])) {
      // fq_nmod_mat_one(mat_p, spinor->fields[prime_idx]);
      fq_nmod_clear(denom_p, spinor->fields[prime_idx]);
      fq_nmod_mat_clear(mat_p, spinor->fields[prime_idx]);
      fq_nmod_mat_clear(mat_p_t, spinor->fields[prime_idx]);
      return spinor_norm_cd_fmpz_mat(spinor, mat, denom);
    }
    else {
      fq_nmod_inv(denom_p, denom_p, spinor->fields[prime_idx]);
      // scaling the reduced matrix - should be replaced bu a one-line function
      for (row = 0; row < fq_nmod_mat_nrows(mat_p, spinor->fields[prime_idx]); row++)
	for (col = 0; col < fq_nmod_mat_ncols(mat_p, spinor->fields[prime_idx]); col++)
	  fq_nmod_mul(fq_nmod_mat_entry(mat_p,row,col), fq_nmod_mat_entry(mat_p,row,col),
		      denom_p, spinor->fields[prime_idx]);
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
    }

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

// When p = 2 the above doesn't work (1 == -1). We write here the slow way of computing (the global) spinor norm until we figure out a better way

// builds a reflection through y with respect to A
// we assume y is given as a row vector, and A is symmetric
void build_reflection(fmpq_mat_t ref, const fmpq_mat_t y, const fmpq_mat_t A)
{
  slong dim;
  fmpq_mat_t yA, A_yt, A_yt_y, yA_yt;
  fmpq_t scalar;

  dim = fmpq_mat_nrows(A);
  assert(dim == fmpq_mat_ncols(A));
  assert(fmpq_mat_nrows(y) == 1);
  assert(fmpq_mat_ncols(y) == dim);
  
  fmpq_mat_init(yA, 1, dim);
  fmpq_mat_init(A_yt, dim, 1);
  fmpq_mat_init(A_yt_y, dim, dim);
  fmpq_mat_init(yA_yt, 1, 1);

  fmpq_mat_one(ref);
  
  fmpq_mat_mul(yA,y,A);
  fmpq_mat_transpose(A_yt, yA);
  fmpq_mat_mul(A_yt_y, A_yt, y);
  fmpq_mat_mul(yA_yt, y, A_yt);

  fmpq_init(scalar);

  fmpq_set_si(scalar, 2, 1);
  fmpq_div(scalar, scalar, fmpq_mat_entry(yA_yt, 0, 0));

  fmpq_mat_scalar_mul_fmpq(A_yt_y, A_yt_y, scalar);
  fmpq_mat_sub(ref, ref, A_yt_y);
  
  fmpq_clear(scalar);

  fmpq_mat_clear(yA_yt);
  fmpq_mat_clear(A_yt_y);
  fmpq_mat_clear(A_yt);
  fmpq_mat_clear(yA);
  return;
}

// here BM is the basis matrix of a subspace,
// and we are looking to return an isotropic vector of A in that subspace

bool get_aniso_vec(fmpq_mat_t x, const fmpq_mat_t BM, const fmpq_mat_t A)
{
  slong dim, dim_B, i, j, k;
  fmpq_mat_t A_res, BMA, BM_t;

  dim = fmpq_mat_nrows(A);
  assert(dim == fmpq_mat_ncols(A));
  assert(dim == fmpq_mat_ncols(BM));
  dim_B = fmpq_mat_nrows(BM);

  fmpq_mat_init(A_res, dim_B, dim_B);
  fmpq_mat_init(BMA, dim_B, dim);
  fmpq_mat_init(BM_t, dim, dim_B);

  fmpq_mat_transpose(BM_t, BM);
  fmpq_mat_mul(BMA, BM, A);
  fmpq_mat_mul(A_res, BMA, BM_t);

  if (fmpq_mat_is_zero(A_res)) {
    fmpq_mat_clear(BM_t);
    fmpq_mat_clear(BMA);
    fmpq_mat_clear(A_res);
    return false;
  }
  
  for (i = 0; i < dim_B; i++) {
    if (!fmpq_is_zero(fmpq_mat_entry(A_res, i, i))) {
      for (j = 0; j < dim; j++)
	fmpq_set(fmpq_mat_entry(x,0,j),fmpq_mat_entry(BM,i,j));
      fmpq_mat_clear(BM_t);
      fmpq_mat_clear(BMA);
      fmpq_mat_clear(A_res);
      return true;
    }
  }

  for (i = 0; i < dim_B; i++)
    for (j = i+1; j < dim_B; j++) {
      if (!fmpq_is_zero(fmpq_mat_entry(A_res, i, j))) {
	for (k = 0; k < dim; k++) {
	  fmpq_set(fmpq_mat_entry(x,0,k),fmpq_mat_entry(BM,i,k));
	  fmpq_add(fmpq_mat_entry(x,0,k), fmpq_mat_entry(x,0,k), fmpq_mat_entry(BM,j,k)); 
	}
	fmpq_mat_clear(BM_t);
	fmpq_mat_clear(BMA);
	fmpq_mat_clear(A_res);
	return true;
      }
    }
  
  fmpq_mat_clear(BM_t);
  fmpq_mat_clear(BMA);
  fmpq_mat_clear(A_res);
  
  return false;
}


// reflections with respect to A - Cartan-Dieudonne
// every row of the returned matrix pts is the coordinates of a reflection point
void find_reflection_pts(fmpq_mat_t pts, const fmpq_mat_t s, const fmpq_mat_t A)
{
  bool aniso;
  slong dim, i, j;
  fmpq_mat_t sA, s_t, sAs_t, s_m_1, fixed, x, xA, At_xt, H, s_H, H_t, HA, A_H, pts_H, t, ts;

  dim = fmpq_mat_nrows(A);
  assert(dim == fmpq_mat_ncols(A));
  assert(dim == fmpq_mat_nrows(s));
  assert(dim == fmpq_mat_ncols(s));
  
  fmpq_mat_init(s_t, dim, dim);
  fmpq_mat_init(sA, dim, dim);
  fmpq_mat_init(sAs_t, dim, dim);

  fmpq_mat_mul(sA, s, A);
  fmpq_mat_transpose(s_t, s);
  fmpq_mat_mul(sAs_t, sA, s_t);
  
  assert(fmpq_mat_equal(A, sAs_t));

  fmpq_mat_clear(sAs_t);
  fmpq_mat_clear(s_t);
  fmpq_mat_clear(sA);
  
  if (fmpq_mat_is_one(s)) {
    // is this a thing?
    fmpq_mat_init(pts, 0, 0); 
    return;
  }

  if ((dim == 1) && fmpq_equal_si(fmpq_mat_entry(s, 0, 0), -1)) {
    fmpq_mat_init(pts, 1, 1);
    fmpq_one(fmpq_mat_entry(pts, 0, 0));
    return;
  }

  fmpq_mat_init(s_m_1, dim, dim);
  fmpq_mat_one(s_m_1);
  fmpq_mat_sub(s_m_1, s, s_m_1);
  
  fmpq_mat_left_kernel(fixed, s_m_1);

  fmpq_mat_init(x, 1, dim);
  aniso = get_aniso_vec(x, fixed, A);

  if (aniso) {
    // case 1
    fmpq_mat_init(xA, 1, dim);
    fmpq_mat_init(At_xt, dim, 1);
    fmpq_mat_mul(xA, x, A);
    fmpq_mat_transpose(At_xt, xA);
    fmpq_mat_left_kernel(H, At_xt);

    fmpq_mat_init(s_H, fmpq_mat_nrows(H), fmpq_mat_nrows(H));
    restrict_mat(s_H, s, H);

    fmpq_mat_init(HA, fmpq_mat_nrows(H), dim);
    fmpq_mat_init(H_t, fmpq_mat_ncols(H), fmpq_mat_nrows(H));
    fmpq_mat_transpose(H_t, H);
    fmpq_mat_mul(HA, H, A);
    fmpq_mat_init(A_H, fmpq_mat_nrows(H), fmpq_mat_nrows(H));
    fmpq_mat_mul(A_H, HA, H_t);

    find_reflection_pts(pts_H, s_H, A_H);

    fmpq_mat_init(pts, fmpq_mat_nrows(pts_H), dim);

    fmpq_mat_mul(pts, pts_H, H);

    fmpq_mat_clear(pts_H);
  
    fmpq_mat_clear(A_H);
    fmpq_mat_clear(HA);
    fmpq_mat_clear(H_t);
  
    fmpq_mat_clear(s_H);
    fmpq_mat_clear(H);
    fmpq_mat_clear(At_xt);
    fmpq_mat_clear(xA);
  }
  else {
    aniso = get_aniso_vec(x, s_m_1, A);
    if (!aniso) {
      // case 3
      assert(dim == 4);
      aniso = get_aniso_vec(x, A, A);
      assert(aniso);
    }
   
    // case 2
    fmpq_mat_init(t, dim, dim);
    fmpq_mat_init(ts, dim, dim);
    build_reflection(t, x, A);
    fmpq_mat_mul(ts, t, s);
    find_reflection_pts(pts_H, ts, A);
    fmpq_mat_init(pts, fmpq_mat_nrows(pts_H)+1, dim);
    for (j = 0; j < dim; j++)
      fmpq_set(fmpq_mat_entry(pts,0,j), fmpq_mat_entry(x,0,j));
    for (i = 0; i < fmpq_mat_nrows(pts_H); i++)
      for (j = 0; j < dim; j++)
	fmpq_set(fmpq_mat_entry(pts,i+1,j), fmpq_mat_entry(pts_H,i,j));
    fmpq_mat_clear(pts_H);
    fmpq_mat_clear(ts);
    fmpq_mat_clear(t);
   
  }

#ifdef DEBUG
  fmpq_mat_init(H, dim, dim);
  fmpq_mat_one(H);
  for (i = 0; i < fmpq_mat_nrows(pts); i++) {
    fmpq_mat_init(t, dim, dim);
    fmpq_mat_init(ts, 1, dim);
    for (j = 0; j < dim; j++)
      fmpq_set(fmpq_mat_entry(ts,0,j), fmpq_mat_entry(pts,i,j));
    build_reflection(t, ts, A);
    fmpq_mat_mul(H, H, t);
    fmpq_mat_clear(ts);
    fmpq_mat_clear(t);
  }
  assert(fmpq_mat_equal(H, s));
  fmpq_mat_clear(H);
#endif // DEBUG
  
  fmpq_mat_clear(x);
  fmpq_mat_clear(fixed);
  fmpq_mat_clear(s_m_1);
  
  return;
}

void spinor_norm_cd(fmpq_t norm, const fmpq_mat_t s, const fmpq_mat_t A)
{
  fmpq_mat_t pts, pA, pAp_t, pts_t, s_t;
  slong num_pts, dim, i;

  dim = fmpq_mat_nrows(A);
  fmpq_mat_init(s_t, dim, dim);
  fmpq_mat_transpose(s_t, s);
  
  fmpq_one(norm);
  
  find_reflection_pts(pts,s_t,A);
  num_pts = fmpq_mat_nrows(pts);  

  // we are overworking here, but it is slightly more convenient at the moment
  
  fmpq_mat_init(pA, num_pts, dim);
  fmpq_mat_init(pts_t, dim, num_pts);
  fmpq_mat_init(pAp_t, num_pts, num_pts);

  fmpq_mat_transpose(pts_t, pts);
  fmpq_mat_mul(pA, pts, A);
  fmpq_mat_mul(pAp_t, pA, pts_t);
  
  for (i = 0; i < num_pts; i++) {
    fmpq_mul(norm, norm, fmpq_mat_entry(pAp_t, i, i));
  }

  fmpq_mat_clear(pts);
  fmpq_mat_clear(pA);
  fmpq_mat_clear(s_t);
  fmpq_mat_clear(pts_t);
  fmpq_mat_clear(pAp_t);
  
  return;
}

// this is not really useful
int rho(const fmpz_t p, const fmpq_mat_t s, const fmpq_mat_t A)
{
  fmpq_t norm;
  int ret;

  fmpq_init(norm);
  
  spinor_norm_cd(norm, s, A);

  ret = fmpq_valuation(norm, p) & 1; // if odd return 1, else 0;

  fmpq_clear(norm);

  return ret;
}

W64 spinor_norm_cd_fmpz_mat(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom)
{
  fmpq_mat_t mat_q, Q_q;
  fmpq_t det;
  fmpq_t spinorm;
  slong prime_idx, n;
  
  W64 val;

  n = fmpz_mat_nrows(mat);
  assert(n == fmpz_mat_ncols(mat));
  
  fmpq_mat_init(mat_q, n, n);
  fmpq_mat_init(Q_q, n, n);
  fmpq_mat_set_fmpz_mat_div_fmpz(mat_q, mat, denom);
  fmpq_mat_set_fmpz_mat(Q_q, spinor->Q);

  fmpq_init(det);
  fmpq_mat_det(det, mat_q);
  if (fmpq_equal_si(det, -1))
    fmpq_mat_neg(mat_q, mat_q);

  fmpq_init(spinorm);
  spinor_norm_cd(spinorm, mat_q, Q_q);
  val = 0;
  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++) {
    val ^= (fmpq_valuation(spinorm, fq_nmod_ctx_prime(spinor->fields[prime_idx])) & 1) << prime_idx;
  }
  fmpq_clear(spinorm);
  fmpq_clear(det);
  fmpq_mat_clear(mat_q);
  fmpq_mat_clear(Q_q);

  return val;
}
