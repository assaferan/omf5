#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <carat/matrix.h>

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>

#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>
#include <flint/nmod_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_poly.h>

#include "arith.h"
#include "fq_nmod_mat.h"
#include "matrix_tools.h"
#include "spinor.h"
#include "typedefs.h"

slong nmod_mat_nullspace_mod4(nmod_mat_t ker, const nmod_mat_t q)
{
  slong nullity, i, j, pivot;
  fq_nmod_mat_t q_2, ker_2, lift_2, trans, trans_t, sol;
  nmod_mat_t lift, lift_q;
  fq_nmod_ctx_t F2;
  fq_nmod_t elt;
  fmpz_t two;

#ifdef DEBUG
  fq_nmod_mat_t test1, test2;
  nmod_mat_t zero, ker_t;
#endif // DEBUG

  fmpz_init(two);
  fmpz_set_si(two, 2);

  fq_nmod_ctx_init(F2, two, 1, "1");
  fq_nmod_mat_init(q_2, QF_RANK, QF_RANK, F2);
  fq_nmod_init(elt, F2);
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++) {
      fq_nmod_set_si(elt, nmod_mat_get_entry(q, i, j), F2);
      fq_nmod_set(fq_nmod_mat_entry(q_2, i, j), elt, F2);
    }

  fq_nmod_mat_kernel(ker_2, q_2, F2);
  nullity = fq_nmod_mat_nrows(ker_2, F2);

  nmod_mat_init(lift, nullity, QF_RANK, 4);
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < nullity; j++)
      nmod_mat_set_entry(lift, j, i, nmod_poly_get_coeff_ui(fq_nmod_mat_entry(ker_2, j, i),0));

  nmod_mat_init(lift_q, nullity, QF_RANK, 4);
  nmod_mat_mul(lift_q, lift, q);
  fq_nmod_mat_init(lift_2, nullity, QF_RANK, F2);
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < nullity; j++) {
      fq_nmod_set_si(elt, nmod_mat_get_entry(lift_q, j, i) / 2, F2);
      fq_nmod_set(fq_nmod_mat_entry(lift_2, j, i), elt, F2);
    }

  // Does that even work?
  // Compute the Echelon form of the matrix.
  fq_nmod_mat_init(trans, QF_RANK, QF_RANK, F2);

#ifdef DEBUG
  fq_nmod_mat_init_set(test1, q_2, F2);
#endif // DEBUG
  
  fq_nmod_mat_rref_trans(q_2, trans, F2);

#ifdef DEBUG
  fq_nmod_mat_init(test2, QF_RANK, QF_RANK, F2);
  fq_nmod_mat_mul(test2, trans, test1, F2);

  assert(fq_nmod_mat_equal(test2, q_2, F2));
  
  fq_nmod_mat_clear(test1, F2);
  fq_nmod_mat_clear(test2, F2);
#endif // DEBUG

  fq_nmod_mat_init(trans_t, QF_RANK, QF_RANK, F2);
  omf5_fq_nmod_mat_transpose(trans_t, trans, F2);
  fq_nmod_mat_mul(lift_2, lift_2, trans_t, F2);

  fq_nmod_mat_init(sol, nullity, QF_RANK, F2);
  
  // The system in reduced echelon form is solved with the variables of the
  // columns which have no pivot set to 0. Here pivot is the number of pivots
  // in the columns before i, so the entry in the row pivot and the column i
  // is not 0 exactly when it is the pivot of that row.
  pivot = 0;
  for (i = 0; i < QF_RANK; i++) {
    if ((pivot < QF_RANK) && !fq_nmod_is_zero(fq_nmod_mat_entry(q_2, pivot, i), F2)) {
      for (j = 0; j < nullity; j++) {
	fq_nmod_set(fq_nmod_mat_entry(sol, j, i), fq_nmod_mat_entry(lift_2, j, pivot), F2);
      }
      pivot++;
    }
    else {
      for (j = 0; j < nullity; j++) {
	fq_nmod_zero(fq_nmod_mat_entry(sol, j, i), F2);
      }
    }
  }
  
  
  for (i = 0; i < nmod_mat_nrows(q); i++)
    for (j = 0; j < nullity; j++) {
      nmod_mat_set_entry(ker, i, j, nmod_poly_get_coeff_ui(fq_nmod_mat_entry(sol, j, i),0) * 2 +
			 nmod_mat_get_entry(lift, j, i));
    }

#ifdef DEBUG
  nmod_mat_init(zero, QF_RANK, QF_RANK, 4);
  nmod_mat_init(ker_t, QF_RANK, QF_RANK, 4);
  nmod_mat_transpose(ker_t, ker);
  nmod_mat_mul(zero, ker_t, q);
  assert(nmod_mat_is_zero(zero));
  nmod_mat_clear(zero);
  nmod_mat_clear(ker_t);
#endif // DEBUG
  
  fq_nmod_clear(elt, F2);
  fq_nmod_mat_clear(ker_2, F2);
  nmod_mat_clear(lift);
  nmod_mat_clear(lift_q);
  fq_nmod_mat_clear(sol, F2);
  fq_nmod_mat_clear(q_2, F2);
  fq_nmod_mat_clear(trans, F2);
  fq_nmod_mat_clear(trans_t, F2);
  fq_nmod_mat_clear(lift_2, F2);

  fq_nmod_ctx_clear(F2);

  fmpz_clear(two);
  
  return nullity;
}

void spinor_init_fmpz(spinor_t spinor, const fmpz_t disc)
{
  fmpz_factor_t bad_primes;
  slong prime_idx, p;

  fmpz_factor_init(bad_primes);
  fmpz_factor(bad_primes, disc);

  fmpz_mat_init(spinor->Q,0,0);
  spinor->num_primes = bad_primes->num;
  spinor->twist = (1LL << (bad_primes->num)) - 1;
  spinor->rads = (nmod_mat_t*)malloc((bad_primes->num) * sizeof(nmod_mat_t));
  spinor->primes = (nmod_t*)malloc((bad_primes->num) * sizeof(nmod_t));
  spinor->pivots = (slong*)malloc((bad_primes->num) * sizeof(slong));
  spinor->exps = (slong*)malloc((bad_primes->num) * sizeof(slong));
  spinor->lifts = (fmpz_mat_t*)malloc((bad_primes->num) * sizeof(fmpz_mat_t));

  for (prime_idx = 0; prime_idx < bad_primes->num; prime_idx++) {
    p = fmpz_get_si(&(bad_primes->p[prime_idx]));
    if (p == 2)
      p = 4;
    nmod_init(&(spinor->primes[prime_idx]), p);
    nmod_mat_init(spinor->rads[prime_idx], 0, 0, spinor->primes[prime_idx].n);
    fmpz_mat_init(spinor->lifts[prime_idx], 0, 0);
    spinor->exps[prime_idx] = bad_primes->exp[prime_idx];
  }

  fmpz_factor_clear(bad_primes);

  return;
}

// The adjugate of the nonsingular matrix q.
static void adjugate(fmpz_mat_t adj, const fmpz_mat_t q)
{
  fmpq_mat_t q_q, q_inv;
  fmpz_t det;
  int invertible;

  fmpz_init(det);
  fmpz_mat_det(det, q);

  fmpq_mat_init(q_q, fmpz_mat_nrows(q), fmpz_mat_ncols(q));
  fmpq_mat_init(q_inv, fmpz_mat_nrows(q), fmpz_mat_ncols(q));
  fmpq_mat_set_fmpz_mat(q_q, q);
  invertible = fmpq_mat_inv(q_inv, q_q);
  assert(invertible);
  (void)invertible;
  fmpq_mat_scalar_mul_fmpz(q_inv, q_inv, det);
  // the entries are integers
  fmpq_mat_get_fmpz_mat(adj, q_inv);

  fmpq_mat_clear(q_inv);
  fmpq_mat_clear(q_q);
  fmpz_clear(det);

  return;
}

void spinor_init_fmpz_mat(spinor_t spinor, const fmpz_mat_t q)
{
  slong prime_idx;
  fmpz_t det;
  nmod_mat_t q_p;
  fmpz_factor_t bad_primes;
  nmod_mat_t ker;
  fmpz_mat_t adj;
  slong nullity, i, j;
  slong p, pivot;

  fmpz_init(det);
  fmpz_mat_det(det, q);
  fmpz_divexact_si(det, det, 2);
  fmpz_factor_init(bad_primes);
  fmpz_factor(bad_primes, det);  
  fmpz_clear(det);

  fmpz_mat_init(adj, fmpz_mat_nrows(q), fmpz_mat_ncols(q));
  adjugate(adj, q);

  fmpz_mat_init_set(spinor->Q,q);
  spinor->num_primes = bad_primes->num;
  spinor->twist = (1LL << (bad_primes->num)) - 1;
  spinor->rads = (nmod_mat_t*)malloc((bad_primes->num) * sizeof(nmod_mat_t));
  spinor->primes = (nmod_t*)malloc((bad_primes->num) * sizeof(nmod_t));
  spinor->pivots = (slong*)malloc((bad_primes->num) * sizeof(slong));
  spinor->exps = (slong*)malloc((bad_primes->num) * sizeof(slong));
  spinor->lifts = (fmpz_mat_t*)malloc((bad_primes->num) * sizeof(fmpz_mat_t));
  
  for (prime_idx = 0; prime_idx < bad_primes->num; prime_idx++) {
    p = fmpz_get_si(&(bad_primes->p[prime_idx]));
    spinor->exps[prime_idx] = bad_primes->exp[prime_idx];

    // The product of q and its adjugate is det(q), so that any column of the adjugate
    // is orthogonal to the lattice modulo the power of p dividing det(q).
    // When the radical modulo p is a line the adjugate has rank 1 modulo p,
    // and we take a column which is not 0 modulo p. Otherwise it is 0 modulo p,
    // and the vector is left as 0.
    fmpz_mat_init(spinor->lifts[prime_idx], fmpz_mat_nrows(q), 1);
    for (j = 0; j < fmpz_mat_ncols(adj); j++) {
      for (i = 0; i < fmpz_mat_nrows(adj); i++)
	if (fmpz_fdiv_ui(fmpz_mat_entry(adj, i, j), p) != 0)
	  break;
      if (i < fmpz_mat_nrows(adj))
	break;
    }
    if (j < fmpz_mat_ncols(adj))
      for (i = 0; i < fmpz_mat_nrows(adj); i++)
	fmpz_set(fmpz_mat_entry(spinor->lifts[prime_idx], i, 0), fmpz_mat_entry(adj, i, j));

    if (p == 2)
      p = 4;
    nmod_init(&(spinor->primes[prime_idx]),p);
    nmod_mat_init_set_fmpz_mat(q_p, q, p);
    nmod_mat_init(ker, QF_RANK, QF_RANK, p);
    // apparentaly, this doesn't work for p = 4 !!!
    // we go around that
    if (p != 4)
      nullity = nmod_mat_nullspace(ker, q_p);
    else
      nullity = nmod_mat_nullspace_mod4(ker, q_p);
    nmod_mat_init(spinor->rads[prime_idx], nullity, QF_RANK, p);
    for (i = 0; i < QF_RANK; i++)
      for (j = 0; j < nullity; j++)
	nmod_mat_set_entry(spinor->rads[prime_idx], j, i, nmod_mat_get_entry(ker, i, j));

    for (pivot = 0; pivot < nmod_mat_ncols(spinor->rads[prime_idx]);) {
      if (nmod_mat_entry(spinor->rads[prime_idx], 0, pivot) % fmpz_get_si(&(bad_primes->p[prime_idx])) != 0)
	break;
      pivot++;
    }
    spinor->pivots[prime_idx] = pivot;
    
    nmod_mat_clear(ker);
    nmod_mat_clear(q_p);
  }
  fmpz_mat_clear(adj);
  fmpz_factor_clear(bad_primes);
  return;
}

void spinor_init_square_matrix(spinor_t spinor, const square_matrix_t q)
{
  fmpz_mat_t q_fmpz;

  fmpz_mat_init_set_square_matrix(q_fmpz, q);
  spinor_init_fmpz_mat(spinor, q_fmpz);
  fmpz_mat_clear(q_fmpz);
  
  return;
}

void spinor_clear(spinor_t spinor)
{
  slong prime_idx;

  fmpz_mat_clear(spinor->Q);
  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++) {
    nmod_mat_clear(spinor->rads[prime_idx]);
    fmpz_mat_clear(spinor->lifts[prime_idx]);
  }
  free(spinor->rads);
  free(spinor->lifts);
  
  free(spinor->primes);
  free(spinor->pivots);
  free(spinor->exps);

  return;
}

// assumes a is a vector of mod p values?
W64 spinor_compute_vals(const spinor_t spinor, const ulong* a)
{
  W64 val, mask;
  slong i;

  val = 0;
  mask = 1;
  for (i = 0; i < spinor->num_primes; i++) {
    if (a[i] != 1) {
#ifdef DEBUG
      assert(a[i] == spinor->primes[i].n-1);
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
  
  assert(val==spinor_norm_cd_fmpz_mat(spinor, mat_fmpz, denom_fmpz));

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

  val = spinor_norm_fmpz_mat(spinor, mat_fmpz, denom_fmpz);

  // This is not valid when p^2 || n
  //  assert(val==spinor_norm_cd_fmpz_mat(spinor, mat_fmpz, denom_fmpz));

  fmpz_mat_clear(mat_fmpz);
  fmpz_clear(denom_fmpz);
  return val;
}

// The prime of the spinor at the given index. The modulus used for the prime 2 is 4.
static ulong spinor_prime(const spinor_t spinor, slong prime_idx)
{
  return (spinor->primes[prime_idx].n == 4) ? 2 : spinor->primes[prime_idx].n;
}

// Whether one of the primes of the spinor divides denom.
static bool spinor_divides(const spinor_t spinor, const fmpz_t denom)
{
  slong prime_idx;

  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++)
    if (fmpz_fdiv_ui(denom, spinor_prime(spinor, prime_idx)) == 0)
      return true;

  return false;
}

// The primes of the spinor are those dividing the discriminant, and the spinor
// has a character with values in {1,-1} at each of them. At a prime q at which
// the isometry is integral, it is the sign by which the isometry of determinant 1
// acts on the radical of the lattice modulo q (modulo 4 when q = 2).
// At a prime q at which the isometry is not integral:
// - When q divides the discriminant exactly once, it is given by the parity of
//   the valuation at q of the spinor norm, which is a character of all the
//   isometries and agrees with the sign above on the integral ones.
// - When q^2 divides the discriminant, the lattice L and a q-neighbor N of it
//   have the same radical line modulo q: if q^k is the power of q dividing
//   the discriminant and e in L has B(e,L) in q^k Z, then e is in N, it is
//   not in qN, and B(e,N) is in q^(k-1) Z. An isometry taking L to N then acts
//   on this line by a sign, and its products with integral isometries act by
//   the product of the signs. The character is this sign, as it is in
//   RadicalSignCharacter of the package ModFrmAlg of Magma. This is the choice
//   for which the L-function of an eigenform whose conductor is divisible by q
//   satisfies its functional equation with the factor 1 - T_q x + q^3 x^2 at q,
//   in the examples in which it was computed. The parity of the valuation of
//   the spinor norm is not used here: when k is even it does not give the
//   character, and when k is odd it gives the opposite sign.

// The sign by which the isometry mat/denom, in lowest terms, acts on the radical
// line at the prime q of the given index, where q divides denom exactly once and
// q^2 divides the discriminant; det_sign is the determinant of the isometry.
// The sign is found with a vector spanning the radical modulo the power of q
// dividing the discriminant, so that its image can be divided by q.
// Returns 1 or -1, and 0 if the isometry does not act on the line by a sign,
// or if q^2 divides denom, or if q^2 does not divide the discriminant,
// or if the radical modulo q is not a line.
static int spinor_radical_sign(const spinor_t spinor, slong prime_idx,
			       const fmpz_mat_t mat, const fmpz_t denom, int det_sign)
{
  fmpz_mat_t image;
  fmpz_t q, unit;
  nmod_t p;
  slong i, n, val, pivot;
  ulong c, unit_p, image_p, lift_p;
  bool is_sign;

  n = fmpz_mat_nrows(mat);
  p = spinor->primes[prime_idx];

  fmpz_init(unit);
  fmpz_init_set_ui(q, spinor_prime(spinor, prime_idx));
  val = fmpz_remove(unit, denom, q);

  // Another vector spanning the radical changes the image by q^k / denom times
  // a vector of the lattice, where q^k is the power of q dividing the discriminant.
  // (For q = 2 it is 2^(k+1), and the image is needed modulo 4.)
  if ((val != 1) || (spinor->exps[prime_idx] < 2)) {
    fmpz_clear(q);
    fmpz_clear(unit);
    return 0;
  }

  fmpz_mat_init(image, n, 1);
  fmpz_mat_mul(image, mat, spinor->lifts[prime_idx]);

  is_sign = true;
  for (i = 0; (i < n) && is_sign; i++)
    if (fmpz_divisible(fmpz_mat_entry(image, i, 0), q))
      fmpz_divexact(fmpz_mat_entry(image, i, 0), fmpz_mat_entry(image, i, 0), q);
    else
      is_sign = false;

  // an entry of the vector spanning the radical which is not 0 modulo q,
  // if the radical is a line
  for (pivot = 0; pivot < n; pivot++)
    if (!fmpz_divisible(fmpz_mat_entry(spinor->lifts[prime_idx], pivot, 0), q))
      break;
  if (pivot == n)
    is_sign = false;

  c = 0;
  if (is_sign) {
    i = pivot;
    unit_p = fmpz_get_nmod(unit, p);
    image_p = fmpz_get_nmod(fmpz_mat_entry(image, i, 0), p);
    lift_p = fmpz_get_nmod(fmpz_mat_entry(spinor->lifts[prime_idx], i, 0), p);
    c = nmod_div(image_p, nmod_mul(unit_p, lift_p, p), p);
    if (det_sign == -1)
      c = nmod_neg(c, p);
    is_sign = (c == 1) || (c == p.n-1);
    // the image is c times the lift
    for (i = 0; (i < n) && is_sign; i++) {
      image_p = fmpz_get_nmod(fmpz_mat_entry(image, i, 0), p);
      if (det_sign == -1)
	image_p = nmod_neg(image_p, p);
      lift_p = fmpz_get_nmod(fmpz_mat_entry(spinor->lifts[prime_idx], i, 0), p);
      if (image_p != nmod_mul(c, nmod_mul(unit_p, lift_p, p), p))
	is_sign = false;
    }
  }

  fmpz_mat_clear(image);
  fmpz_clear(q);
  fmpz_clear(unit);

  if (!is_sign)
    return 0;

  return (c == 1) ? 1 : -1;
}

// TODO - use the modular arithmetic built in inside nmod_t instead of gcdext
// Here mat/denom is the isometry in lowest terms,
// and we're computing the spinor norm at all spinor primes
static W64 spinor_norm_lowest_terms(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom)
{
  nmod_mat_t mat_p, mat_p_t, rad_mat;
  slong idx, row, col, pivot, prime_idx, n;
  fmpz_t det; // for some reason det is not implemented for fq_nmod_mat. we reduce the det mod p instead
  ulong denom_p, det_p;
  ulong* evs;
  W64 norm, zas;
  bool has_zas;
  int sign;
  nmod_t p;

  n = fmpz_mat_nrows(mat);
  assert(n == fmpz_mat_ncols(mat));
  
  evs = (ulong*)malloc((spinor->num_primes)*sizeof(ulong));

  zas = 0;
  has_zas = false;

  fmpz_init(det);
  fmpz_mat_det(det, mat);
  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++) {
    p = spinor->primes[prime_idx];
    if (fmpz_fdiv_ui(denom, spinor_prime(spinor, prime_idx)) == 0) {
      // The isometry is not integral at this prime.
      evs[prime_idx] = 1;
      if (spinor->exps[prime_idx] == 1) {
	if (!has_zas) {
	  zas = spinor_norm_zas_fmpz_mat(spinor, mat, denom);
	  has_zas = true;
	}
	if ((zas >> prime_idx) & 1)
	  evs[prime_idx] = p.n-1;
      }
      else {
	sign = spinor_radical_sign(spinor, prime_idx, mat, denom,
				   fmpz_sgn(det) * fmpz_sgn(denom));
	if (sign == 0) {
	  if (nmod_mat_nrows(spinor->rads[prime_idx]) != 1)
	    fprintf(stderr, "Error: the radical of the lattice modulo %lu is not a line, and the Hecke operator at %lu is not supported for such a lattice\n",
		    spinor_prime(spinor, prime_idx), spinor_prime(spinor, prime_idx));
	  else
	    fprintf(stderr, "Error: an isometry does not act by a sign on the radical of the lattice at the prime %lu, which divides its denominator and whose square divides the discriminant\n",
		    spinor_prime(spinor, prime_idx));
	  exit(EXIT_FAILURE);
	}
	if (sign == -1)
	  evs[prime_idx] = p.n-1;
      }
      continue;
    }
    // The denominator is invertible, also modulo 4 when the prime is 2.
    denom_p = fmpz_get_nmod(denom, p);
    
    assert(denom_p < p.n);

    nmod_mat_init_set_fmpz_mat(mat_p, mat, p.n);
    nmod_mat_init(mat_p_t, n, n, p.n);
      
    // scaling the reduced matrix - should be replaced bu a one-line function
    for (row = 0; row < nmod_mat_nrows(mat_p); row++)
      for (col = 0; col < nmod_mat_ncols(mat_p); col++)
	nmod_mat_set_entry(mat_p, row, col, nmod_div(nmod_mat_entry(mat_p, row, col), denom_p, p));

    // !! TODO - we could instead multiply the end result by the determinant?
    det_p = fmpz_get_nmod(det, p);  
    // scaling the determinant
    for (idx = 0; idx < n; idx++)
      det_p = nmod_div(det_p, denom_p, p);
    // If the matrix has determinant -1, we modify it to have a matrix in SO
    if (det_p != 1) {
      assert(det_p == p.n-1);
      nmod_mat_neg(mat_p, mat_p);
    }
    
#ifdef DEBUG_LEVEL_FULL
    printf("rad = \n");
    nmod_mat_print_pretty(spinor->rads[prime_idx]);
    printf("\n");
    printf("s_p =  \n");
    nmod_mat_print_pretty(mat_p);
    printf("\n");
    printf("scale (inverse) = %lu (mod %lu)\n", denom_p, p.n);
#endif // DEBUG_LEVEL_FULL
    
    nmod_mat_init(rad_mat, nmod_mat_nrows(spinor->rads[prime_idx]),
		  nmod_mat_ncols(mat_p), p.n);
    // we need to transpose, given the direction of our isometries
    nmod_mat_transpose(mat_p_t, mat_p);
    nmod_mat_mul(rad_mat, spinor->rads[prime_idx], mat_p_t);

    // for now assume rad is a single vector (we choose our lattices this way)
    // !! TODO !! -  modify to determinant so it will work in the general case
    assert(nmod_mat_nrows(spinor->rads[prime_idx]) == 1);

    pivot = spinor->pivots[prime_idx];
    
    // !! TODO !! - we should always set the pivot in the kernel to be 1
    // so that no arithmetic will be needed here.
    evs[prime_idx] = nmod_div(nmod_mat_entry(rad_mat, 0, pivot),
			      nmod_mat_entry(spinor->rads[prime_idx], 0, pivot), p);

    assert((evs[prime_idx] == 1) || (evs[prime_idx] == p.n-1));

    nmod_mat_clear(rad_mat);
    nmod_mat_clear(mat_p_t);
    nmod_mat_clear(mat_p);
  }
  fmpz_clear(det);

  norm = spinor_compute_vals(spinor, evs);

  free(evs);
  
  return norm;
}

// Here mat/denom is the isometry, and we're computing the spinor norm at all spinor primes
W64 spinor_norm_fmpz_mat(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom)
{
  fmpz_mat_t mat_red;
  fmpz_t denom_red, g;
  W64 norm;

  // The fraction is not always in lowest terms, and whether the isometry is integral
  // at a prime is read from the denominator. This is only needed when
  // the denominator is divisible by a prime of the spinor.
  if (!spinor_divides(spinor, denom))
    return spinor_norm_lowest_terms(spinor, mat, denom);

  fmpz_init(g);
  fmpz_mat_content(g, mat);
  fmpz_gcd(g, g, denom);

  fmpz_mat_init(mat_red, fmpz_mat_nrows(mat), fmpz_mat_ncols(mat));
  fmpz_init(denom_red);
  fmpz_mat_scalar_divexact_fmpz(mat_red, mat, g);
  fmpz_divexact(denom_red, denom, g);

  norm = spinor_norm_lowest_terms(spinor, mat_red, denom_red);

  fmpz_clear(denom_red);
  fmpz_mat_clear(mat_red);
  fmpz_clear(g);

  return norm;
}

// Here below we write the spinor norm according to the Cartan-Dieudonne decomposition

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

// this is zassenhaus' method of computing the spinor norm
W64 spinor_norm_zas_fmpz_mat(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom)
{
  fmpq_mat_t mat_q, Q_q, one;
  fmpq_t det;
  fmpq_t spinorm;
  slong prime_idx, n;
  fmpz_t p;
  
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
  fmpq_mat_init(one, n, n);
  fmpq_mat_one(one);
  fmpq_mat_add(mat_q, mat_q, one);
  fmpq_mat_clear(one);

  fmpq_init(spinorm);
  fmpq_mat_det(spinorm, mat_q);
  fmpq_mul_si(spinorm, spinorm, 2);
  if (fmpq_is_zero(spinorm)) {
    fmpq_clear(spinorm);
    fmpq_clear(det);
    fmpq_mat_clear(mat_q);
    fmpq_mat_clear(Q_q);
    return spinor_norm_cd_fmpz_mat(spinor, mat, denom);
  }
  val = 0;
  fmpz_init(p);
  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++) {
    fmpz_set_ui(p, spinor_prime(spinor, prime_idx));
    val ^= ((W64)(fmpq_valuation(spinorm, p) & 1)) << prime_idx;
  }
  fmpz_clear(p);
  fmpq_clear(spinorm);
  fmpq_clear(det);
  fmpq_mat_clear(mat_q);
  fmpq_mat_clear(Q_q);

  return val;
}

W64 spinor_norm_cd_fmpz_mat(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom)
{
  fmpq_mat_t mat_q, Q_q;
  fmpq_t det;
  fmpq_t spinorm;
  slong prime_idx, n;
  fmpz_t p;
  
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
  fmpz_init(p);
  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++) {
    fmpz_set_ui(p, spinor_prime(spinor, prime_idx));
    val ^= ((W64)(fmpq_valuation(spinorm, p) & 1)) << prime_idx;
  }
  fmpz_clear(p);
  fmpq_clear(spinorm);
  fmpq_clear(det);
  fmpq_mat_clear(mat_q);
  fmpq_mat_clear(Q_q);

  return val;
}
