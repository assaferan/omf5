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
  // The vectors found are in the kernel modulo 4 when the kernel modulo 2 is a
  // line. Otherwise the system modulo 4 may have no solution, and the caller
  // then uses the dimension only.
  if (nullity == 1) {
    nmod_mat_init(zero, QF_RANK, QF_RANK, 4);
    nmod_mat_init(ker_t, QF_RANK, QF_RANK, 4);
    nmod_mat_transpose(ker_t, ker);
    nmod_mat_mul(zero, ker_t, q);
    assert(nmod_mat_is_zero(zero));
    nmod_mat_clear(zero);
    nmod_mat_clear(ker_t);
  }
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
  spinor->odd_on_rad_2 = false;
  fmpz_mat_init(spinor->deep_2, 0, 0);
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

// The rows of rad are a basis of the radical modulo 2 of the lattice whose
// Gram matrix is q. Returns whether the quadratic form is odd at a vector of
// the radical. The parity of the form at v depends only on v modulo 2 when v
// is in the radical modulo 2, and the form modulo 2 is additive on the
// radical, so that it is enough to look at the vectors of the basis.
static bool spinor_radical_2_is_odd(const nmod_mat_t rad, const fmpz_mat_t q)
{
  nmod_mat_t q_4, rad_t, rad_q, prod;
  slong i, n, r;
  bool is_odd;

  n = nmod_mat_ncols(rad);
  r = nmod_mat_nrows(rad);

  // The value of the quadratic form at v is half of that of the Gram matrix.
  nmod_mat_init_set_fmpz_mat(q_4, q, 4);
  nmod_mat_init(rad_t, n, r, 4);
  nmod_mat_init(rad_q, r, n, 4);
  nmod_mat_init(prod, r, r, 4);
  nmod_mat_transpose(rad_t, rad);
  nmod_mat_mul(rad_q, rad, q_4);
  nmod_mat_mul(prod, rad_q, rad_t);

  is_odd = false;
  for (i = 0; i < r; i++)
    if (nmod_mat_entry(prod, i, i) != 0)
      is_odd = true;

  nmod_mat_clear(prod);
  nmod_mat_clear(rad_q);
  nmod_mat_clear(rad_t);
  nmod_mat_clear(q_4);

  return is_odd;
}

// Sets the columns of deep to vectors which span, modulo 2, the deepest part
// of the radical at 2 of the lattice whose Gram matrix has the adjugate adj,
// when it is a line or a plane; deep has no columns otherwise.
// A vector v is in the radical modulo 2^k when its products with the lattice
// are divisible by 2^k, and the deepest part is the image modulo 2 of the
// radical modulo 2^k for the largest k at which this image is not 0. It is the
// last component of a Jordan decomposition at 2, modulo 2.
// The product of the Gram matrix and its adjugate is the determinant, so the
// columns of the adjugate divided by the power of 2 in its content are in this
// radical, and they span its image modulo 2.
static void spinor_deep_part_2(fmpz_mat_t deep, const fmpz_mat_t adj)
{
  fmpz_mat_t a;
  nmod_mat_t a_2, cols;
  fmpz_t content, unit, two;
  slong i, j, k, n, rank;

  n = fmpz_mat_nrows(adj);
  fmpz_init(content);
  fmpz_init(unit);
  fmpz_init_set_ui(two, 2);
  fmpz_mat_init(a, n, n);

  fmpz_mat_content(content, adj);
  fmpz_remove(unit, content, two);
  fmpz_divexact(content, content, unit);
  fmpz_mat_scalar_divexact_fmpz(a, adj, content);
  nmod_mat_init_set_fmpz_mat(a_2, a, 2);
  rank = nmod_mat_rank(a_2);
  if (rank > 2)
    rank = 0;

  fmpz_mat_clear(deep);
  fmpz_mat_init(deep, n, rank);
  // columns of the matrix which are linearly independent modulo 2
  nmod_mat_init(cols, n, rank, 2);
  k = 0;
  for (j = 0; (j < n) && (k < rank); j++) {
    for (i = 0; i < n; i++)
      nmod_mat_set_entry(cols, i, k, nmod_mat_entry(a_2, i, j));
    if (nmod_mat_rank(cols) == k + 1) {
      for (i = 0; i < n; i++)
	fmpz_set(fmpz_mat_entry(deep, i, k), fmpz_mat_entry(a, i, j));
      k++;
    }
  }

  nmod_mat_clear(cols);
  nmod_mat_clear(a_2);
  fmpz_mat_clear(a);
  fmpz_clear(two);
  fmpz_clear(unit);
  fmpz_clear(content);

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
  slong nullity, i, j, factor_idx;
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
  spinor->odd_on_rad_2 = false;
  fmpz_mat_init(spinor->deep_2, fmpz_mat_nrows(q), 0);
  spinor->rads = (nmod_mat_t*)malloc((bad_primes->num) * sizeof(nmod_mat_t));
  spinor->primes = (nmod_t*)malloc((bad_primes->num) * sizeof(nmod_t));
  spinor->pivots = (slong*)malloc((bad_primes->num) * sizeof(slong));
  spinor->exps = (slong*)malloc((bad_primes->num) * sizeof(slong));
  spinor->lifts = (fmpz_mat_t*)malloc((bad_primes->num) * sizeof(fmpz_mat_t));
  
  // The primes of the spinor are those dividing the discriminant, except for
  // the prime 2 when the radical modulo 2 is the whole lattice modulo 2.
  prime_idx = 0;
  for (factor_idx = 0; factor_idx < bad_primes->num; factor_idx++) {
    p = fmpz_get_si(&(bad_primes->p[factor_idx]));
    spinor->exps[prime_idx] = bad_primes->exp[factor_idx];

    // The product of q and its adjugate is det(q), so that any column of the adjugate
    // is orthogonal to the lattice modulo the power of p dividing det(q).
    // When the radical modulo p is a line the adjugate has rank 1 modulo p,
    // and we take a column which is not 0 modulo p. Otherwise it is 0 modulo p,
    // and the vector is left as 0. It is used only when the radical is a line.
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

    nmod_mat_clear(ker);
    nmod_mat_clear(q_p);

    if (p == 4) {
      // The radical modulo 2 is of odd dimension. When it is a line, the
      // character at 2 is the sign by which an isometry acts on the radical
      // modulo 4. When it is of dimension 3, only the radical modulo 2 is
      // used. If the quadratic form is even on the radical, the character is
      // the sign of a permutation of the quotient of the lattice modulo 2 by
      // the radical, a plane over the field of 2 elements. If the form is odd
      // at a vector of the radical, the character is that of the spinor norm
      // when 16 does not divide the discriminant. When 16 divides it, the
      // deepest part of the radical at 2 is a line or a plane, and the
      // character is the sign by which an isometry acts on it.
      // Otherwise the radical is the whole lattice modulo 2, there is no
      // character at 2 here, and the conductors are odd.
      if (nullity == 3) {
	for (j = 0; j < nullity; j++)
	  for (i = 0; i < QF_RANK; i++)
	    nmod_mat_set_entry(spinor->rads[prime_idx], j, i,
			       nmod_mat_entry(spinor->rads[prime_idx], j, i) % 2);
	spinor->odd_on_rad_2 = spinor_radical_2_is_odd(spinor->rads[prime_idx], q);
	if (spinor->odd_on_rad_2 && (spinor->exps[prime_idx] >= 4))
	  spinor_deep_part_2(spinor->deep_2, adj);
      }
      else if (nullity != 1) {
	nmod_mat_clear(spinor->rads[prime_idx]);
	fmpz_mat_clear(spinor->lifts[prime_idx]);
	continue;
      }
    }
    else
      // The basis of the radical is in reduced row echelon form, so that the
      // coordinates of a vector of the radical are its entries at the pivots.
      nmod_mat_rref(spinor->rads[prime_idx]);

    for (pivot = 0; pivot < nmod_mat_ncols(spinor->rads[prime_idx]);) {
      if (nmod_mat_entry(spinor->rads[prime_idx], 0, pivot) % fmpz_get_si(&(bad_primes->p[factor_idx])) != 0)
	break;
      pivot++;
    }
    spinor->pivots[prime_idx] = pivot;

    prime_idx++;
  }
  spinor->num_primes = prime_idx;
  spinor->twist = (1LL << prime_idx) - 1;
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
  fmpz_mat_clear(spinor->deep_2);
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

// Whether p is one of the primes of the spinor.
bool spinor_has_prime(const spinor_t spinor, ulong p)
{
  slong prime_idx;

  for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++)
    if (spinor_prime(spinor, prime_idx) == p)
      return true;

  return false;
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

// The spinor has a character with values in {1,-1} at each of its primes, which
// are those dividing the discriminant. Let q be one of them, L the lattice,
// and R(L) the radical of L modulo q, which is the image in L/qL of the
// vectors e in L with B(e,L) in qZ.
//
// At an isometry of determinant 1 which is integral at q, the character is the
// determinant of the isometry on R(L). When R(L) is a line this is the sign
// by which it acts on the line. For q = 2 the sign is read modulo 4, and this
// is the character when R(L) is a line. The case q = 2 with R(L) of dimension
// 3 is described after the case of a line.
//
// At an isometry g which is not integral at q:
// - When q divides the discriminant exactly once, the character is given by
//   the parity of the valuation at q of the spinor norm, which is a character
//   of all the isometries and agrees with the sign above on the integral ones.
// - When q^2 divides the discriminant, g takes L to a q-neighbor N of it:
//   there is x in L, not in the radical, such that N is spanned by x/q and by
//   the vectors of L orthogonal to x modulo q. The radicals of L and of N are
//   identified as follows. If y in L has B(x,y) = 1 modulo q and e in L has
//   B(e,L) in qZ, then w = e - B(e,x) y is e modulo qL, it is in N, and B(w,N)
//   is in qZ. The class of w modulo qN depends only on the class of e modulo
//   qL, and this is an isomorphism of R(L) with R(N), which does not depend on
//   x or on y. The isometry g takes R(L) to R(N), and the character at g is the
//   determinant of this map followed by the inverse of the identification.
//   It is multiplied by the characters of integral isometries when g is
//   multiplied by them on either side, and it is 1 at an isometry which is the
//   identity on the orthogonal complement of a hyperbolic plane.
//   When R(L) is a line, let q^k be the power of q dividing the discriminant
//   and let e in L have B(e,L) in q^k Z, so that e spans the radical of L
//   modulo q^k. Then e is in N, it is not in qN, and the radical of N modulo
//   q^k is spanned by a vector which is e modulo q^(k-1) L. So g takes e to
//   c e modulo q^(k-1) L for a sign c, and the character at g is c, as it is
//   in RadicalSignCharacter of the package ModFrmAlg of Magma. This is the
//   choice for which the L-function of an eigenform whose conductor is
//   divisible by q satisfies its functional equation with the factor
//   1 - T_q x + q^3 x^2 at q, in the examples in which it was computed.
//   The parity of the valuation of the spinor norm is not used here: when k
//   is even it does not give the character, and when k is odd it gives the
//   opposite sign.
// - The Hecke operator T_{q^2} has the neighbors N which contain X/q for a
//   plane X in L, isotropic modulo q, which does not meet the radical modulo
//   q; then qL is in N and qN is in L. When R(L) is a line the sign c is
//   defined for the isometries g to these neighbors in the same way, and it
//   is the character at g. When R(L) is not a line and q is odd there are no
//   such neighbors: L modulo its radical has rank at most 3, and no isotropic
//   plane.
//
// Let q = 2. The radical R(L) of the bilinear form modulo 2 is of odd
// dimension. When it is a line, an isometry g integral at 2 acts on L/R(L),
// a quadratic space of dimension 4 over the field of 2 elements, and its sign
// on the radical modulo 4 is (-1)^rank(g - 1) there, the Dickson invariant.
// When R(L) is of dimension 3, the quadratic form Q modulo 2 is additive on
// R(L), and the character depends on whether it vanishes there.
// - Q vanishes on R(L) modulo 2. The character at g is (-1)^rank(g - 1) for
//   the action of g on L/R(L), a plane over the field of 2 elements, which is
//   the sign of the permutation of its three nonzero vectors. Then Q modulo 2
//   is defined on L/R(L), which is a hyperbolic plane when L has 2-neighbors,
//   and has two isotropic vectors other than 0. An isometry integral at 2
//   fixes both or exchanges them, and the character is 1 or -1 accordingly.
//   Let g be an isometry to a 2-neighbor N of L, where N is the 2-neighbor at
//   x and g^(-1) L the 2-neighbor at x'. The classes of x and x' in L/R(L) are
//   isotropic and not 0, and the character at g is -1 when they are the same
//   and 1 when they are not. It is multiplied by the characters of integral
//   isometries when g is multiplied by them on either side: on the left they
//   move the class of x, and on the right that of x'. It is 1 at an isometry
//   which is the identity on the orthogonal complement of a hyperbolic plane
//   spanned by x and y, and takes x to x/2 and y to 2y: there x' = y.
// - Q does not vanish on R(L) modulo 2. The character at g is 1 when the
//   valuation at 2 of the spinor norm of g is even, and -1 when it is odd, as
//   when 2 divides the discriminant exactly once. At an isometry to a
//   2-neighbor it is the opposite, so that it is 1 at an isometry which is
//   the identity on the orthogonal complement of a hyperbolic plane, whose
//   spinor norm is 2.
//   When 4 divides the discriminant exactly, this is again the sign of a
//   permutation of a plane, in the lattices in which the two were compared,
//   and at an isometry integral at 2 this sign is used, because it is much
//   cheaper than the spinor norm (spinor_radical_rank_sign_2).
//   The plane is R0(L), the vectors of R(L) at which Q is even, and Q/2 is a
//   quadratic form modulo 2 on it, which is not degenerate. The character at
//   g is (-1)^rank(g - 1) for the action of g on R0(L). At an isometry g to
//   the 2-neighbor N at x, the planes R0(L) and R0(N) are identified as the
//   radicals are for an odd prime: if y in L has B(x,y) odd and e in L has its
//   class in R0(L), then w = e - B(e,x) y is e modulo 2L, it is in N, B(w,N)
//   is in 2Z and Q(w) is even. The character at g is (-1)^rank(h - 1) for the
//   map h of R0(L) which is g followed by the inverse of the identification.
//   When 8 divides the discriminant the action on R0(L) is not enough: when
//   8 divides it exactly, Q/2 is additive on R0(L), it tells the three
//   nonzero vectors apart together with the radical of the bilinear form
//   modulo 4, and every g is the identity on R0(L).
//   The action on L/R(L) is not used in this case. Let x be isotropic in L
//   over the 2-adic integers with N the 2-neighbor at x, let y be isotropic
//   with B(x,y) = 1, and let w be orthogonal to both, with Q(w) odd; its
//   class is in R(L). The isometry
//   v -> v + B(v,x) w - B(v,w) x - Q(w) B(v,x) x
//   takes L to L and N to N. It takes y to y + w - Q(w) x, so it acts on
//   L/R(L) by a transposition, and it is the identity on N/R(N), which is
//   spanned by x/2 and 2y. So a function on the isometries to the 2-neighbors
//   which is multiplied on both sides by the sign of the permutation of
//   L/R(L) is 0.
// When R(L) is of dimension 3 no plane in L modulo 2 which is isotropic
// misses the radical, and T_4 is 0. When R(L) is the whole of L modulo 2
// there is no character at 2 here.

// The coordinates of the rows of vecs in the basis of the radical at the odd
// prime of the given index. Returns false if a row is not in the radical.
static bool spinor_radical_coordinates(nmod_mat_t coords, const spinor_t spinor,
				       slong prime_idx, const nmod_mat_t vecs)
{
  nmod_mat_t prod;
  slong i, j, pivot, n, r;
  bool in_radical;

  n = nmod_mat_ncols(spinor->rads[prime_idx]);
  r = nmod_mat_nrows(spinor->rads[prime_idx]);

  // The basis is in reduced row echelon form.
  for (j = 0; j < r; j++) {
    for (pivot = 0; pivot < n; pivot++)
      if (nmod_mat_entry(spinor->rads[prime_idx], j, pivot) != 0)
	break;
    assert(pivot < n);
    for (i = 0; i < nmod_mat_nrows(vecs); i++)
      nmod_mat_set_entry(coords, i, j, nmod_mat_entry(vecs, i, pivot));
  }

  nmod_mat_init(prod, nmod_mat_nrows(vecs), n, spinor->primes[prime_idx].n);
  nmod_mat_mul(prod, coords, spinor->rads[prime_idx]);
  in_radical = nmod_mat_equal(prod, vecs);
  nmod_mat_clear(prod);

  return in_radical;
}

// The character at the odd prime q of the given index at the isometry
// mat/denom, in lowest terms, where q divides denom exactly once, and the
// isometry takes the lattice to a q-neighbor of it; det_sign is the
// determinant of the isometry. Returns 1 or -1, and 0 if q^2 divides denom
// or if the isometry is not one to a q-neighbor at a vector outside the radical.
// It is used when the radical modulo q is not a line.
static int spinor_neighbor_sign(const spinor_t spinor, slong prime_idx,
				const fmpz_mat_t mat, const fmpz_t denom, int det_sign)
{
  fmpz_mat_t x, gx, inv, w, image;
  fmpz_t q, unit, inv_denom, inv_unit, b, q_pow;
  nmod_mat_t images, coords;
  nmod_t p;
  slong i, j, k, col, n, r, val, inv_val;
  ulong c, scale, det;
  int sign;
  bool ok;

  n = fmpz_mat_nrows(mat);
  p = spinor->primes[prime_idx];
  r = nmod_mat_nrows(spinor->rads[prime_idx]);

  fmpz_init(unit);
  fmpz_init_set_ui(q, p.n);
  val = fmpz_remove(unit, denom, q);
  ok = (val == 1);

  // The neighbor is spanned by the columns of mat/denom. It is the q-neighbor
  // at x for any column x of mat which is not 0 modulo q.
  fmpz_mat_init(x, n, 1);
  fmpz_mat_init(gx, n, 1);
  for (col = 0; col < n; col++) {
    for (i = 0; i < n; i++)
      if (!fmpz_divisible(fmpz_mat_entry(mat, i, col), q))
	break;
    if (i < n)
      break;
  }
  if (col == n)
    ok = false;

  // y is c times the vector of the basis at the index k
  c = 0;
  k = 0;
  if (ok) {
    for (i = 0; i < n; i++)
      fmpz_set(fmpz_mat_entry(x, i, 0), fmpz_mat_entry(mat, i, col));
    fmpz_mat_mul(gx, spinor->Q, x);
    for (k = 0; k < n; k++)
      if (!fmpz_divisible(fmpz_mat_entry(gx, k, 0), q))
	break;
    if (k == n)
      ok = false; // x is in the radical
    else
      c = nmod_inv(fmpz_get_nmod(fmpz_mat_entry(gx, k, 0), p), p);
  }

  // The inverse of the isometry is denom * inv / inv_denom. It takes the
  // radical of the neighbor to that of the lattice.
  fmpz_mat_init(inv, n, n);
  fmpz_init(inv_denom);
  fmpz_init(inv_unit);
  fmpz_init(q_pow);
  inv_val = 0;
  if (ok)
    ok = fmpz_mat_inv(inv, inv_denom, mat);
  if (ok) {
    inv_val = fmpz_remove(inv_unit, inv_denom, q);
    ok = (inv_val >= 1);
  }

  fmpz_init(b);
  fmpz_mat_init(w, n, 1);
  fmpz_mat_init(image, n, 1);
  nmod_mat_init(images, r, n, p.n);
  if (ok) {
    fmpz_pow_ui(q_pow, q, inv_val - 1);
    scale = nmod_div(fmpz_get_nmod(unit, p), fmpz_get_nmod(inv_unit, p), p);
    for (j = 0; (j < r) && ok; j++) {
      // w = e - B(e,x) y for the vector e of the basis of the radical
      for (i = 0; i < n; i++)
	fmpz_set_ui(fmpz_mat_entry(w, i, 0), nmod_mat_entry(spinor->rads[prime_idx], j, i));
      fmpz_zero(b);
      for (i = 0; i < n; i++)
	fmpz_addmul(b, fmpz_mat_entry(w, i, 0), fmpz_mat_entry(gx, i, 0));
      fmpz_mul_ui(b, b, c);
      fmpz_sub(fmpz_mat_entry(w, k, 0), fmpz_mat_entry(w, k, 0), b);
      // its image under the inverse of the isometry, which is in the lattice
      fmpz_mat_mul(image, inv, w);
      for (i = 0; (i < n) && ok; i++)
	if (fmpz_divisible(fmpz_mat_entry(image, i, 0), q_pow)) {
	  fmpz_divexact(fmpz_mat_entry(image, i, 0), fmpz_mat_entry(image, i, 0), q_pow);
	  nmod_mat_set_entry(images, j, i,
			     nmod_mul(fmpz_get_nmod(fmpz_mat_entry(image, i, 0), p), scale, p));
	}
	else
	  ok = false;
    }
  }

  // The determinant is that of the inverse, which is the same sign.
  sign = 0;
  if (ok) {
    nmod_mat_init(coords, r, r, p.n);
    if (spinor_radical_coordinates(coords, spinor, prime_idx, images)) {
      det = nmod_mat_det(coords);
      // the opposite of an isometry of determinant -1 is of determinant 1
      if ((det_sign == -1) && (r % 2 == 1))
	det = nmod_neg(det, p);
      if (det == 1)
	sign = 1;
      else if (det == p.n-1)
	sign = -1;
    }
    nmod_mat_clear(coords);
  }

  nmod_mat_clear(images);
  fmpz_mat_clear(image);
  fmpz_mat_clear(w);
  fmpz_clear(b);
  fmpz_clear(q_pow);
  fmpz_clear(inv_unit);
  fmpz_clear(inv_denom);
  fmpz_mat_clear(inv);
  fmpz_mat_clear(gx);
  fmpz_mat_clear(x);
  fmpz_clear(q);
  fmpz_clear(unit);

  return sign;
}

// The sign by which the isometry mat/denom, in lowest terms, acts on the radical
// line at the prime q of the given index, where q divides denom exactly once and
// q^2 divides the discriminant; det_sign is the determinant of the isometry.
// It is used when the radical modulo q is a line, for the isometries to the
// neighbors of the Hecke operators T_q and T_{q^2}. For q = 2 the sign is
// needed modulo 4.
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

// The character at 2 at an isometry which is integral at 2, when the radical
// modulo 2 is of dimension 3; mat_p is the isometry modulo 4, or its negative.
// It is (-1)^rank(g - 1) for the action of the isometry g on the quotient of
// the lattice modulo 2 by the radical. It is used when the quadratic form is
// even on the radical. Returns 1 or -1.
static int spinor_quotient_sign_2(const spinor_t spinor, slong prime_idx, const nmod_mat_t mat_p)
{
  nmod_mat_t ext;
  slong i, j, n, r, rank;

  n = nmod_mat_nrows(mat_p);
  r = nmod_mat_nrows(spinor->rads[prime_idx]);

  // The columns are those of g - 1 and the vectors of the basis of the radical,
  // and they span the preimage of the image of g - 1 on the quotient.
  nmod_mat_init(ext, n, n + r, 2);
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      nmod_mat_set_entry(ext, i, j, (nmod_mat_entry(mat_p, i, j) + ((i == j) ? 1 : 0)) % 2);
    for (j = 0; j < r; j++)
      nmod_mat_set_entry(ext, i, n + j, nmod_mat_entry(spinor->rads[prime_idx], j, i) % 2);
  }
  rank = nmod_mat_rank(ext) - r;
  nmod_mat_clear(ext);

  return (rank % 2 == 0) ? 1 : -1;
}

// The character at 2 at an isometry which is integral at 2, when the radical
// modulo 2 is of dimension 3, the quadratic form is odd at a vector of it and
// 4 divides the discriminant exactly. The rows of rad_mat are the images of
// the vectors of the basis of the radical, up to a sign.
// It is (-1)^rank(g - 1) for the action of the isometry g on the radical
// modulo 2. Here the Gram matrix on the radical is twice an odd matrix S of
// odd determinant. The vectors of the radical at which the quadratic form is
// even are a plane modulo 2, the radical is the sum of this plane and of the
// line orthogonal to it for S, and g fixes the vector of this line, so the
// character is also the sign of the permutation of the three nonzero vectors
// of the plane. In the lattices in which it was compared, this is the
// character of the spinor norm: it is 1 when the valuation at 2 of the spinor
// norm is even. Returns 1 or -1.
static int spinor_radical_rank_sign_2(const spinor_t spinor, slong prime_idx, const nmod_mat_t rad_mat)
{
  nmod_mat_t diff;
  slong i, j, rank;

  nmod_mat_init(diff, nmod_mat_nrows(rad_mat), nmod_mat_ncols(rad_mat), 2);
  for (i = 0; i < nmod_mat_nrows(rad_mat); i++)
    for (j = 0; j < nmod_mat_ncols(rad_mat); j++)
      nmod_mat_set_entry(diff, i, j, (nmod_mat_entry(rad_mat, i, j) +
				      nmod_mat_entry(spinor->rads[prime_idx], i, j)) % 2);
  rank = nmod_mat_rank(diff);
  nmod_mat_clear(diff);

  return (rank % 2 == 0) ? 1 : -1;
}

// Sets x to a column of mat modulo 2 which is not 0, when mat modulo 2 is of
// rank 1 and x is not in the radical modulo 2 of the Gram matrix q_2.
// Returns whether this is so. x is a column vector modulo 2.
static bool spinor_neighbor_vector_2(nmod_mat_t x, const fmpz_mat_t mat, const nmod_mat_t q_2)
{
  nmod_mat_t mat_2, qx;
  slong i, col, n;
  bool ok;

  n = fmpz_mat_nrows(mat);
  nmod_mat_init_set_fmpz_mat(mat_2, mat, 2);
  ok = (nmod_mat_rank(mat_2) == 1);

  if (ok) {
    for (col = 0; col < n; col++) {
      for (i = 0; i < n; i++)
	if (nmod_mat_entry(mat_2, i, col) != 0)
	  break;
      if (i < n)
	break;
    }
    for (i = 0; i < n; i++)
      nmod_mat_set_entry(x, i, 0, nmod_mat_entry(mat_2, i, col));
    nmod_mat_init(qx, n, 1, 2);
    nmod_mat_mul(qx, q_2, x);
    ok = !nmod_mat_is_zero(qx);
    nmod_mat_clear(qx);
  }

  nmod_mat_clear(mat_2);

  return ok;
}

// The character at 2 at the isometry mat/denom, in lowest terms, where 2
// divides denom exactly once and the isometry takes the lattice to a
// 2-neighbor of it, when the radical modulo 2 is of dimension 3 and the
// quadratic form vanishes on it modulo 2.
// The lattice goes to its 2-neighbor at x under the isometry, and to its
// 2-neighbor at x' under the inverse. The character is -1 when x and x' are
// the same modulo the radical, and 1 otherwise.
// Returns 1 or -1, and 0 if 4 divides denom or if the isometry or its inverse
// is not one to a 2-neighbor at a vector outside the radical.
static int spinor_neighbor_sign_2(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom)
{
  fmpz_mat_t inv;
  fmpz_t inv_denom, content, unit, two;
  nmod_mat_t q_2, x, x_inv, diff;
  slong n, val;
  int sign;
  bool ok;

  n = fmpz_mat_nrows(mat);

  fmpz_init(unit);
  fmpz_init_set_ui(two, 2);
  ok = (fmpz_remove(unit, denom, two) == 1);

  // The inverse of the isometry is denom * inv / inv_denom, and inv / content
  // is its numerator in lowest terms. The power of 2 in its denominator is
  // that of inv_denom / (denom * content).
  fmpz_mat_init(inv, n, n);
  fmpz_init(inv_denom);
  fmpz_init(content);
  if (ok)
    ok = fmpz_mat_inv(inv, inv_denom, mat);
  if (ok) {
    fmpz_mat_content(content, inv);
    fmpz_mat_scalar_divexact_fmpz(inv, inv, content);
    val = fmpz_remove(unit, inv_denom, two);
    val -= fmpz_remove(unit, content, two);
    ok = (val == 2);
  }

  nmod_mat_init_set_fmpz_mat(q_2, spinor->Q, 2);
  nmod_mat_init(x, n, 1, 2);
  nmod_mat_init(x_inv, n, 1, 2);
  nmod_mat_init(diff, n, 1, 2);
  if (ok)
    ok = spinor_neighbor_vector_2(x, mat, q_2);
  if (ok)
    ok = spinor_neighbor_vector_2(x_inv, inv, q_2);

  sign = 0;
  if (ok) {
    // whether x + x' is in the radical
    nmod_mat_add(x, x, x_inv);
    nmod_mat_mul(diff, q_2, x);
    sign = nmod_mat_is_zero(diff) ? -1 : 1;
  }

  nmod_mat_clear(diff);
  nmod_mat_clear(x_inv);
  nmod_mat_clear(x);
  nmod_mat_clear(q_2);
  fmpz_clear(content);
  fmpz_clear(inv_denom);
  fmpz_mat_clear(inv);
  fmpz_clear(two);
  fmpz_clear(unit);

  return sign;
}

// Whether the isometry mat/denom, in lowest terms, where denom is even, can be
// one which takes the lattice to a 2-neighbor of it at a vector outside the
// radical: 2 divides denom exactly once, mat is of rank 1 modulo 2, and its
// column modulo 2 is not in the radical.
static bool spinor_is_neighbor_isometry_2(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom)
{
  nmod_mat_t q_2, x;
  bool ok;

  if (fmpz_fdiv_ui(denom, 4) != 2)
    return false;

  nmod_mat_init_set_fmpz_mat(q_2, spinor->Q, 2);
  nmod_mat_init(x, fmpz_mat_nrows(mat), 1, 2);
  ok = spinor_neighbor_vector_2(x, mat, q_2);
  nmod_mat_clear(x);
  nmod_mat_clear(q_2);

  return ok;
}

// The character at 2 at the isometry g = mat/denom, in lowest terms, which is
// integral at 2 or takes the lattice to a 2-neighbor of it, when the radical
// modulo 2 is of dimension 3, the quadratic form is odd at a vector of it and
// 16 divides the discriminant. It is the sign by which g acts on the deepest
// part of the radical at 2, read from the images under the inverse of g of the
// vectors which span it.
// When this part is a line spanned by e, the image of e is e or -e modulo 4
// times the lattice, and the character is this sign, after the isometry is
// multiplied by -1 if its determinant is -1, whose sign is det_sign. When it is
// a plane, the images of the vectors are in the plane modulo 2 times the
// lattice, and the character is the sign of the permutation of the three
// nonzero vectors of the plane.
// A 2-neighbor is the same as the lattice outside a hyperbolic plane H of it,
// in which H is replaced by a lattice between 2H and H/2. The radicals modulo
// 2^k of the two lattices then differ only by vectors of 2^(k-1) H, where k is
// at least 2, and at least 3 when the deepest part is a line. Hence the two
// lattices have the same deepest part, modulo 2 times either of them, and
// modulo 4 times either of them when it is a line. So the same formula holds
// at an isometry to a 2-neighbor, and it is 1 at an isometry which is the
// identity on the orthogonal complement of a hyperbolic plane.
// Returns 1 or -1, and 0 if the images are not as above.
static int spinor_deep_sign_2(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom, int det_sign)
{
  fmpz_mat_t inv, image;
  fmpz_t inv_denom, odd, two, t;
  slong i, j, k, n, r, val;
  ulong x, e, unit, v[2];
  int idx[2], third, sign;
  bool ok, plus, minus, same;

  n = fmpz_mat_nrows(mat);
  r = fmpz_mat_ncols(spinor->deep_2);

  fmpz_mat_init(inv, n, n);
  fmpz_mat_init(image, n, r);
  fmpz_init(inv_denom);
  fmpz_init(odd);
  fmpz_init(t);
  fmpz_init_set_ui(two, 2);

  // The inverse of the isometry is denom * inv / inv_denom. The images are
  // vectors of the lattice at 2 when the power of 2 in the denominator
  // divides them, and the odd part of the denominator is a unit.
  ok = fmpz_mat_inv(inv, inv_denom, mat);
  if (ok) {
    fmpz_mat_mul(image, inv, spinor->deep_2);
    fmpz_mat_scalar_mul_fmpz(image, image, denom);
    val = fmpz_remove(odd, inv_denom, two);
    // an odd number is its own inverse modulo 4
    unit = fmpz_fdiv_ui(odd, 4);
    for (j = 0; ok && (j < r); j++)
      for (i = 0; ok && (i < n); i++) {
	fmpz_fdiv_q_2exp(t, fmpz_mat_entry(image, i, j), val);
	fmpz_mul_2exp(odd, t, val);
	ok = fmpz_equal(odd, fmpz_mat_entry(image, i, j));
	fmpz_mul_ui(fmpz_mat_entry(image, i, j), t, unit);
      }
  }

  sign = 0;
  if (ok && (r == 1)) {
    plus = minus = true;
    for (i = 0; i < n; i++) {
      x = fmpz_fdiv_ui(fmpz_mat_entry(image, i, 0), 4);
      e = fmpz_fdiv_ui(fmpz_mat_entry(spinor->deep_2, i, 0), 4);
      if (x != e)
	plus = false;
      if ((x + e) % 4 != 0)
	minus = false;
    }
    // e is not 0 modulo 2, so at most one of the two holds
    if (plus)
      sign = det_sign;
    else if (minus)
      sign = -det_sign;
  }
  else if (ok && (r == 2)) {
    // The three nonzero vectors of the plane are those of the two columns
    // and their sum, at the indices 0, 1 and 2.
    for (j = 0; j < 2; j++) {
      idx[j] = -1;
      for (k = 0; k < 3; k++) {
	same = true;
	for (i = 0; i < n; i++) {
	  v[0] = fmpz_fdiv_ui(fmpz_mat_entry(spinor->deep_2, i, 0), 2);
	  v[1] = fmpz_fdiv_ui(fmpz_mat_entry(spinor->deep_2, i, 1), 2);
	  x = (k < 2) ? v[k] : (v[0] + v[1]) % 2;
	  if (fmpz_fdiv_ui(fmpz_mat_entry(image, i, j), 2) != x)
	    same = false;
	}
	if (same)
	  idx[j] = k;
      }
    }
    if ((idx[0] >= 0) && (idx[1] >= 0) && (idx[0] != idx[1])) {
      // The image of the third vector is the sum of the two images. A
      // permutation of three vectors is odd when it fixes exactly one.
      third = 3 - idx[0] - idx[1];
      sign = ((idx[0] == 0) + (idx[1] == 1) + (third == 2) == 1) ? -1 : 1;
    }
  }

  fmpz_clear(two);
  fmpz_clear(t);
  fmpz_clear(odd);
  fmpz_clear(inv_denom);
  fmpz_mat_clear(image);
  fmpz_mat_clear(inv);

  return sign;
}

// TODO - use the modular arithmetic built in inside nmod_t instead of gcdext
// Here mat/denom is the isometry in lowest terms,
// and we're computing the spinor norm at all spinor primes
static W64 spinor_norm_lowest_terms(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom)
{
  nmod_mat_t mat_p, mat_p_t, rad_mat, coords;
  slong idx, row, col, pivot, prime_idx, n, r;
  fmpz_t det; // for some reason det is not implemented for fq_nmod_mat. we reduce the det mod p instead
  ulong denom_p, det_p;
  ulong* evs;
  W64 norm, zas;
  bool has_zas, in_radical;
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
    if ((p.n == 4) && (fmpz_mat_ncols(spinor->deep_2) != 0)) {
      // The sign on the deepest part of the radical at 2, by the same formula
      // at an isometry which is integral at 2 and at one to a 2-neighbor.
      if (fmpz_is_even(denom) && !spinor_is_neighbor_isometry_2(spinor, mat, denom)) {
	fprintf(stderr, "Error: the character at the prime 2 is not defined at an isometry whose denominator is even: the isometry does not take the lattice to a neighbor of it at 2\n");
	exit(EXIT_FAILURE);
      }
      sign = spinor_deep_sign_2(spinor, mat, denom, fmpz_sgn(det) * fmpz_sgn(denom));
      if (sign == 0) {
	fprintf(stderr, "Error: the character at the prime 2 is not defined at an isometry: it does not act by a sign on the deepest part of the radical at 2\n");
	exit(EXIT_FAILURE);
      }
      evs[prime_idx] = (sign == 1) ? 1 : p.n-1;
      continue;
    }
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
      else if ((p.n == 4) && (nmod_mat_nrows(spinor->rads[prime_idx]) != 1)) {
	// The radical modulo 2 is of dimension 3, and the isometry has to take
	// the lattice to a 2-neighbor of it.
	if (!spinor->odd_on_rad_2)
	  sign = spinor_neighbor_sign_2(spinor, mat, denom);
	else if (!spinor_is_neighbor_isometry_2(spinor, mat, denom))
	  sign = 0;
	else {
	  // The character is that of the spinor norm, and the sign is such
	  // that it is 1 at an isometry which is the identity on the orthogonal
	  // complement of a hyperbolic plane, whose spinor norm is 2.
	  if (!has_zas) {
	    zas = spinor_norm_zas_fmpz_mat(spinor, mat, denom);
	    has_zas = true;
	  }
	  sign = ((zas >> prime_idx) & 1) ? 1 : -1;
	}
	if (sign == 0) {
	  fprintf(stderr, "Error: the character at the prime 2 is not defined at an isometry whose denominator is even: the isometry does not take the lattice to a neighbor of it at 2\n");
	  exit(EXIT_FAILURE);
	}
	if (sign == -1)
	  evs[prime_idx] = p.n-1;
      }
      else {
	// For the prime 2 the radical modulo 2 is a line. For an odd prime the
	// two signs are the same at an isometry to a q-neighbor when the radical
	// is a line, and only the first is defined at an isometry to a neighbor
	// of T_{q^2}.
	if (nmod_mat_nrows(spinor->rads[prime_idx]) == 1)
	  sign = spinor_radical_sign(spinor, prime_idx, mat, denom,
				     fmpz_sgn(det) * fmpz_sgn(denom));
	else
	  sign = spinor_neighbor_sign(spinor, prime_idx, mat, denom,
				      fmpz_sgn(det) * fmpz_sgn(denom));
	if (sign == 0) {
	  fprintf(stderr, "Error: the character at the prime %lu is not defined at an isometry whose denominator it divides: the isometry does not take the lattice to a neighbor of it at %lu or at its square\n",
		  spinor_prime(spinor, prime_idx), spinor_prime(spinor, prime_idx));
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

    r = nmod_mat_nrows(spinor->rads[prime_idx]);
    if ((p.n == 4) && (r == 3)) {
      if (!spinor->odd_on_rad_2)
	sign = spinor_quotient_sign_2(spinor, prime_idx, mat_p);
      else if (spinor->exps[prime_idx] == 2)
	sign = spinor_radical_rank_sign_2(spinor, prime_idx, rad_mat);
      else {
	// the parity of the valuation at 2 of the spinor norm
	if (!has_zas) {
	  zas = spinor_norm_zas_fmpz_mat(spinor, mat, denom);
	  has_zas = true;
	}
	sign = ((zas >> prime_idx) & 1) ? -1 : 1;
      }
      evs[prime_idx] = (sign == 1) ? 1 : p.n-1;
    }
    else if (p.n == 4) {
      // the radical modulo 2 is a line
      assert(r == 1);
      pivot = spinor->pivots[prime_idx];
      evs[prime_idx] = nmod_div(nmod_mat_entry(rad_mat, 0, pivot),
				nmod_mat_entry(spinor->rads[prime_idx], 0, pivot), p);
      assert((evs[prime_idx] == 1) || (evs[prime_idx] == p.n-1));
    }
    else {
      // the determinant of the isometry on the radical
      nmod_mat_init(coords, r, r, p.n);
      in_radical = spinor_radical_coordinates(coords, spinor, prime_idx, rad_mat);
      evs[prime_idx] = nmod_mat_det(coords);
      nmod_mat_clear(coords);
      if (!in_radical || ((evs[prime_idx] != 1) && (evs[prime_idx] != p.n-1))) {
	fprintf(stderr, "Error: an isometry which is integral at the prime %lu does not act on the radical of the lattice modulo %lu with determinant 1 or -1\n",
		p.n, p.n);
	exit(EXIT_FAILURE);
      }
    }

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
