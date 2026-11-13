/*****************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : mass.c
 *
 * Description: Functions for computing mass formula
 *
 *****************************************************************
 */

// System dependencies

#include <assert.h>

// Internal dependencies

#include <flint/fmpz.h>
#include <flint/fmpz_factor.h>

#include "arith.h"
#include "mass.h"

#if QF_RANK % 2 == 0
#error "The mass formula is only implemented for lattices of odd rank"
#endif

/* The mass is computed as in
   [CS] J. H. Conway and N. J. A. Sloane, Low-dimensional lattices. IV.
        The mass formula, Proc. R. Soc. Lond. A 419 (1988), 259-286.
   It is the product of the standard mass, which depends only on the rank,
   and of the p-mass divided by the standard p-mass for the primes p
   dividing twice the determinant. */

/* The diagonal factor M_p of a Jordan constituent, as a function of
   its species ([CS], Table 2). */
static void diagonal_factor(fmpq_t M, slong species, const fmpz_t p)
{
  fmpq_t p_inv, p_k, tmp;
  slong n, s, k;

  fmpq_one(M);
  if (species == 0)
    return;

  fmpq_init(p_inv);
  fmpq_init(p_k);
  fmpq_init(tmp);

  n = (species > 0) ? species : -species;
  s = (n + 1) / 2;

  fmpq_set_fmpz(p_inv, p);
  fmpq_inv(p_inv, p_inv);

  // M = 2 * prod_{k = 1}^{s-1} (1 - p^{-2k}), which is then inverted
  fmpq_set_si(M, 2, 1);
  for (k = 1; k < s; k++) {
    fmpq_pow_si(p_k, p_inv, 2*k);
    fmpq_set_si(tmp, 1, 1);
    fmpq_sub(tmp, tmp, p_k);
    fmpq_mul(M, M, tmp);
  }
  if (n % 2 == 0) {
    fmpq_pow_si(p_k, p_inv, s);
    fmpq_set_si(tmp, 1, 1);
    if (species > 0)
      fmpq_sub(tmp, tmp, p_k);
    else
      fmpq_add(tmp, tmp, p_k);
    fmpq_mul(M, M, tmp);
  }
  fmpq_inv(M, M);

  fmpq_clear(tmp);
  fmpq_clear(p_k);
  fmpq_clear(p_inv);

  return;
}

/* The residue modulo 8 of a rational number which is a unit at 2 */
static int unit_mod_8(const fmpq_t u)
{
  // The denominator is odd, so that it is its own inverse modulo 8.
  return (fmpz_fdiv_ui(fmpq_numref(u), 8) * fmpz_fdiv_ui(fmpq_denref(u), 8)) % 8;
}

/* The p-mass of the lattice whose gram matrix is q, divided by the standard
   p-mass ([CS], equations (3) and (6)).
   A Jordan constituent is given by its scale p^e, its rank, and
   for odd p - the class of its determinant modulo squares,
   for p = 2 - its type (I if it is odd, II if it is even), and its octane value,
              which determines its determinant modulo squares, up to sign. */
static void local_mass(fmpq_t m, const square_matrix_t q, const fmpz_t p)
{
  jordan_data_t jordan;
  fmpq_t det, tmp, M;
  fmpz_t p_e;
  slong num_scales, i, j, k, e_min, e_max, n, t, o, ct, exp_two, species;
  slong* rank;
  bool* is_odd;
  int* octane;
  bool is_two, is_free;

  fmpq_init(det);
  fmpq_init(tmp);
  fmpq_init(M);
  fmpz_init(p_e);
  jordan_data_init(jordan, QF_RANK);

  jordan_decomposition(jordan, q, p);

  is_two = fmpz_equal_si(p, 2);

  // The constituents are indexed by the exponent of their scale. There are two
  // more on each side, of rank 0, so that all the constituents which are not
  // listed are of rank 0, and are surrounded by constituents of rank 0.
  e_min = jordan->exponents[0];
  e_max = jordan->exponents[jordan->num_blocks-1];
  num_scales = e_max - e_min + 5;

  rank = (slong*)calloc(num_scales, sizeof(slong));
  is_odd = (bool*)calloc(num_scales, sizeof(bool));
  octane = (int*)calloc(num_scales, sizeof(int));

  fmpq_one(m);
  ct = 0;
  for (i = 0; i < jordan->num_blocks; i++) {
    k = jordan->exponents[i] - e_min + 2;
    n = fmpq_mat_nrows(jordan->grams[i]);
    // jordan_decomposition returns a single block for every scale,
    // in increasing order of the exponents.
    assert(rank[k] == 0);
    rank[k] = n;

    for (j = 0; j < i; j++)
      ct += (jordan->exponents[i] - jordan->exponents[j]) * n * fmpq_mat_nrows(jordan->grams[j]);

    // the determinant of the constituent, divided by its scale
    fmpz_pow_ui(p_e, p, jordan->exponents[i] * n);
    fmpq_mat_det(det, jordan->grams[i]);
    fmpq_div_fmpz(det, det, p_e);

    if (is_two) {
      // The oddity of the constituent is the sum of the diagonal entries
      // that are units. It is the trace of a diagonal form if there is one,
      // since the oddity of an even form of rank 2 is 0.
      // This relies on the form of the blocks which jordan_decomposition
      // returns at 2: an orthogonal sum of forms of rank 1 and of even forms
      // of rank 2, since it splits a vector off whenever the smallest
      // valuation is attained on the diagonal.
      for (j = 0; j < n; j++) {
	fmpq_div_2exp(tmp, fmpq_mat_entry(jordan->grams[i], j, j), jordan->exponents[i]);
	if (fmpz_is_odd(fmpq_numref(tmp))) {
	  is_odd[k] = true;
	  octane[k] += unit_mod_8(tmp);
	}
      }
      o = unit_mod_8(det);
      if ((o == 3) || (o == 5))
	octane[k] += 4;
      octane[k] %= 8;
    }
    else {
      // The species is n, and when n is even it has the sign for which
      // (-1)^(n/2) * det is a square.
      species = n;
      if (n % 2 == 0) {
	if (n % 4 == 2)
	  fmpq_neg(det, det);
	if (!fmpq_is_local_square(det, p))
	  species = -n;
      }
      diagonal_factor(M, species, p);
      fmpq_mul(m, m, M);
    }
  }

  // the cross terms
  assert(ct % 2 == 0);
  fmpz_pow_ui(p_e, p, ct / 2);
  fmpq_mul_fmpz(m, m, p_e);

  if (is_two) {
    // The species of the constituents ([CS], Table 1), including the ones of
    // rank 0 next to the first and last. A constituent is bound if one of its
    // neighbours is of type I, and free otherwise.
    exp_two = 0;
    for (k = 1; k + 1 < num_scales; k++) {
      is_free = !(is_odd[k-1] || is_odd[k+1]);
      n = rank[k];
      o = octane[k];
      t = (is_odd[k] && (n % 2 == 0)) ? (n/2 - 1) : (n/2);
      if (is_free && ((o == 0) || (o == 1) || (o == 7)))
	species = 2*t;
      else if (is_free && ((o == 3) || (o == 4) || (o == 5)))
	species = -2*t;
      else
	species = 2*t + 1;
      diagonal_factor(M, species, p);
      fmpq_mul(m, m, M);

      // The type factor is 2^(n(I,I) - n(II)), where n(I,I) is the number of
      // pairs of adjacent constituents of type I, and n(II) is the sum of the
      // ranks of the constituents of type II.
      if (is_odd[k] && is_odd[k+1])
	exp_two++;
      if (!is_odd[k])
	exp_two -= n;
    }
    if (exp_two >= 0)
      fmpq_mul_2exp(m, m, exp_two);
    else
      fmpq_div_2exp(m, m, -exp_two);
  }

  // The standard p-mass. In odd rank it is the diagonal factor of
  // a single constituent of odd rank.
  diagonal_factor(M, QF_RANK, p);
  fmpq_div(m, m, M);

  free(octane);
  free(is_odd);
  free(rank);
  jordan_data_clear(jordan);
  fmpz_clear(p_e);
  fmpq_clear(M);
  fmpq_clear(tmp);
  fmpq_clear(det);

  return;
}

/* compute the mass of the genus represented by q */
void get_mass(fmpq_t mass, const square_matrix_t q)
{
  fmpz_t det;
  fmpz_mat_t Q;
  fmpz_factor_t fac;
  fmpq_t zeta, m_p;
  slong i, j, k;

  fmpz_init(det);
  fmpq_init(zeta);
  fmpq_init(m_p);
  fmpz_factor_init(fac);
  fmpz_mat_init(Q, QF_RANK, QF_RANK);

  // The standard mass. In odd rank n it is
  // 2 * prod_{k=1}^{(n-1)/2} |zeta(1-2k)|, where zeta(1-2k) = -B_{2k} / (2k).
  fmpq_set_si(mass, 2, 1);
  for (k = 1; k <= QF_RANK / 2; k++) {
    bernoulli_number(zeta, 2*k);
    fmpq_abs(zeta, zeta);
    fmpq_mul(mass, mass, zeta);
    fmpq_set_si(zeta, 1, 2*k);
    fmpq_mul(mass, mass, zeta);
  }

  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      fmpz_set_si(fmpz_mat_entry(Q, i, j), q[i][j]);

  // the primes dividing twice the determinant
  fmpz_mat_det(det, Q);
  fmpz_mul_2exp(det, det, 1);
  fmpz_factor(fac, det);

  for (i = 0; i < fac->num; i++) {
    local_mass(m_p, q, &(fac->p[i]));
    fmpq_mul(mass, mass, m_p);
  }

  fmpz_mat_clear(Q);
  fmpz_factor_clear(fac);
  fmpq_clear(m_p);
  fmpq_clear(zeta);
  fmpz_clear(det);

  return;
}
