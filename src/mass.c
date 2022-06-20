#include <assert.h>

#include "arith.h"
#include "mass.h"
#include "qf_inv.h"

void local_factor(fmpq_t f, const fmpq_mat_t g, const fmpz_t p)
{
  fmpz_t p_sqr;
  fmpq_t p_i, one, tmp, d;
  ulong i, m, r;
  slong sign;

  fmpz_init(p_sqr);
  fmpq_init(p_i);
  fmpq_init(one);
  fmpq_init(tmp);
  m = fmpq_mat_ncols(g);
  
  fmpz_mul(p_sqr, p, p);
  fmpq_one(f);
  fmpq_one(one);
  fmpq_div_fmpz(p_i, one, p_sqr);

  for (i = 2; i+2 <= m; i+= 2) {
    fmpq_sub(tmp, one, p_i);
    fmpq_mul(f, f, tmp);
    fmpq_div_fmpz(p_i, p_i, p_sqr);
  }
  if (m % 2 == 1) {
    if (m != 1) {
      fmpq_sub(tmp, one, p_i);
      fmpq_mul(f, f, tmp);
    }
    fmpq_clear(tmp);
    fmpq_clear(one);
    fmpq_clear(p_i);
    fmpz_clear(p_sqr);
    return;
  }
  r = m / 2;
  sign = (r % 2 == 0) ? 1 : -1;
  fmpq_init(d);
  fmpq_mat_det(d,g);
  if (sign == -1)
    fmpq_neg(d,d);
  
  if (((fmpq_valuation(d,p)) % 2) == 0) {
    fmpq_one(p_i);
    for (i = 0; i < r; i++)
      fmpq_div_fmpz(p_i, p_i, p);

    if (fmpq_is_local_square(d,p)) {
      fmpq_sub(tmp, one, p_i);
      fmpq_mul(f, f, tmp);
    }
    else {
      fmpq_add(tmp, one, p_i);
      fmpq_mul(f, f, tmp);
    }
  }
  fmpq_clear(d);
  fmpq_clear(tmp);
  fmpq_clear(one);
  fmpq_clear(p_i);
  fmpz_clear(p_sqr);
  return;
}

void diagonal_join(fmpq_mat_t joined, const jordan_data_t jordan)
{
  ulong i, j, k;
  ulong joined_j, joined_k;

  joined_j = joined_k = 0;
  for (i = 0; i < jordan->num_blocks; i++) {
    for (j = 0; j < fmpq_mat_nrows(jordan->grams[i]); j++)
      for (k = 0; k < fmpq_mat_ncols(jordan->grams[i]); k++)
	fmpq_set(fmpq_mat_entry(joined, joined_j + j, joined_k + k), fmpq_mat_entry(jordan->grams[i], j, k));
    joined_j += fmpq_mat_nrows(jordan->grams[i]);
    joined_k += fmpq_mat_ncols(jordan->grams[i]);
  }

  return;
}

void _combine(fmpq_t mass, const matrix_TYP* q, fmpz_t p)
{
  assert(!fmpz_equal_si(p,2));

  jordan_data_t jordan;
  fmpq_t e, f, tmp1, tmp2, lf, n_rat, denom, p_e;
  fmpz_mat_t q_mat;
  fmpq_mat_t diag;
  fmpz_t d, e_int;
  ulong* ms;
  ulong i, j, m, t, v, n;

  n = q->rows;
  fmpq_init(e);
  fmpq_init(f);
  fmpq_init(tmp1);
  fmpq_init(tmp2);
  fmpq_init(lf);
  fmpz_mat_init(q_mat, n, n);
  fmpq_mat_init(diag, n, n);
  fmpz_init(d);
  fmpz_init(e_int);
  fmpq_init(p_e);
  fmpq_init(denom);
  jordan_data_init(jordan, n);
  
  jordan_decomposition(jordan, q, p);

  fmpq_one(f);
  fmpq_zero(e);

  ms = (ulong*)malloc(jordan->num_blocks*sizeof(ulong));
  m = 0;
  for (i = 0; i < jordan->num_blocks; i++) {
    ms[i] = fmpq_mat_ncols(jordan->grams[i]);
    m += ms[i];
  }

  for (i = 0; i < jordan->num_blocks; i++) {
    t = (i == 0) ? 0 : jordan->exponents[i-1];
    fmpq_set_si(tmp1, (jordan->exponents[i]-t)*(m+1)*m, 2);
    fmpq_set_si(tmp2, jordan->exponents[i]*(n+1)*ms[i],2);
    fmpq_add(e, e, tmp1);
    fmpq_sub(e, e, tmp2);
    local_factor(lf,jordan->grams[i], p);
    fmpq_mul(f, f, lf);
    m -= ms[i];
  }
  
  // !! We might run into trouble at 2 here
  // check if we need disc or half-disc
  // size_t v = Math<R>::valuation(q.discriminant(), p);
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      fmpz_set_si(fmpz_mat_entry(q_mat, i, j), q->array.SZ[i][j]);

  fmpz_mat_det(d, q_mat);
  v = fmpz_valuation(d, p);

  if ((n % 2 == 0) && (v % 2 == 1)) {
    fmpq_set_si(n_rat, n-1, 2);
    fmpq_add(e, e, n_rat);
  }

  assert(fmpq_is_integral(e));

  fmpz_divexact(e_int, fmpq_numref(e), fmpq_denref(e));
  fmpq_one(p_e);
  fmpq_mul_fmpz(p_e, p_e, p);
  fmpq_pow_si(p_e, p_e, fmpz_get_si(e_int));

  diagonal_join(diag, jordan);
  local_factor(lf, diag, p);
  fmpq_mul(denom, f, p_e); 
  fmpq_div(lf, lf, denom);
  fmpq_div_2exp(lf, lf, jordan->num_blocks-1);

  fmpq_mul(mass, mass, lf);

  fmpq_clear(denom);
  fmpq_clear(tmp1);
  fmpq_clear(tmp2);
  fmpq_clear(e);
  fmpq_clear(f);
  fmpq_clear(lf);
  fmpz_mat_clear(q_mat);
  fmpq_mat_clear(diag);
  fmpz_clear(d);
  fmpz_clear(e_int);
  fmpq_clear(p_e);
  jordan_data_clear(jordan);
  
  return;
}

void get_mass(fmpq_t mass, const matrix_TYP* q)
{
  int n = q->rows;
  int r = n / 2;
  qf_inv_t witt;
  fmpz_t det, two, two_i, mass_two, six, disc, r_Z;
  fmpz_mat_t hasse, B, Q;
  int num_hasse, num_B;
  fmpz_factor_t fac;
  int i, j, cmp, val2, w;
  int has_two;
  fmpq_t bernoulli, mass_r, mass_r_1, mass_2r_1;

  fmpz_init(det);
  fmpz_init(disc);
  fmpz_init(two_i);
  fmpz_init(mass_two);
  fmpq_init(bernoulli);
  fmpz_init_set_si(r_Z, -r);
  fmpz_init_set_si(two,2);
  fmpz_init_set_si(six,6);
  fmpz_factor_init(fac);

  fmpz_mat_init(Q, q->rows, q->cols);
  for (i = 0; i < q->rows; i++)
    for (j = 0; j < q->cols; j++)
      fmpz_set_si(fmpz_mat_entry(Q, i, j), q->array.SZ[i][j]);
  
  invariants(det, witt, Q);

  fmpz_mat_init(hasse, 1, witt->num_stored);
  num_hasse = witt_to_hasse(hasse, det, witt);
  fmpz_factor(fac, det);
  fmpz_mat_init(B, 1, num_hasse + fac->num);

  num_B = i = j = 0;
  has_two = FALSE;
  // TODO : Is hasse sorted at this point?
  // merge sort into B
  while ((i < fac->num) && (j < num_hasse)) {
    cmp = fmpz_cmp(&(fac->p[i]), fmpz_mat_entry(hasse, 0, j));
    if (cmp < 0) {
      if (!fmpz_equal_si(&(fac->p[i]), 2))
	fmpz_set(fmpz_mat_entry(B, 0, num_B++), &(fac->p[i]));
      i++;
    }
    if (cmp > 0) {
      if (!fmpz_equal_si(fmpz_mat_entry(hasse, 0, j), 2))
	fmpz_set(fmpz_mat_entry(B, 0, num_B++), fmpz_mat_entry(hasse, 0, j));
      else
	has_two = TRUE;
      j++;
    }
    if (cmp == 0) {
      if (!fmpz_equal_si(&(fac->p[i]), 2))
	fmpz_set(fmpz_mat_entry(B, 0, num_B++), &(fac->p[i]));
      else
	has_two = TRUE;
      i++;
      j++;
    }
  }

  // originally, leave hasse only with 2, if it had a 2.
  
  val2 = fmpz_valuation(det,two);

  // mass from infinity and 2
  fmpq_set_si(mass, 1, 1 << r);
  
  for (i = 1; i < n / 2 + n % 2; i++) {
    fmpz_set_si(two_i, -2*i);
    bernoulli_number(bernoulli, 2*i);
    fmpq_div_fmpz(bernoulli, bernoulli, two_i);
    fmpq_mul(mass, mass, bernoulli);
  }
     
  if (n % 2 == 1) {	 
    if (val2 % 2 == 1) {
      fmpz_set_si(mass_two, (1 << r) + (has_two ? -1 : 1));
      fmpq_mul_fmpz(mass, mass, mass_two); 
      fmpq_div_fmpz(mass, mass, two);
      has_two = FALSE;
    }
    else if (has_two) {
      fmpz_set_si(mass_two, (1 << (n-1))-1);
      fmpq_mul_fmpz(mass, mass, mass_two);
      fmpq_div_fmpz(mass, mass, six);
    }
  }
  else {

    if (r % 2 == 1)
      fmpz_neg(disc, det);
    else
      fmpz_set(disc, det);

    if (fmpz_is_square(disc)) {
      bernoulli_number(bernoulli, r);
      fmpq_div_fmpz(bernoulli, bernoulli, r_Z);;
      fmpq_mul(mass, mass, bernoulli);
    } else {
      bernoulli_number_chi(bernoulli, r, disc);
      fmpq_div_fmpz(bernoulli, bernoulli, r_Z);
      fmpq_mul(mass, mass, bernoulli);
      if (r % 2 == 0) {
	bernoulli_number(bernoulli, r);
	fmpq_div_fmpz(bernoulli, bernoulli, r_Z);
	fmpq_mul(mass, mass, bernoulli);
      }

      if (val2 % 2 == 1) {
	fmpq_div_2exp(mass, mass, 1);
	has_two = FALSE;
      }
    }
    if (has_two) {
      // checking if disc is a local square at 2
      w = 1;
      fmpz_fdiv_q_2exp(disc, disc, val2);
      fmpz_mod_ui(disc, disc, 8); 
      
      if ((val2 % 2 != 1) && (fmpz_is_one(disc)))
	w = -1;

      fmpq_init(mass_r);
      fmpq_init(mass_r_1);
      fmpq_init(mass_2r_1);
      fmpq_mul_2exp(mass_r_1,mass,r-1);
      fmpq_mul_2exp(mass_r, mass_r_1, 1);
      fmpq_mul_2exp(mass_2r_1, mass_r, r-1);
      fmpq_add(mass, mass, mass_2r_1);
      if (w == 1) {
	fmpq_add(mass, mass, mass_r);
	fmpq_add(mass, mass, mass_r_1);
      } else {
	fmpq_sub(mass, mass, mass_r);
	fmpq_sub(mass, mass, mass_r_1);
      }
      fmpq_clear(mass_2r_1);
      fmpq_clear(mass_r_1);
      fmpq_clear(mass_r);
      fmpq_div_fmpz(mass, mass, six);
    }
  }
  
  // odd places which are not unimodular or have Witt invariant -1.
  for (i = 0; i < num_B; i++) {
    _combine(mass, q, fmpz_mat_entry(B, 0, i));
  }

  fmpq_abs(mass, mass);
  
  fmpz_mat_clear(B);
  fmpz_mat_clear(hasse);
  fmpz_factor_clear(fac);
  fmpz_clear(six);
  fmpz_clear(two);
  fmpz_clear(r_Z);
  fmpq_clear(bernoulli);
  fmpz_clear(mass_two);
  fmpz_clear(two_i);
  fmpz_clear(disc);
  fmpz_clear(det);
  
  return;

  /* for (Integer<R> p : B) */
  /*   mass *= Genus<R,n>::_combine(q,p); */
     
  /* return mass.abs(); */
}
