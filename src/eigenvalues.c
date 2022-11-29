#include <assert.h>

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpq_poly.h>

#include "eigenvalues.h"
#include "fmpq_poly.h"
#include "hecke.h"
#include "nf_mat.h"

void eigenvalues_init(eigenvalues_t evs, slong num, slong dim)
{
  slong i;
  
  evs->num = num;
  evs->dim = dim;
  evs->nfs = (nf_t*)malloc(num * sizeof(nf_t));
  evs->eigenvals = (nf_elem_t*)malloc(num * sizeof(nf_elem_t));
  evs->eigenvecs = (nf_elem_t**)malloc(num * sizeof(nf_elem_t*));
  for (i = 0; i < num; i++)
    evs->eigenvecs[i] = (nf_elem_t*)malloc(dim * sizeof(nf_elem_t));
  evs->lift_type = (aut_type*)malloc(num * sizeof(aut_type));
			       
  return;
}

void eigenvalues_init_set_mat(eigenvalues_t evs, matrix_TYP* mat)
{
  fmpz_mat_t M;
  fmpz_poly_t cp;
  fmpz_poly_factor_t cp_fac;
  fmpq_poly_t nf_poly;
  int row, col, i, j;
  nf_mat_t M_K, ker, lambda;
  
  fmpz_mat_init(M, mat->rows, mat->cols);

  for (row = 0; row < mat->rows; row++)
    for (col = 0; col  < mat->cols; col++)
      *fmpz_mat_entry(M, row, col) = mat->array.SZ[row][col];

  fmpz_poly_init(cp);
  fmpz_mat_charpoly(cp, M);
#ifdef DEBUG_LEVEL_FULL
  fmpz_poly_print_pretty(cp, "x");
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  fmpz_poly_factor_init(cp_fac);
  fmpz_poly_factor_zassenhaus(cp_fac, cp);

#ifdef DEBUG_LEVEL_FULL
  for (i = 0; i < cp_fac->num; i++) {
    fmpz_poly_print_pretty(&(cp_fac->p[i]), "x");
    printf("\n");
  }
#endif // DEBUG_LEVEL_FULL

  eigenvalues_init(evs, cp_fac->num, mat->rows);
  
  for (i = 0; i < cp_fac->num; i++) {
    fmpq_poly_init(nf_poly);
    fmpq_poly_set_fmpz_poly(nf_poly, &(cp_fac->p[i]));
    nf_init(evs->nfs[i], nf_poly);
    nf_elem_init(evs->eigenvals[i], evs->nfs[i]);
    nf_elem_gen(evs->eigenvals[i], evs->nfs[i]);
    nf_mat_init(M_K, &(evs->nfs[i]), mat->rows, mat->cols);
    nf_mat_init(lambda, &(evs->nfs[i]), mat->rows, mat->cols);
    fmpz_mat_get_nf_mat(M_K, M);
    nf_mat_one(lambda);
    nf_mat_scalar_mul_nf(lambda, lambda, evs->eigenvals[i]);
    nf_mat_sub(M_K, M_K, lambda);
    nf_mat_kernel(ker, M_K);
#ifdef DEBUG
    printf("ker = \n");
    nf_mat_print_pretty(ker);
#endif // DEBUG
    for (j = 0; j < mat->rows; j++) {
      nf_elem_init(evs->eigenvecs[i][j], evs->nfs[i]);
      nf_elem_set(evs->eigenvecs[i][j], *nf_mat_entry(ker, 0, j), evs->nfs[i]);
    }
    nf_mat_clear(ker);
    nf_mat_clear(M_K);
    nf_mat_clear(lambda);
  }

  fmpz_mat_clear(M);
  fmpz_poly_factor_clear(cp_fac);
  fmpz_poly_clear(cp);
  
  return;
}

void eigenvalues_clear(eigenvalues_t evs)
{
  int i, j;

  for (i = 0; i < evs->num; i++) {
    for (j = 0; j < evs->dim; j++) {
      nf_elem_clear(evs->eigenvecs[i][j], evs->nfs[i]);
    }
    free(evs->eigenvecs[i]);
    nf_elem_clear(evs->eigenvals[i], evs->nfs[i]);
    nf_clear(evs->nfs[i]);
  }

  free(evs->eigenvecs);
  free(evs->eigenvals);
  free(evs->nfs);
  free(evs->lift_type);
  
  return;
}

bool nf_elem_is_square(const nf_elem_t x, const nf_t K)
{
  // we do it in a weird way because of the lack of functionality of nf_elem

  slong i;
  fmpq_mat_t x_reg;
  fmpq_poly_t min_poly, x2;
  fmpq_poly_factor_t fac;
  bool is_sqr = false;

  fmpq_poly_init(x2);
  fmpq_poly_init(min_poly);
  fmpq_mat_init(x_reg, fmpq_poly_degree(K->pol), fmpq_poly_degree(K->pol));
  fmpq_poly_set_str(x2, "3  0 0 1");
  
  nf_elem_rep_mat(x_reg, x, K);
  fmpq_mat_minpoly(min_poly, x_reg);
  fmpq_poly_compose(min_poly, min_poly, x2);
  fmpq_poly_factor(fac, min_poly);

  for (i = 0; i < fac->num; i++) {
    if (fac->num > 1) { // Theorem 1.4 from Landau's paper
      is_sqr = true;
      break;
    }
  }

  fmpq_poly_factor_clear(fac);
  fmpq_poly_clear(x2);
  fmpq_poly_clear(min_poly);
  fmpq_mat_clear(x_reg);

  return is_sqr;
}

bool nf_elem_is_square_fast(const nf_elem_t x, const nf_t K)
{
  // we do it in a weird way because of the lack of functionality of nf_elem

  slong i, n;
  slong coeff_int, p_int;
  fmpz_t den, coeff, p, prod, bound;
  nf_elem_t x_int;
  fmpz_poly_t f_int;
  nmod_poly_t f_int_mod_p, g_int_mod_p, r_int_mod_p;
  nmod_poly_factor_t fac;
  fq_nmod_ctx_t F; // the finite field
  fq_nmod_t r;

  n = fmpq_poly_degree(K->pol);
  
  fmpz_init(den);
  fmpz_init(coeff);
  fmpz_init(p);
  fmpz_init(bound);
  nf_elem_init(x_int, K);
  fmpz_poly_init(f_int);
  nmod_poly_factor_init(fac);
  
  nf_elem_get_den(den, x, K);
  fmpz_mul(den, den, den);
  nf_elem_scalar_mul_fmpz(x_int, x, den, K);

  fmpz_one(bound);
  for (i = 0; i < n; i ++) {
    nf_elem_get_coeff_fmpz(coeff, x_int, i, K);
    if (fmpz_cmpabs(bound, coeff) < 0)
      fmpz_abs(bound, coeff);
  }

  fmpz_mul(bound, bound, bound);
  // bound^2 is not enough!! Consider 21/2*(7 + sqrt(13))
  fmpz_mul(bound, bound, bound);
  fmpz_init(prod);
  fmpz_one(prod);

  // This is only probabilistic -
  // should really go only up to bound, reconstruct using CRT and check
  while (fmpz_cmp(prod, bound) < 0) {
    // construct the finite field mod p
    fmpz_nextprime(p, p, true);
    fmpz_mul(prod, prod, p);
    p_int = fmpz_get_si(p);
    fmpq_poly_get_numerator(f_int, K->pol);
    nmod_poly_init2(f_int_mod_p, p_int, n+1);
    for (i = 0; i <= n; i++) {
      fmpz_poly_get_coeff_fmpz(coeff, f_int, i);
      fmpz_mod_ui(coeff, coeff, p_int);
      coeff_int = fmpz_get_ui(coeff);
      /*
      coeff_int = fmpz_poly_get_coeff_si(f_int, i);
      if (coeff_int < 0)
	coeff_int += p_int * ((-coeff_int) / p_int + 1);
      */
      nmod_poly_set_coeff_ui(f_int_mod_p, i, coeff_int);
    }
    
    nmod_poly_init2(g_int_mod_p, p_int, n);
    for (i = 0; i < n; i++) {
      nf_elem_get_coeff_fmpz(coeff, x_int, i, K);
      fmpz_mod_ui(coeff, coeff, p_int);
      coeff_int = fmpz_get_ui(coeff);
      /*
      coeff_int = fmpz_get_si(coeff);
      if (coeff_int < 0)
	coeff_int += p_int * ((-coeff_int) / p_int + 1);
      */
      nmod_poly_set_coeff_ui(g_int_mod_p, i, coeff_int);
    }

    nmod_poly_factor(fac, f_int_mod_p);
    for (i = 0; i < fac->num; i++) {
      nmod_poly_init(r_int_mod_p, p_int);
      nmod_poly_rem(r_int_mod_p, g_int_mod_p, &(fac->p[i]));
      fq_nmod_ctx_init_modulus(F, &(fac->p[i]), "x");
      fq_nmod_init(r, F);
      fq_nmod_set(r, r_int_mod_p, F);
      if (!fq_nmod_is_square(r, F)) {
	fq_nmod_clear(r, F);
	fq_nmod_ctx_clear(F);
	nmod_poly_clear(r_int_mod_p);
	nmod_poly_clear(g_int_mod_p);
	nmod_poly_clear(f_int_mod_p);
	nf_elem_clear(x_int, K);
	nmod_poly_factor_clear(fac);
	fmpz_poly_clear(f_int);
	fmpz_clear(p);
	fmpz_clear(coeff);
	fmpz_clear(den);
	return false;
      }
      fq_nmod_clear(r, F);
      fq_nmod_ctx_clear(F);
      nmod_poly_clear(r_int_mod_p);
    }

    nmod_poly_clear(g_int_mod_p);
    nmod_poly_clear(f_int_mod_p);
  }
  fmpz_clear(prod);

  nf_elem_clear(x_int, K);
  nmod_poly_factor_clear(fac);
  fmpz_poly_clear(f_int);
  fmpz_clear(p);
  fmpz_clear(coeff);
  fmpz_clear(den);
  
  return true;
}

aut_type ev_is_lpoly_reducible(const eigenvalues_t evs, slong ev_idx, slong p, slong c, const genus_t genus)
{
  int k;
  nf_elem_t a[2];
  nf_elem_t disc; // discriminant of the quadratic controlling lifts
  nf_elem_t root;
  bool is_lift;
  aut_type lift_type;
  
  for (k = 1; k < 3; k++) {
    get_hecke_ev_nbr_data_all_conductors(a[k-1], genus, evs, p, k, ev_idx, c);
  }

#ifdef DEBUG
  for (k = 1; k < 3; k++) {
    printf("a[%d] = ", k-1); 
    nf_elem_print_pretty(a[k-1], evs->nfs[ev_idx], "a");
    printf("\n");
  }
#endif // DEBUG
  
  nf_elem_init(disc, evs->nfs[ev_idx]);

  // disc = a[0]^2 - 4*p*(a[1] + 1 - p^2)
  nf_elem_mul(disc, a[0], a[0], evs->nfs[ev_idx]);
  nf_elem_add_si(a[1], a[1], 1 - p*p, evs->nfs[ev_idx]);
  nf_elem_scalar_mul_si(a[1], a[1], 4*p, evs->nfs[ev_idx]);
  nf_elem_sub(disc, disc, a[1], evs->nfs[ev_idx]); 

#ifdef DEBUG
  printf("disc = "); 
  nf_elem_print_pretty(disc, evs->nfs[ev_idx], "a");
  printf(" over ");
  nf_print(evs->nfs[ev_idx]);
  printf("\n");
#endif // DEBUG

  // !! TODO - this is not good enough, we might miss forms!
  // e.g. 889. otoh, nf_elem_is_square is too slow.
  // Should do CRT reconstruction
  is_lift = nf_elem_is_square_fast(disc, evs->nfs[ev_idx]);
  assert(is_lift == nf_elem_is_square(disc, evs->nfs[ev_idx]));

  lift_type = is_lift ? Y : G;

  // if a lift, find which type of lift
  if (is_lift) {
    nf_elem_init(root, evs->nfs[ev_idx]);
    // check whether (2*p*(p+1) - a[0])^2 == disc for Saito-Kurokawa, i.e. p(p+1) is a root
    nf_elem_sub_si(root, a[0], 2*p*(p+1), evs->nfs[ev_idx]);
    nf_elem_mul(root, root, root, evs->nfs[ev_idx]);
    if (nf_elem_equal(root, disc, evs->nfs[ev_idx]))
      lift_type = P;
    
    // check whether Siegel-Eisenstein - a[0] = p^3 + p^2 + p + 1
    if ( (lift_type == P) &&
	 (nf_elem_equal_si(a[0], (p*p*p + p*p + p + 1), evs->nfs[ev_idx])) )
      lift_type = F;
  }
  
  nf_elem_clear(disc, evs->nfs[ev_idx]);
  for (k = 1; k < 3; k++) {
    nf_elem_clear(a[k-1], evs->nfs[ev_idx]);
  }
  
  return lift_type;
}

void eigenvalues_set_lifts(eigenvalues_t evs, slong prec, slong c, const genus_t genus)
{
  slong cnt, i;
  fmpz_t p;
  aut_type lift_type;
  
  fmpz_init_set_si(p,2); // p = 2
  
  for (cnt = 0; cnt < prec; cnt++) {
    while (fmpz_divisible(genus->disc, p))
      fmpz_nextprime(p, p, true);
    for (i = 0; i < evs->num; i++) {
      lift_type = ev_is_lpoly_reducible(evs, i, fmpz_get_si(p), c, genus);
      if (cnt == 0) 
	evs->lift_type[i] = lift_type;
      else
	if (evs->lift_type[i] != lift_type)
	  evs->lift_type[i] = G;
    }
  }

  fmpz_clear(p);
  
  return;
}

slong ev_get_pivot(const nf_elem_t* evec, const nf_t K, const genus_t genus, slong cond)
{
  slong pivot, i, max_aut;

  max_aut = 1;
  pivot = -1;
  for (i = 0; i < genus->dims[cond]; i++) {
    if (!nf_elem_is_zero(evec[i], K)) {
      if (max_aut < genus->num_auts[cond][i]) {
	max_aut = genus->num_auts[cond][i];
	pivot = i;
      }
    }
  }
  
  return pivot;
}
