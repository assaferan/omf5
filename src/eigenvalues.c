/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : eigenvalues.c
 *
 * Description: Compute eigenforms and their eigenvalues
 *              in spaces of orthogonal modular forms.
 *
 **********************************************************
 */

// System dependencies

#include <assert.h>

// Required packages dependencies

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>

// Self dependencies

#include <flint/nmod_poly.h>
#include <flint/nmod_poly_factor.h>

#include "eigenvalues.h"
#include "fmpq_poly.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "nf_mat.h"

// initialize the eigenvalues struct - allocate according memory
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
  evs->is_old = (bool*)malloc(num * sizeof(bool));
  for (i = 0; i < num; i++)
    evs->is_old[i] = false;
  // These are initialized only for the forms with is_old.
  evs->old_basis = (fmpq_mat_t*)malloc(num * sizeof(fmpq_mat_t));
  evs->old_hecke = (fmpq_mat_t*)malloc(num * sizeof(fmpq_mat_t));

  return;
}

// compute eigenvecotrs and eigenvalues for the matrix mat
void eigenvalues_init_set_mat(eigenvalues_t evs, matrix_TYP* mat)
{
  fmpz_mat_t M;
  fmpz_poly_t cp;
  fmpz_poly_factor_t cp_fac;
  fmpq_poly_t nf_poly;
  int row, col, i, j;
  nf_mat_t M_K, ker, lambda;

  // start by creating an fmpz_mat matrix from mat
  // !! TODO - maybe we now have a function that does that,
  // and we can replace the following lines by a function?
  fmpz_mat_init(M, mat->rows, mat->cols);

  for (row = 0; row < mat->rows; row++)
    for (col = 0; col  < mat->cols; col++)
      *fmpz_mat_entry(M, row, col) = mat->array.SZ[row][col];

  // We then factor the characteristic polynomial
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

  // and initialize the eigenvalues data according to the number of factors
  // (we only save one for each Galois orbit)
  eigenvalues_init(evs, cp_fac->num, mat->rows);
  
  for (i = 0; i < cp_fac->num; i++) {
    // We create a number field , K, defined by the irreducible factor
    fmpq_poly_init(nf_poly);
    fmpq_poly_set_fmpz_poly(nf_poly, &(cp_fac->p[i]));
    nf_init(evs->nfs[i], nf_poly);
    // and set the eigenvalue to be the generator of the number field,
    // i.e. a root of this polynomial
    nf_elem_init(evs->eigenvals[i], evs->nfs[i]);
    nf_elem_gen(evs->eigenvals[i], evs->nfs[i]);
    // We now proceed to compute the eigenvector.
    // For that we need to base change our matrix to K
    nf_mat_init(M_K, &(evs->nfs[i]), mat->rows, mat->cols);
    nf_mat_init(lambda, &(evs->nfs[i]), mat->rows, mat->cols);
    fmpz_mat_get_nf_mat(M_K, M);
    // lambda is the scalar matrix whose value is our eigenvalue 
    nf_mat_one(lambda);
    nf_mat_scalar_mul_nf(lambda, lambda, evs->eigenvals[i]);
    nf_mat_sub(M_K, M_K, lambda);
    // Kernel(M - lambda I) is the eigenspace
    nf_mat_kernel(ker, M_K);
#ifdef DEBUG
    printf("ker = \n");
    nf_mat_print_pretty(ker);
#endif // DEBUG
    // We add all basis vectors as eigenvectors
    // When we have multiplicity one, these will always be one dimensional,
    // and we will obtain the unique eigenvector with eigenvalue lambda.
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

// Clear the memory allocated for eigenvalues data
void eigenvalues_clear(eigenvalues_t evs)
{
  int i, j;

  // The field, the eigenvector and the eigenvalue are initialized for every
  // form, whatever its type is. This includes the old forms, and the forms
  // which were not selected by their indices.
  for (i = 0; i < evs->num; i++) {
    for (j = 0; j < evs->dim; j++) {
      nf_elem_clear(evs->eigenvecs[i][j], evs->nfs[i]);
    }
    free(evs->eigenvecs[i]);
    nf_elem_clear(evs->eigenvals[i], evs->nfs[i]);
    nf_clear(evs->nfs[i]);
    if (evs->is_old[i]) {
      fmpq_mat_clear(evs->old_basis[i]);
      fmpq_mat_clear(evs->old_hecke[i]);
    }
  }

  free(evs->eigenvecs);
  free(evs->eigenvals);
  free(evs->nfs);
  free(evs->lift_type);
  free(evs->is_old);
  free(evs->old_basis);
  free(evs->old_hecke);
  
  return;
}

// Test whether T acts on the space of an old form as an element of its field
bool ev_is_scalar_on_old_space(const eigenvalues_t evs, slong ev_idx, const fmpq_mat_t T)
{
  slong i, j, pivot, dim_W;
  fmpq_mat_t T_W, WT, both, q_A, tmp;
  nf_elem_t* evec_T;
  nf_elem_t lambda;
  fmpq_poly_t q;
  fmpq_t coeff;
  bool is_scalar;

  assert(evs->is_old[ev_idx]);
  assert((fmpq_mat_nrows(T) == evs->dim) && (fmpq_mat_ncols(T) == evs->dim));

  dim_W = fmpq_mat_nrows(evs->old_basis[ev_idx]);

  // The space W of the form is preserved by T when the rows of W T add
  // nothing to the rank of the basis of W. This is checked first, since
  // restrict_mat is for a space which is preserved.
  fmpq_mat_init(WT, dim_W, evs->dim);
  fmpq_mat_mul(WT, evs->old_basis[ev_idx], T);
  fmpq_mat_init(both, 2 * dim_W, evs->dim);
  for (i = 0; i < dim_W; i++)
    for (j = 0; j < evs->dim; j++) {
      fmpq_set(fmpq_mat_entry(both, i, j), fmpq_mat_entry(evs->old_basis[ev_idx], i, j));
      fmpq_set(fmpq_mat_entry(both, dim_W + i, j), fmpq_mat_entry(WT, i, j));
    }
  fmpq_mat_init(tmp, 2 * dim_W, evs->dim);
  is_scalar = (fmpq_mat_rref(tmp, both) == dim_W);
  fmpq_mat_clear(tmp);
  fmpq_mat_clear(both);
  fmpq_mat_clear(WT);

  pivot = -1;
  for (i = 0; (i < evs->dim) && (pivot == -1); i++)
    if (!nf_elem_is_zero(evs->eigenvecs[ev_idx][i], evs->nfs[ev_idx]))
      pivot = i;

  if ((!is_scalar) || (pivot == -1))
    return false;

  fmpq_mat_init(T_W, dim_W, dim_W);
  restrict_mat(T_W, T, evs->old_basis[ev_idx]);

  // If T is an element of the field on W, the vector of the form is an
  // eigenvector of T, and the element is the eigenvalue. So we take the
  // number lambda which the vector gives at one of its entries,
  evec_T = (nf_elem_t*)malloc(evs->dim * sizeof(nf_elem_t));
  for (i = 0; i < evs->dim; i++)
    nf_elem_init(evec_T[i], evs->nfs[ev_idx]);
  nf_elem_vec_mul_fmpq_mat(evec_T, (const nf_elem_t*)evs->eigenvecs[ev_idx], T, evs->nfs[ev_idx]);
  nf_elem_init(lambda, evs->nfs[ev_idx]);
  nf_elem_div(lambda, evec_T[pivot], evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]);
  for (i = 0; i < evs->dim; i++)
    nf_elem_clear(evec_T[i], evs->nfs[ev_idx]);
  free(evec_T);

  // write it as a polynomial q in the generator of the field,
  fmpq_poly_init(q);
  nf_elem_get_fmpq_poly(q, lambda, evs->nfs[ev_idx]);
  nf_elem_clear(lambda, evs->nfs[ev_idx]);

  // and compare T on W with q(A), where A is the matrix on W of the
  // generator. This is the Horner scheme: fmpq_poly_evaluate_fmpq_mat is of
  // no use here, since it rescales its result.
  fmpq_mat_init(q_A, dim_W, dim_W);
  fmpq_mat_init(tmp, dim_W, dim_W);
  fmpq_init(coeff);
  for (i = fmpq_poly_degree(q); i >= 0; i--) {
    fmpq_mat_mul(tmp, q_A, evs->old_hecke[ev_idx]);
    fmpq_mat_swap(q_A, tmp);
    fmpq_poly_get_coeff_fmpq(coeff, q, i);
    for (j = 0; j < dim_W; j++)
      fmpq_add(fmpq_mat_entry(q_A, j, j), fmpq_mat_entry(q_A, j, j), coeff);
  }
  is_scalar = fmpq_mat_equal(T_W, q_A);

  fmpq_clear(coeff);
  fmpq_mat_clear(tmp);
  fmpq_mat_clear(q_A);
  fmpq_poly_clear(q);
  fmpq_mat_clear(T_W);

  return is_scalar;
}

// Check if an element of a number field is a square
// !! TODO - this should be in a different file
// returns true if and only if x is a square in K
bool nf_elem_is_square(const nf_elem_t x, const nf_t K)
{
  // we do it in a weird way because of the lack of functionality of nf_elem
  // Based on S. Landau, Factoring Polynomials over Algebraic Number Fields, SIAM Journal on Computing, 1985.
  // where we simply try to factor t^2 - x
  //
  // Let B = K[y] / (y^2 - x). When x is not 0, B is a field if x is not
  // a square in K, and the product of two copies of K if it is.
  // An element of B whose characteristic polynomial over Q is squarefree
  // generates B over Q, and then B is a field if and only if this
  // polynomial is irreducible.
  // We take the elements s*a + y, where a is the generator of K, for
  // s = 0, 1, 2, ...: all but finitely many of them are such elements.
  // The characteristic polynomial of s*a + y is the norm of (t - s*a)^2 - x.
  // For s = 0 it is f(t^2), where f is the characteristic polynomial of x,
  // and it is squarefree when x generates K.

  slong i, j, n, s;
  fmpq_mat_t a_reg, x_reg, theta;
  fmpq_poly_t char_poly;
  fmpz_poly_t char_poly_int;
  fmpz_poly_factor_t fac;
  nf_elem_t a;
  bool is_sqr = false;
  bool is_squarefree = false;

  if (nf_elem_is_zero(x, K))
    return true;

  n = fmpq_poly_degree(K->pol);

  nf_elem_init(a, K);
  fmpq_mat_init(a_reg, n, n);
  fmpq_mat_init(x_reg, n, n);
  fmpq_mat_init(theta, 2*n, 2*n);
  fmpq_poly_init(char_poly);
  fmpz_poly_init(char_poly_int);

  // a_reg and x_reg are the matrices representing multiplication by a and by x
  nf_elem_gen(a, K);
  nf_elem_rep_mat(a_reg, a, K);
  nf_elem_rep_mat(x_reg, x, K);

  for (s = 0; !is_squarefree; s++) {
    // theta is the matrix representing multiplication by s*a + y on B,
    // in the basis 1, a, ..., a^(n-1), y, a*y, ..., a^(n-1)*y
    fmpq_mat_zero(theta);
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
	fmpq_mul_si(fmpq_mat_entry(theta, i, j), fmpq_mat_entry(a_reg, i, j), s);
	fmpq_set(fmpq_mat_entry(theta, n+i, n+j), fmpq_mat_entry(theta, i, j));
	fmpq_set(fmpq_mat_entry(theta, n+i, j), fmpq_mat_entry(x_reg, i, j));
      }
      fmpq_one(fmpq_mat_entry(theta, i, n+i));
    }
    fmpq_mat_charpoly(char_poly, theta);
    is_squarefree = fmpq_poly_is_squarefree(char_poly);
  }

  // The content of the polynomial is not among the factors
  fmpq_poly_get_numerator(char_poly_int, char_poly);
  fmpz_poly_factor_init(fac);
  fmpz_poly_factor(fac, char_poly_int);

  is_sqr = (fac->num > 1);

  fmpz_poly_factor_clear(fac);
  fmpz_poly_clear(char_poly_int);
  fmpq_poly_clear(char_poly);
  fmpq_mat_clear(theta);
  fmpq_mat_clear(x_reg);
  fmpq_mat_clear(a_reg);
  nf_elem_clear(a, K);

  return is_sqr;
}

// Check if an element of a number field is a square
// !! TODO - this should be in a different file
// returns false only if x is not a square in K
// The above function might be occassionaly slow.
// This is a faster approach, which reduces x modulo small primes.
// An element which is not a square can be a square modulo all of them.
bool nf_elem_is_square_fast(const nf_elem_t x, const nf_t K)
{
  // we do it in a weird way because of the lack of functionality of nf_elem

  slong i, n;
  slong coeff_int, p_int;
  fmpz_t den, coeff, p, prod, bound, disc;
  bool p_sq_divides_disc;
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

  // bound is the maximal absolute value of a coefficient in the defining polynomial
  fmpz_one(bound);
  for (i = 0; i < n; i ++) {
    nf_elem_get_coeff_fmpz(coeff, x_int, i, K);
    if (fmpz_cmpabs(bound, coeff) < 0)
      fmpz_abs(bound, coeff);
  }

  fmpz_mul(bound, bound, bound);
  // bound^2 is not enough!! Consider 21/2*(7 + sqrt(13))
  fmpz_mul(bound, bound, bound);
  // Are we sure that bound^3 is enough?
  fmpz_init(prod);
  fmpz_one(prod);

  // The discriminant of the defining polynomial. The index of Z[a] in the
  // ring of integers of K can only be divisible by p when p^2 divides it.
  fmpz_init(disc);
  fmpq_poly_get_numerator(f_int, K->pol);
  fmpz_poly_discriminant(disc, f_int);

  // This is only probabilistic -
  // should really go only up to bound, reconstruct using CRT and check
  while (fmpz_cmp(prod, bound) < 0) {
    // construct the finite field mod p
    fmpz_nextprime(p, p, true);
    fmpz_mul(prod, prod, p);
    p_int = fmpz_get_si(p);
    p_sq_divides_disc = fmpz_divisible_si(disc, p_int * p_int);
    // We now want to check whether it is square in each of the extensions
    // of our number field over p
    nmod_poly_init2(f_int_mod_p, p_int, n+1);
    // Construct the defining polynomial with coefficients in F_p
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

    // Construct the polynomial represeting the element x over F_p
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

    // Factor f_int_mod_p to construct the different extensions of F_p
    nmod_poly_factor(fac, f_int_mod_p);
    for (i = 0; i < fac->num; i++) {
      // At a repeated factor, when p^2 divides the discriminant, Z[a] might
      // not be maximal at p. Then the residue field of K can be larger
      // than F, and x can be a square in K without being a square in F,
      // so we cannot test there. (e.g. -4 in Q[a] / (a^2+9), at p = 3)
      if ((fac->exp[i] > 1) && p_sq_divides_disc)
	continue;
      // Find the image of x in this finite field
      nmod_poly_init(r_int_mod_p, p_int);
      nmod_poly_rem(r_int_mod_p, g_int_mod_p, &(fac->p[i]));
      // construct the extension F = F_q 
      fq_nmod_ctx_init_modulus(F, &(fac->p[i]), "x");
      fq_nmod_init(r, F);
      fq_nmod_set(r, r_int_mod_p, F);
      // If x is not a square in F, then it is not a square. 
      if (!fq_nmod_is_square(r, F)) {
	fq_nmod_clear(r, F);
	fq_nmod_ctx_clear(F);
	nmod_poly_clear(r_int_mod_p);
	nmod_poly_clear(g_int_mod_p);
	nmod_poly_clear(f_int_mod_p);
	nf_elem_clear(x_int, K);
	nmod_poly_factor_clear(fac);
	fmpz_poly_clear(f_int);
	fmpz_clear(disc);
	fmpz_clear(prod);
	fmpz_clear(bound);
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
  fmpz_clear(disc);
  fmpz_clear(prod);
  fmpz_clear(bound);

  nf_elem_clear(x_int, K);
  nmod_poly_factor_clear(fac);
  fmpz_poly_clear(f_int);
  fmpz_clear(p);
  fmpz_clear(coeff);
  fmpz_clear(den);

  // If we are here, x is a square in K_P for all P above p for all primes up to bound
  
  return true;
}

// This function checks whether the L-polynomial is reducible (indicating that the form is a lift)
// and returns the corresponding automorphic type.
// The relevant form is indicated by ev_idx (so f is the form in that index in evs)
// p is the prime at which to compute the L-polynomial,
// c is the conductor
// and genus is the genus
aut_type ev_is_lpoly_reducible(const eigenvalues_t evs, slong ev_idx, slong p, slong c, const genus_t genus)
{
  int k;
  nf_elem_t a[2];
  nf_elem_t disc; // discriminant of the quadratic controlling lifts
  nf_elem_t root;
  bool is_lift, is_square_fast;
  aut_type lift_type;

  // Since we are in rank 5, there are only two Hecke operators T_{p,1} and T_{p,2} at every prime.
  // We store the Hecke eigenvalues of f for these operators
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

  // It appears that the reduciblity of the L-polynomial is controlled
  // by a quadratic equation, so one has only to check whether
  // a certain discriminant is a square
  
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

  // nf_elem_is_square_fast finds most of the elements which are not squares
  // not to be squares, and it is never wrong when it does.
  // It can find an element to be a square when it is not, so
  // an element which passes it is checked by nf_elem_is_square,
  // which factors a polynomial of twice the degree of the field.
  is_square_fast = nf_elem_is_square_fast(disc, evs->nfs[ev_idx]);
  is_lift = is_square_fast && nf_elem_is_square(disc, evs->nfs[ev_idx]);
  assert(is_square_fast || !nf_elem_is_square(disc, evs->nfs[ev_idx]));

  // If is a lift setting it to a Yoshida lift temporarily, otherwise type G.
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
    nf_elem_clear(root, evs->nfs[ev_idx]);
  }
  
  nf_elem_clear(disc, evs->nfs[ev_idx]);
  for (k = 1; k < 3; k++) {
    nf_elem_clear(a[k-1], evs->nfs[ev_idx]);
  }
  
  return lift_type;
}

// Setting lift types for all eigenforms
void eigenvalues_set_lifts(eigenvalues_t evs, slong prec, slong c, const genus_t genus)
{
  slong cnt, i;
  fmpz_t p;
  aut_type lift_type;
  
  fmpz_init_set_si(p,2); // p = 2
  
  for (cnt = 0; cnt < prec; cnt++) {
    // skipping bad primes
    while (fmpz_divisible(genus->disc, p))
      fmpz_nextprime(p, p, true);
    for (i = 0; i < evs->num; i++) {
      // If this is an oldform, we don't try to further analyze
      if (evs->lift_type[i] == O)
	continue;
      // A form which is not a lift at one prime is not a lift,
      // so there is no need to consider it at further primes.
      if ((cnt > 0) && (evs->lift_type[i] == G))
	continue;
      lift_type = ev_is_lpoly_reducible(evs, i, fmpz_get_si(p), c, genus);
      // If this is the first pass, initialize according to result.
      // If we get a different lift type at different primes then it cannot be a lift.
      if (cnt == 0) 
	evs->lift_type[i] = lift_type;
      else
	if (evs->lift_type[i] != lift_type)
	  evs->lift_type[i] = G;
    }
    fmpz_nextprime(p, p, true);
  }

  fmpz_clear(p);
  
  return;
}

// Get the pivot for an eigenvector
slong ev_get_pivot(const nf_elem_t* evec, const nf_t K, const genus_t genus, slong cond)
{
  slong pivot, i, max_aut;

  max_aut = 1;
  pivot = -1;
  for (i = 0; i < genus->dims[cond]; i++) {
    // Look for a nonzero entry
    if (!nf_elem_is_zero(evec[i], K)) {
      // with maximal number of automorphisms
      if (max_aut < genus->num_auts[cond][i]) {
	max_aut = genus->num_auts[cond][i];
	pivot = i;
      }
    }
  }
  
  return pivot;
}
