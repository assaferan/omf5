#include <assert.h>

#include <time.h>

#include "genus.h"
#include "hash.h"
#include "isometry.h"
#include "nbr_data.h"
#include "neighbor.h"
#include "square_matrix.h"
#include "weight.h"

// Just update the matrix acting at each index
int process_isotropic_vector_isometries(neighbor_manager_t nbr_man,
					square_matrix_t** large_hecke,
					const genus_t genus,
					double* theta_time, double* isom_time,
					double* total_time, slong* num_isom, int gen_idx)
{
  int i, j;
  clock_t cputime;
  square_matrix_t nbr, aut_tr, orbit_sum;
  isometry_t s_nbr, s_inv;
  isometry_t hash_isom;
  bool non_singular;
  int vec_cmp;
  slong orb_size, stab_size;
  isometry_t* orbit_isom;
  vector_t* orbit;
  bool found;
  vector_t g_vec;
  
  assert(genus->genus_reps->num_stored > 0);

  orbit = (vector_t*)malloc(nbr_man->num_auts * sizeof(vector_t));
  orbit_isom = (isometry_t*)malloc(nbr_man->num_auts * sizeof(isometry_t));
  isometry_init(orbit_isom[0]);

  // right now we grow the orbit.
  // Could also compute on the fly and in the end divide by the stabilizer's size
  // weird that it isn't already reduced?
  vector_mod_p(nbr_man->iso_vec, nbr_man->p);
  // this might be unnecessary, but we want to be on the safe side for now
  normalize_mod_p(nbr_man->iso_vec, nbr_man->p);
  vector_set(orbit[0], nbr_man->iso_vec);
  // Applying the automorphism group to the isotropic vector
  stab_size = 0;
  orb_size = 1;
  for (i = 0; i < nbr_man->num_auts; i++) {
    square_matrix_mul_vec_left(g_vec, nbr_man->iso_vec, nbr_man->auts[i]);
    vector_mod_p(g_vec, nbr_man->p);
    normalize_mod_p(g_vec, nbr_man->p);
    vec_cmp = vector_cmp(g_vec, nbr_man->iso_vec);
    if (vec_cmp < 0) {
      free(orbit);
      free(orbit_isom);
      return 0;
    }
    if (vec_cmp == 0)
      stab_size++;
    if (vec_cmp > 0) {
      found = false;
      for (j = 0; j < orb_size; j++)
	if (vector_cmp(g_vec, orbit[j]) == 0) {
	  found = true;
	  break;
	}
      if (!found) {
	vector_set(orbit[orb_size], g_vec);
	square_matrix_transpose(aut_tr, nbr_man->auts[i]);
	isometry_init_set_square_matrix(orbit_isom[orb_size], aut_tr, 1);
	orb_size++;
      }
    }
  }
  assert(orb_size == nbr_man->num_auts / stab_size);
  (void)stab_size; // only used in the assertion
  
  non_singular = nbr_process_build_nb_and_isom(nbr, s_nbr, nbr_man);

  if (!non_singular) {
    free(orbit);
    free(orbit_isom);
    return 0;
  }
  
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
  
  cputime = clock();
  i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

  assert( (i >= 0) && (i <= genus->dims[0]) );

  // TODO - compute only the upper half, and complete using hermitian property
  // Determine which subspaces this representative contributes.

  // TODO - reduce the number of actual computations required here

  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
  for (j = 0; j < orb_size; j++)
    assert(isometry_is_isom(orbit_isom[j], genus->genus_reps->keys[gen_idx], genus->genus_reps->keys[gen_idx]));
  assert(isometry_is_isom(genus->isoms[gen_idx], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[gen_idx]));
  
  isometry_muleq_left(s_nbr, genus->isoms[gen_idx]);
  for (j = 0; j < orb_size; j++)
     isometry_muleq_left(orbit_isom[j], genus->isoms[gen_idx]);
  
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], nbr));
  for (j = 0; j < orb_size; j++)
     assert(isometry_is_isom(orbit_isom[j], genus->genus_reps->keys[0], genus->genus_reps->keys[gen_idx]));

  isometry_inv(s_inv, genus->isoms[gen_idx]);
  for (j = 0; j < orb_size; j++)
     isometry_muleq_right(orbit_isom[j], s_inv);
  
  for (j = 0; j < orb_size; j++)
    assert(isometry_is_isom(orbit_isom[j], genus->genus_reps->keys[0],
			    genus->genus_reps->keys[0]));
  
  assert(isometry_is_isom(hash_isom, nbr, genus->genus_reps->keys[i]));
  
  isometry_muleq_right(s_nbr, hash_isom);
  
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], genus->genus_reps->keys[i]));
  assert(isometry_is_isom(genus->isoms[i], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[i]));
  
  isometry_inv(s_inv, genus->isoms[i]);
  
  assert(isometry_is_isom(s_inv, genus->genus_reps->keys[i],
			  genus->genus_reps->keys[0]));
  
  isometry_muleq_right(s_nbr, s_inv);
  
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0],
			  genus->genus_reps->keys[0]));

  // !! TODO - handle denominators !!!
  square_matrix_init(orbit_sum);
  square_matrix_zero(orbit_sum);
  for (j = 0; j < orb_size; j++) {
    square_matrix_add(orbit_sum, orbit_sum, orbit_isom[j]->s);
  }
  square_matrix_muleq_right(orbit_sum, s_nbr->s);
  square_matrix_add(large_hecke[gen_idx][i],large_hecke[gen_idx][i], orbit_sum); 

  isometry_clear(s_nbr);
  isometry_clear(s_inv);
  isometry_clear(hash_isom);

  square_matrix_clear(nbr);

  free(orbit);
  free(orbit_isom);
  
  return orb_size;
}


slong hecke_col_isometries(square_matrix_t** large_hecke,
			   int p, int gen_idx, const genus_t genus)
{
  square_matrix_t Q;
  neighbor_manager_t nbr_man;
  slong lc;
  int i;
  slong num_isom;
  double theta_time, isom_time, total_time;

  num_isom = lc = 0;
  theta_time = isom_time = total_time = 0;

  lc = 0;
  
  square_matrix_set(Q,genus->genus_reps->keys[gen_idx]);

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL

  for (i = 0; i < p; i++) {
    nbr_process_init(nbr_man, Q, p, i, genus->isoms[gen_idx]);

    while (!(nbr_process_has_ended(nbr_man))) {
      lc += process_isotropic_vector_isometries(nbr_man, large_hecke, genus,
						&theta_time, &isom_time, &total_time,
						&num_isom, gen_idx);
      nbr_process_advance(nbr_man);
    }
    nbr_process_clear(nbr_man);
  }

#ifdef DEBUG
  printf("theta_time = %f, isom_time = %f, total_time = %f, num_isom = %ld / %ld \n",
	 theta_time/lc, isom_time/lc, total_time, num_isom, lc);
#endif // DEBUG

  return lc;
}


square_matrix_t** hecke_matrices_isometries(const genus_t genus, int p)
{
  square_matrix_t** hecke;
 
  int gen_idx;
  slong i, j;

  hecke = (square_matrix_t**)malloc(genus->dims[0] * sizeof(square_matrix_t*));

  for (i = 0; i < genus->dims[0]; i++) {
    hecke[i] = (square_matrix_t*)malloc(genus->dims[0] * sizeof(square_matrix_t));
    for (j = 0; j < genus->dims[0]; j++)
      square_matrix_zero(hecke[i][j]);
  }
  
  for (gen_idx = 0; gen_idx < genus->dims[0]; gen_idx++)
    hecke_col_isometries(hecke, p, gen_idx, genus);
  
  return hecke;
}

/* Algebraic modular forms of weight (a, b): see weight.h. */

#include <flint/fmpz_mat.h>
#include <flint/fmpz_vec.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <flint/ulong_extras.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_vec.h>
#include <flint/fmpz_mpoly.h>

#include "aut_grp.h"
#include "spinor.h"

// genus.c
int popcnt(W64 x);


// Fills exps (if not NULL) with the exponents of the monomials of degree a in QF_RANK variables,
// in a fixed order, and returns their number.
static slong weight_monomials(ulong* exps, int a)
{
  ulong e[QF_RANK];
  slong n, i, j;
  int rest;

  // the compositions of a into QF_RANK parts, by a counter over the first QF_RANK - 1 parts
  for (i = 0; i < QF_RANK; i++)
    e[i] = 0;
  n = 0;
  while (1) {
    rest = a;
    for (i = 0; i < QF_RANK - 1; i++)
      rest -= (int)e[i];
    if (rest >= 0) {
      if (exps != NULL) {
	for (j = 0; j < QF_RANK - 1; j++)
	  exps[QF_RANK*n + j] = e[j];
	exps[QF_RANK*n + QF_RANK - 1] = (ulong)rest;
      }
      n++;
    }
    // next counter value, each part at most a
    for (i = 0; i < QF_RANK - 1; i++) {
      if ((int)e[i] < a) {
	e[i]++;
	break;
      }
      e[i] = 0;
    }
    if (i == QF_RANK - 1)
      break;
  }
  return n;
}

static slong weight_monomial_index(const ulong* exps, slong n, const ulong* e)
{
  slong k, i;

  for (k = 0; k < n; k++) {
    for (i = 0; (i < QF_RANK) && (exps[QF_RANK*k + i] == e[i]); i++);
    if (i == QF_RANK)
      return k;
  }
  return -1;
}

// Sets B to a basis, as rows in reduced row echelon form, of the rational nullspace of A, and piv to
// the pivot of each row. Returns its dimension. B is initialized here.
static slong weight_nullspace_rref(fmpq_mat_t B, slong** piv, const fmpz_mat_t A)
{
  fmpz_mat_t X, Y;
  fmpz_t den;
  slong nullity, i, j, n;

  n = fmpz_mat_ncols(A);
  fmpz_mat_init(X, n, n);
  nullity = fmpz_mat_nullspace(X, A);
  // the reduced row echelon form, which is unique, by the multimodular algorithm
  fmpz_mat_init(Y, nullity, n);
  for (i = 0; i < nullity; i++)
    for (j = 0; j < n; j++)
      fmpz_swap(fmpz_mat_entry(Y, i, j), fmpz_mat_entry(X, j, i));
  fmpz_init(den);
  if (nullity > 0)
    fmpz_mat_rref_mul(Y, den, Y);
  fmpq_mat_init(B, nullity, n);
  for (i = 0; i < nullity; i++)
    for (j = 0; j < n; j++)
      fmpq_set_fmpz_frac(fmpq_mat_entry(B, i, j), fmpz_mat_entry(Y, i, j), den);
  fmpz_mat_clear(Y);
  fmpz_clear(den);
  *piv = (slong*)malloc((nullity > 0 ? nullity : 1) * sizeof(slong));
  for (i = 0; i < nullity; i++) {
    for (j = 0; fmpq_is_zero(fmpq_mat_entry(B, i, j)); j++);
    (*piv)[i] = j;
  }
  fmpz_mat_clear(X);
  return nullity;
}

// Sets the integral matrix A to the rows of the rational matrix M times the lcm of their denominators.
static void weight_clear_denominators(fmpz_mat_t A, const fmpq_mat_t M)
{
  fmpq_mat_get_fmpz_mat_rowwise(A, NULL, M);
  return;
}

// The monomials of bidegree (a, b) in x_1, ..., x_QF_RANK and y_1, ..., y_QF_RANK: the monomial
// kx * ny + ky is the product of the monomials kx of degree a in x and ky of degree b in y, of the
// orders of weight_monomials. Fills exps (if not NULL), WEIGHT_NVARS exponents for each, the x first.
static slong weight_bimonomials(ulong* exps, int a, int b)
{
  slong nx, ny, kx, ky, i;
  ulong *ex, *ey;

  if ((a < 0) || (b < 0))
    return 0;
  nx = weight_monomials(NULL, a);
  ny = weight_monomials(NULL, b);
  if (exps != NULL) {
    ex = (ulong*)malloc(QF_RANK * nx * sizeof(ulong));
    ey = (ulong*)malloc(QF_RANK * ny * sizeof(ulong));
    weight_monomials(ex, a);
    weight_monomials(ey, b);
    for (kx = 0; kx < nx; kx++)
      for (ky = 0; ky < ny; ky++)
	for (i = 0; i < QF_RANK; i++) {
	  exps[WEIGHT_NVARS*(kx*ny + ky) + i] = ex[QF_RANK*kx + i];
	  exps[WEIGHT_NVARS*(kx*ny + ky) + QF_RANK + i] = ey[QF_RANK*ky + i];
	}
    free(ex);
    free(ey);
  }
  return nx * ny;
}

// The monomials of a bidegree, with the index of a monomial of it.
typedef struct
{
  slong nx, ny;
  ulong *ex, *ey;
} weight_target_struct;

static void weight_target_init(weight_target_struct* t, int a, int b)
{
  t->nx = weight_monomials(NULL, a);
  t->ny = weight_monomials(NULL, b);
  t->ex = (ulong*)malloc(QF_RANK * (t->nx > 0 ? t->nx : 1) * sizeof(ulong));
  t->ey = (ulong*)malloc(QF_RANK * (t->ny > 0 ? t->ny : 1) * sizeof(ulong));
  weight_monomials(t->ex, a);
  weight_monomials(t->ey, b);
  return;
}

static void weight_target_clear(weight_target_struct* t)
{
  free(t->ex);
  free(t->ey);
  return;
}

static slong weight_target_index(const weight_target_struct* t, const ulong* e)
{
  slong kx, ky;

  kx = weight_monomial_index(t->ex, t->nx, e);
  ky = weight_monomial_index(t->ey, t->ny, e + QF_RANK);
  // checked in every build: a negative index would write outside the matrix of the conditions
  if ((kx < 0) || (ky < 0)) {
    fprintf(stderr, "Error: a monomial of the weight is not in its target bidegree\n");
    abort();
  }
  return kx * t->ny + ky;
}

// The number of monomials of degree d in QF_RANK = 5 variables, C(d + 4, 4), and of bidegree (a, b).
static double weight_num_monomials(double d)
{
  return (d < 0) ? 0 : (d + 1) * (d + 2) * (d + 3) * (d + 4) / 24;
}

static double weight_num_bimonomials(double a, double b)
{
  return weight_num_monomials(a) * weight_num_monomials(b);
}

double weight_rep_matrix_size(int a, int b)
{
  double nmons, nrows, dim, blocks;

  // in closed form, so that a large weight is rejected at once: the basis, dim x nmons, and the blocks
  // of the conditions, about 1/16 of the rows and of the columns each
  nmons = weight_num_bimonomials(a, b);
  nrows = ((a >= 2) ? weight_num_bimonomials(a - 2.0, b) : 0) +
    ((b >= 1) ? weight_num_bimonomials(a - 1.0, b - 1.0) + weight_num_bimonomials(a + 1.0, b - 1.0) : 0) +
    ((b >= 2) ? weight_num_bimonomials(a, b - 2.0) : 0);
  dim = (2.0*a + 3) * (2.0*b + 1) * (a - b + 1.0) * (a + b + 2.0) / 6;
  blocks = (nrows + nmons) * nmons / 256;
  return (dim * nmons > blocks) ? dim * nmons : blocks;
}

// A (n x n) is the matrix of the substitution x_r -> sum_j h_rj x_j on the n monomials of degree d of
// exps, of the order of weight_monomials: A[k][i] is the coefficient of the monomial i in the image of
// the monomial k.
static void weight_monomial_action(fmpz_mat_t A, const ulong* exps, slong n, int d, const fmpz_mat_t h)
{
  fmpz_mpoly_ctx_t ctx;
  fmpz_mpoly_struct lin[QF_RANK];
  fmpz_mpoly_t img, pw, var;
  fmpz_t c;
  ulong e[QF_RANK];
  slong k, r, j, t, i;

  fmpz_mpoly_ctx_init(ctx, QF_RANK, ORD_LEX);
  fmpz_mpoly_init(img, ctx);
  fmpz_mpoly_init(pw, ctx);
  fmpz_mpoly_init(var, ctx);
  fmpz_init(c);
  for (r = 0; r < QF_RANK; r++) {
    fmpz_mpoly_init(lin + r, ctx);
    for (j = 0; j < QF_RANK; j++) {
      fmpz_mpoly_gen(var, j, ctx);
      fmpz_mpoly_scalar_mul_fmpz(var, var, fmpz_mat_entry(h, r, j), ctx);
      fmpz_mpoly_add(lin + r, lin + r, var, ctx);
    }
  }
  fmpz_mat_zero(A);
  for (k = 0; k < n; k++) {
    fmpz_mpoly_one(img, ctx);
    for (r = 0; r < QF_RANK; r++) {
      if (exps[QF_RANK*k + r] == 0) continue;
      fmpz_mpoly_pow_ui(pw, lin + r, exps[QF_RANK*k + r], ctx);
      fmpz_mpoly_mul(img, img, pw, ctx);
    }
    for (t = 0; t < fmpz_mpoly_length(img, ctx); t++) {
      fmpz_mpoly_get_term_exp_ui(e, img, t, ctx);
      i = weight_monomial_index(exps, n, e);
      // checked in every build: the image of a monomial of degree d has degree d
      if (i < 0) {
	fprintf(stderr, "Error: the image of a monomial of degree %d is not of degree %d\n", d, d);
	abort();
      }
      fmpz_mpoly_get_term_coeff_fmpz(c, img, t, ctx);
      fmpz_set(fmpz_mat_entry(A, k, i), c);
    }
  }
  for (r = 0; r < QF_RANK; r++)
    fmpz_mpoly_clear(lin + r, ctx);
  fmpz_mpoly_clear(img, ctx);
  fmpz_mpoly_clear(pw, ctx);
  fmpz_mpoly_clear(var, ctx);
  fmpz_clear(c);
  fmpz_mpoly_ctx_clear(ctx);
  return;
}

// The images of the monomial of bidegree (a, b) with exponents mon under the operators which define the
// representation, in coordinates in which the inverse Gram matrix is diag(c): Delta_xx, Delta_xy, Delta_yy
// and E, with Delta_uv = sum_i c_i d/du_i d/dv_i and E = sum_i x_i d/dy_i. Sets rows[t] to the row of the
// t-th image, numbered as in weight_rep_init, and coeffs[t] to its coefficient; returns their number.
static slong weight_conditions(slong* rows, fmpz* coeffs, const ulong* mon, int a, int b, const fmpz* c,
			       const weight_target_struct* t_xx, const weight_target_struct* t_xy,
			       const weight_target_struct* t_yy, const weight_target_struct* t_e,
			       slong off_xy, slong off_yy, slong off_e)
{
  ulong e[WEIGHT_NVARS];
  slong i, r, n;

  n = 0;
  for (i = 0; i < QF_RANK; i++) {
    // Delta_xx
    if ((a >= 2) && (mon[i] >= 2)) {
      for (r = 0; r < WEIGHT_NVARS; r++)
	e[r] = mon[r];
      e[i] -= 2;
      rows[n] = weight_target_index(t_xx, e);
      fmpz_mul_ui(coeffs + n, c + i, mon[i] * (mon[i] - 1));
      n++;
    }
    // Delta_yy
    if ((b >= 2) && (mon[QF_RANK + i] >= 2)) {
      for (r = 0; r < WEIGHT_NVARS; r++)
	e[r] = mon[r];
      e[QF_RANK + i] -= 2;
      rows[n] = off_yy + weight_target_index(t_yy, e);
      fmpz_mul_ui(coeffs + n, c + i, mon[QF_RANK + i] * (mon[QF_RANK + i] - 1));
      n++;
    }
    if ((b < 1) || (mon[QF_RANK + i] == 0))
      continue;
    // Delta_xy
    if (mon[i] >= 1) {
      for (r = 0; r < WEIGHT_NVARS; r++)
	e[r] = mon[r];
      e[i]--;
      e[QF_RANK + i]--;
      rows[n] = off_xy + weight_target_index(t_xy, e);
      fmpz_mul_ui(coeffs + n, c + i, mon[i] * mon[QF_RANK + i]);
      n++;
    }
    // E
    for (r = 0; r < WEIGHT_NVARS; r++)
      e[r] = mon[r];
    e[i]++;
    e[QF_RANK + i]--;
    rows[n] = off_e + weight_target_index(t_e, e);
    fmpz_set_ui(coeffs + n, mon[QF_RANK + i]);
    n++;
  }
  return n;
}

void weight_rep_init(weight_rep_t W, int a, int b, const square_matrix_t q0)
{
  slong i, j, k, r, t, nrows, off_xy, off_yy, off_e, nc, nr, nullity, cnt, rank, nterms, dim_expected;
  slong *col_sig, *cols, *row_loc, *row_list, *rows;
  fmpz *c, *coeffs;
  fmpq_mat_t Q, M, Minv, G;
  fmpq_t mu, dot;
  fmpz_mat_t Mi, Ax, By, Nz, Fs, F2, Res, Blk, X, Rz;
  nmod_mat_t Rp;
  fmpz_t g;
  ulong prime;
  bool done;
  fmpz_t den;
  weight_target_struct t_xx, t_xy, t_yy, t_e;
  int sig;

  assert((0 <= b) && (b <= a));
  W->a = a;
  W->b = b;
  W->nmons = weight_bimonomials(NULL, a, b);
  W->exps = (ulong*)malloc(WEIGHT_NVARS * W->nmons * sizeof(ulong));
  weight_bimonomials(W->exps, a, b);
  W->nx = weight_monomials(NULL, a);
  W->ny = weight_monomials(NULL, b);
  W->ex = (ulong*)malloc(QF_RANK * W->nx * sizeof(ulong));
  W->ey = (ulong*)malloc(QF_RANK * W->ny * sizeof(ulong));
  weight_monomials(W->ex, a);
  weight_monomials(W->ey, b);
  dim_expected = (2*a + 3) * (2*b + 1) * (a - b + 1) * (a + b + 2) / 6;

  // The polynomials f of bidegree (a, b) which are harmonic for q0^(-1) in x, in y and in both, that is
  // Delta_xx f = Delta_xy f = Delta_yy f = 0 with Delta_uv = sum_{i,j} qinv_ij d/du_i d/dv_j, and are
  // highest weight vectors for GL_2, E f = 0 with E = sum_i x_i d/dy_i. They form the irreducible
  // representation of highest weight (a, b) of O(q0) (Howe duality for O(5) x GL_2 on the harmonic
  // polynomials).
  //
  // They are found in the coordinates x = M x', y = M y', with M^T q0 M = D diagonal (Gram-Schmidt over
  // the rationals): f(M x', M y') satisfies the same conditions with qinv replaced by D^(-1), and E is
  // unchanged. There each operator keeps the parities of e_x_i + e_y_i, i = 1, ..., QF_RANK, of the
  // exponents of a monomial, so the conditions split into blocks, one for each parity vector.
  fmpq_mat_init(Q, QF_RANK, QF_RANK);
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      fmpq_set_si(fmpq_mat_entry(Q, i, j), q0[i][j], 1);
  fmpq_mat_init(M, QF_RANK, QF_RANK);
  fmpq_mat_one(M);
  fmpq_mat_init(G, QF_RANK, QF_RANK);
  fmpq_init(mu);
  fmpq_init(dot);
  // the column k of M is e_k minus its projections on the previous columns, for the form q0
  for (k = 0; k < QF_RANK; k++)
    for (j = 0; j < k; j++) {
      // mu = <e_k, m_j> / <m_j, m_j>, with <u, v> = u^T q0 v
      fmpq_zero(mu);
      fmpq_zero(dot);
      for (i = 0; i < QF_RANK; i++)
	fmpq_addmul(mu, fmpq_mat_entry(Q, k, i), fmpq_mat_entry(M, i, j));
      fmpq_mat_mul(G, Q, M);
      for (i = 0; i < QF_RANK; i++)
	fmpq_addmul(dot, fmpq_mat_entry(M, i, j), fmpq_mat_entry(G, i, j));
      fmpq_div(mu, mu, dot);
      for (i = 0; i < QF_RANK; i++)
	fmpq_submul(fmpq_mat_entry(M, i, k), mu, fmpq_mat_entry(M, i, j));
    }
  // D = M^T q0 M, and c_i proportional to 1 / D_ii, integral
  fmpq_mat_mul(G, Q, M);
  fmpq_mat_transpose(Q, M);
  fmpq_mat_mul(G, Q, G);
  c = _fmpz_vec_init(QF_RANK);
  for (i = 0; i < QF_RANK; i++) {
    // checked in every build: q0 is positive definite, and D diagonal
    for (j = 0; j < QF_RANK; j++)
      if ((i != j) && !fmpq_is_zero(fmpq_mat_entry(G, i, j))) {
	fprintf(stderr, "Error: the Gram-Schmidt basis of the weight is not orthogonal\n");
	abort();
      }
    if (fmpq_sgn(fmpq_mat_entry(G, i, i)) <= 0) {
      fprintf(stderr, "Error: the quadratic form of the weight is not positive definite\n");
      abort();
    }
  }
  // 1 / D_ii = den_i / num_i; times the product of the numerators
  for (i = 0; i < QF_RANK; i++) {
    fmpz_set(c + i, fmpq_denref(fmpq_mat_entry(G, i, i)));
    for (j = 0; j < QF_RANK; j++)
      if (j != i)
	fmpz_mul(c + i, c + i, fmpq_numref(fmpq_mat_entry(G, j, j)));
  }

  // the rows of the conditions, numbered as Delta_xx f of bidegree (a - 2, b), then Delta_xy f of
  // (a - 1, b - 1), Delta_yy f of (a, b - 2) and E f of (a + 1, b - 1)
  weight_target_init(&t_xx, a - 2 >= 0 ? a - 2 : 0, b);
  weight_target_init(&t_xy, a - 1 >= 0 ? a - 1 : 0, b - 1 >= 0 ? b - 1 : 0);
  weight_target_init(&t_yy, a, b - 2 >= 0 ? b - 2 : 0);
  weight_target_init(&t_e, a + 1, b - 1 >= 0 ? b - 1 : 0);
  off_xy = (a >= 2) ? t_xx.nx * t_xx.ny : 0;
  off_yy = off_xy + ((b >= 1) ? t_xy.nx * t_xy.ny : 0);
  off_e = off_yy + ((b >= 2) ? t_yy.nx * t_yy.ny : 0);
  nrows = off_e + ((b >= 1) ? t_e.nx * t_e.ny : 0);

  col_sig = (slong*)malloc((W->nmons > 0 ? W->nmons : 1) * sizeof(slong));
  for (k = 0; k < W->nmons; k++) {
    sig = 0;
    for (i = 0; i < QF_RANK; i++)
      sig |= (int)((W->exps[WEIGHT_NVARS*k + i] + W->exps[WEIGHT_NVARS*k + QF_RANK + i]) & 1) << i;
    col_sig[k] = sig;
  }
  cols = (slong*)malloc((W->nmons > 0 ? W->nmons : 1) * sizeof(slong));
  row_loc = (slong*)malloc((nrows > 0 ? nrows : 1) * sizeof(slong));
  row_list = (slong*)malloc((nrows > 0 ? nrows : 1) * sizeof(slong));
  for (r = 0; r < nrows; r++)
    row_loc[r] = -1;
  rows = (slong*)malloc(4 * QF_RANK * sizeof(slong));
  coeffs = _fmpz_vec_init(4 * QF_RANK);

  // the kernel, in the coordinates x', y', as the rows of Nz
  fmpz_mat_init(Nz, dim_expected, W->nmons);
  cnt = 0;
  for (sig = 0; sig < (1 << QF_RANK); sig++) {
    nc = 0;
    for (k = 0; k < W->nmons; k++)
      if (col_sig[k] == sig)
	cols[nc++] = k;
    if (nc == 0)
      continue;
    nr = 0;
    for (j = 0; j < nc; j++) {
      nterms = weight_conditions(rows, coeffs, W->exps + WEIGHT_NVARS*cols[j], a, b, c,
				 &t_xx, &t_xy, &t_yy, &t_e, off_xy, off_yy, off_e);
      for (t = 0; t < nterms; t++)
	if (row_loc[rows[t]] < 0) {
	  row_loc[rows[t]] = nr;
	  row_list[nr++] = rows[t];
	}
    }
    fmpz_mat_init(X, nc, nc);
    if (nr == 0) {
      nullity = nc;
      fmpz_mat_one(X);
    }
    else {
      fmpz_mat_init(Blk, nr, nc);
      for (j = 0; j < nc; j++) {
	nterms = weight_conditions(rows, coeffs, W->exps + WEIGHT_NVARS*cols[j], a, b, c,
				   &t_xx, &t_xy, &t_yy, &t_e, off_xy, off_yy, off_e);
	for (t = 0; t < nterms; t++)
	  fmpz_add(fmpz_mat_entry(Blk, row_loc[rows[t]], j), fmpz_mat_entry(Blk, row_loc[rows[t]], j), coeffs + t);
      }
      nullity = fmpz_mat_nullspace(X, Blk);
      fmpz_mat_clear(Blk);
    }
    for (r = 0; r < nr; r++)
      row_loc[row_list[r]] = -1;
    // checked in every build: at most the dimension of the representation
    if (cnt + nullity > dim_expected) {
      fprintf(stderr, "Error: the representation of weight (%d, %d) has dimension more than %ld\n", a, b,
	      (long)dim_expected);
      abort();
    }
    for (i = 0; i < nullity; i++) {
      for (j = 0; j < nc; j++)
	fmpz_set(fmpz_mat_entry(Nz, cnt, cols[j]), fmpz_mat_entry(X, j, i));
      cnt++;
    }
    fmpz_mat_clear(X);
  }
  // the dimension of the representation of highest weight (a, b) of SO(5), checked in every build
  if (cnt != dim_expected) {
    fprintf(stderr, "Error: the representation of weight (%d, %d) has dimension %ld, not %ld\n", a, b,
	    (long)cnt, (long)dim_expected);
    abort();
  }
  W->dim = cnt;

  // back to x, y: f(x, y) = f'(M^(-1) x, M^(-1) y). With Mi = den M^(-1) integral, and Ax, By the matrices
  // of Mi on the monomials of degree a and b, the coefficient of x^i y^j is, up to den^(a + b),
  // sum_{kx, ky} f'[kx, ky] Ax[kx][i] By[ky][j], computed as Ax^T (f' By) for all the rows at once.
  fmpq_mat_init(Minv, QF_RANK, QF_RANK);
  fmpq_mat_inv(Minv, M);
  fmpz_mat_init(Mi, QF_RANK, QF_RANK);
  fmpz_init(den);
  fmpq_mat_get_fmpz_mat_matwise(Mi, den, Minv);
  fmpz_mat_init(Ax, W->nx, W->nx);
  fmpz_mat_init(By, W->ny, W->ny);
  weight_monomial_action(Ax, W->ex, W->nx, a, Mi);
  weight_monomial_action(By, W->ey, W->ny, b, Mi);
  fmpz_mat_init(Fs, W->dim * W->nx, W->ny);
  for (i = 0; i < W->dim; i++)
    for (k = 0; k < W->nmons; k++)
      fmpz_set(fmpz_mat_entry(Fs, i * W->nx + k / W->ny, k % W->ny), fmpz_mat_entry(Nz, i, k));
  fmpz_mat_clear(Nz);
  fmpz_mat_init(F2, W->dim * W->nx, W->ny);
  fmpz_mat_mul(F2, Fs, By);
  fmpz_mat_clear(Fs);
  // F2 rearranged as nx x (dim ny), then multiplied by Ax^T
  fmpz_mat_init(Fs, W->nx, (W->dim > 0 ? W->dim : 1) * W->ny);
  for (i = 0; i < W->dim; i++)
    for (k = 0; k < W->nx; k++)
      for (j = 0; j < W->ny; j++)
	fmpz_swap(fmpz_mat_entry(Fs, k, i * W->ny + j), fmpz_mat_entry(F2, i * W->nx + k, j));
  fmpz_mat_clear(F2);
  fmpz_mat_transpose(Ax, Ax);
  fmpz_mat_init(Res, W->nx, (W->dim > 0 ? W->dim : 1) * W->ny);
  fmpz_mat_mul(Res, Ax, Fs);
  fmpz_mat_clear(Fs);
  fmpz_mat_init(Rz, W->dim, W->nmons);
  for (i = 0; i < W->dim; i++)
    for (k = 0; k < W->nx; k++)
      for (j = 0; j < W->ny; j++)
	fmpz_swap(fmpz_mat_entry(Rz, i, k * W->ny + j), fmpz_mat_entry(Res, k, i * W->ny + j));
  fmpz_mat_clear(Res);

  // the basis in reduced row echelon form, which is unique. Its pivots are found modulo a prime, and the
  // rows are then the solution X of P X = Rz, with P the columns of Rz at the pivots. X is checked to be
  // in reduced row echelon form, which proves it is the basis; otherwise, another prime is tried.
  for (i = 0; i < W->dim; i++) {
    _fmpz_vec_content(den, fmpz_mat_entry(Rz, i, 0), W->nmons);
    if (!fmpz_is_zero(den) && !fmpz_is_one(den))
      _fmpz_vec_scalar_divexact_fmpz(fmpz_mat_entry(Rz, i, 0), fmpz_mat_entry(Rz, i, 0), W->nmons, den);
  }
  W->pivots = (slong*)malloc((W->dim > 0 ? W->dim : 1) * sizeof(slong));
  fmpz_mat_init(X, W->dim, W->nmons);
  fmpz_mat_init(Blk, W->dim, W->dim);
  prime = UWORD(1) << 60;
  done = (W->dim == 0);
  while (!done) {
    prime = n_nextprime(prime, 1);
    nmod_mat_init(Rp, W->dim, W->nmons, prime);
    fmpz_mat_get_nmod_mat(Rp, Rz);
    rank = nmod_mat_rref(Rp);
    if (rank == W->dim) {
      for (i = 0; i < W->dim; i++) {
	for (k = 0; nmod_mat_entry(Rp, i, k) == 0; k++);
	W->pivots[i] = k;
      }
      for (i = 0; i < W->dim; i++)
	for (j = 0; j < W->dim; j++)
	  fmpz_set(fmpz_mat_entry(Blk, i, j), fmpz_mat_entry(Rz, i, W->pivots[j]));
      if (fmpz_mat_solve(X, den, Blk, Rz)) {
	done = true;
	for (i = 0; (i < W->dim) && done; i++) {
	  for (k = 0; (k < W->pivots[i]) && done; k++)
	    done = fmpz_is_zero(fmpz_mat_entry(X, i, k));
	  done = done && fmpz_equal(fmpz_mat_entry(X, i, W->pivots[i]), den);
	  for (j = 0; (j < W->dim) && done; j++)
	    done = (j == i) || fmpz_is_zero(fmpz_mat_entry(X, i, W->pivots[j]));
	}
      }
    }
    nmod_mat_clear(Rp);
  }
  fmpz_mat_clear(Blk);
  fmpz_mat_swap(Rz, X);
  fmpz_mat_clear(X);
  // the row i of the basis is the row i of Rz divided by den; with g the gcd of den and the content of
  // the row, its denominator is den / g and its numerators the row divided by g
  W->basis_den = _fmpz_vec_init(W->dim > 0 ? W->dim : 1);
  fmpz_init(g);
  for (i = 0; i < W->dim; i++) {
    _fmpz_vec_content(g, fmpz_mat_entry(Rz, i, 0), W->nmons);
    fmpz_gcd(g, g, den);
    fmpz_divexact(W->basis_den + i, den, g);
    _fmpz_vec_scalar_divexact_fmpz(fmpz_mat_entry(Rz, i, 0), fmpz_mat_entry(Rz, i, 0), W->nmons, g);
  }
  fmpz_clear(g);
  // the basis for weight_rep_action, the row i as the columns i * ny, ..., i * ny + ny - 1
  fmpz_mat_init(W->basis_z, W->nx, (W->dim > 0 ? W->dim : 1) * W->ny);
  for (i = 0; i < W->dim; i++)
    for (k = 0; k < W->nmons; k++)
      fmpz_swap(fmpz_mat_entry(W->basis_z, k / W->ny, i * W->ny + k % W->ny), fmpz_mat_entry(Rz, i, k));

  fmpz_mat_clear(Rz);
  fmpz_mat_clear(Ax);
  fmpz_mat_clear(By);
  fmpz_mat_clear(Mi);
  fmpz_clear(den);
  fmpq_mat_clear(Minv);
  fmpq_mat_clear(Q);
  fmpq_mat_clear(M);
  fmpq_mat_clear(G);
  fmpq_clear(mu);
  fmpq_clear(dot);
  _fmpz_vec_clear(c, QF_RANK);
  _fmpz_vec_clear(coeffs, 4 * QF_RANK);
  free(rows);
  free(col_sig);
  free(cols);
  free(row_loc);
  free(row_list);
  weight_target_clear(&t_xx);
  weight_target_clear(&t_xy);
  weight_target_clear(&t_yy);
  weight_target_clear(&t_e);
  return;
}

void weight_rep_clear(weight_rep_t W)
{
  _fmpz_vec_clear(W->basis_den, W->dim > 0 ? W->dim : 1);
  fmpz_mat_clear(W->basis_z);
  free(W->ex);
  free(W->ey);
  free(W->pivots);
  free(W->exps);
  return;
}

void weight_rep_action(fmpq_mat_t R, const weight_rep_t W, const isometry_t s)
{
  fmpz_mat_t A, B, Ax, Z, Zx, Bp, P, S, Sinv;
  fmpz_t det, scale, sum;
  slong k, m, kx, ky, py, npx, cnt, j;
  slong *xrow, *xpos, *piv_of;
  int sign;

  // det(s) = det(S) / denom^QF_RANK, which is 1 or -1; the factor det(s)^(a + b)
  fmpz_mat_init_set_square_matrix(S, s->s);
  fmpz_init(det);
  fmpz_mat_det(det, S);
  sign = ((fmpz_sgn(det) < 0) != (s->denom < 0)) ? -1 : 1;
  if (((W->a + W->b) % 2) == 0)
    sign = 1;

  // f(h x, h y) for h = s^(-1) = s_inv / inv_denom. With f = sum_{kx, ky} f[kx, ky] x^kx y^ky and A, B
  // the matrices of s_inv on the monomials of degree a and b, the coefficient of x^px y^py in
  // f(s_inv x, s_inv y) is sum_ky B[ky][py] sum_kx A[kx][px] f[kx, ky]. The coordinates of f(h x, h y)
  // are its coefficients at the pivots (px, py), divided by inv_denom^(a + b). The inner sums are one
  // product Z = Ax^T basis_z, Ax the columns of A at the monomials px of the pivots.
  fmpz_mat_init_set_square_matrix(Sinv, s->s_inv);
  fmpz_mat_init(A, W->nx, W->nx);
  fmpz_mat_init(B, W->ny, W->ny);
  weight_monomial_action(A, W->ex, W->nx, W->a, Sinv);
  weight_monomial_action(B, W->ey, W->ny, W->b, Sinv);
  fmpz_mat_clear(Sinv);
  xrow = (slong*)malloc(W->nx * sizeof(slong));
  xpos = (slong*)malloc((W->dim > 0 ? W->dim : 1) * sizeof(slong));
  for (kx = 0; kx < W->nx; kx++)
    xrow[kx] = -1;
  npx = 0;
  for (m = 0; m < W->dim; m++) {
    kx = W->pivots[m] / W->ny;
    if (xrow[kx] < 0)
      xrow[kx] = npx++;
    xpos[m] = xrow[kx];
  }
  fmpz_mat_init(Ax, (npx > 0) ? npx : 1, W->nx);
  for (kx = 0; kx < W->nx; kx++)
    if (xrow[kx] >= 0)
      for (k = 0; k < W->nx; k++)
	fmpz_set(fmpz_mat_entry(Ax, xrow[kx], k), fmpz_mat_entry(A, k, kx));
  fmpz_mat_init(Z, (npx > 0) ? npx : 1, (W->dim > 0 ? W->dim : 1) * W->ny);
  if (W->dim > 0)
    fmpz_mat_mul(Z, Ax, W->basis_z);

  // then, for each monomial px of the pivots, the sums over ky as one product: the row of Z at px as a
  // dim x ny matrix, times the columns of B at the monomials py of the pivots (px, py)
  fmpz_init(scale);
  fmpz_init(sum);
  fmpz_set_si(scale, s->inv_denom);
  fmpz_pow_ui(scale, scale, (ulong)(W->a + W->b));
  if (sign < 0)
    fmpz_neg(scale, scale);
  piv_of = (slong*)malloc((W->dim > 0 ? W->dim : 1) * sizeof(slong));
  fmpz_mat_init(Zx, (W->dim > 0 ? W->dim : 1), W->ny);
  for (kx = 0; kx < npx; kx++) {
    cnt = 0;
    for (m = 0; m < W->dim; m++)
      if (xpos[m] == kx)
	piv_of[cnt++] = m;
    fmpz_mat_init(Bp, W->ny, cnt);
    for (j = 0; j < cnt; j++) {
      py = W->pivots[piv_of[j]] % W->ny;
      for (ky = 0; ky < W->ny; ky++)
	fmpz_set(fmpz_mat_entry(Bp, ky, j), fmpz_mat_entry(B, ky, py));
    }
    for (k = 0; k < W->dim; k++)
      for (ky = 0; ky < W->ny; ky++)
	fmpz_set(fmpz_mat_entry(Zx, k, ky), fmpz_mat_entry(Z, kx, k * W->ny + ky));
    fmpz_mat_init(P, W->dim, cnt);
    fmpz_mat_mul(P, Zx, Bp);
    for (k = 0; k < W->dim; k++) {
      fmpz_mul(det, W->basis_den + k, scale);
      for (j = 0; j < cnt; j++)
	fmpq_set_fmpz_frac(fmpq_mat_entry(R, piv_of[j], k), fmpz_mat_entry(P, k, j), det);
    }
    fmpz_mat_clear(P);
    fmpz_mat_clear(Bp);
  }
  fmpz_mat_clear(Zx);
  free(piv_of);

  free(xrow);
  free(xpos);
  fmpz_mat_clear(A);
  fmpz_mat_clear(B);
  fmpz_mat_clear(Ax);
  fmpz_mat_clear(Z);
  fmpz_mat_clear(S);
  fmpz_clear(det);
  fmpz_clear(scale);
  fmpz_clear(sum);
  return;
}

// The value of the character of the conductor index c at the spinor norm bits.
static int weight_char(W64 bits, slong c)
{
  return (popcnt(bits & (W64)c) & 1) ? -1 : 1;
}

void weight_space_init(weight_space_t S, const genus_t genus, const weight_rep_t W)
{
  slong c, i, g, r, col, h;
  aut_grp_t grp;
  isometry_t s;
  fmpq_mat_t* R;
  W64* bits;
  fmpq_mat_t A;
  fmpz_mat_t Z;
  slong ngens;

  h = genus->genus_reps->num_stored;
  S->num_conductors = genus->num_conductors;
  S->num_classes = h;
  S->H = (fmpq_mat_struct**)malloc(S->num_conductors * sizeof(fmpq_mat_struct*));
  S->pivots = (slong***)malloc(S->num_conductors * sizeof(slong**));
  S->offsets = (slong**)malloc(S->num_conductors * sizeof(slong*));
  S->dims = (slong*)malloc(S->num_conductors * sizeof(slong));
  for (c = 0; c < S->num_conductors; c++) {
    S->H[c] = (fmpq_mat_struct*)malloc(h * sizeof(fmpq_mat_struct));
    S->pivots[c] = (slong**)malloc(h * sizeof(slong*));
    S->offsets[c] = (slong*)malloc(h * sizeof(slong));
    S->dims[c] = 0;
  }

  for (i = 0; i < h; i++) {
    // the automorphisms of L_i, taken to the first lattice as in ignore_set
    aut_grp_init_square_matrix(grp, genus->genus_reps->keys[i]);
    ngens = grp->num_gens;
    R = (fmpq_mat_t*)malloc((ngens > 0 ? ngens : 1) * sizeof(fmpq_mat_t));
    bits = (W64*)malloc((ngens > 0 ? ngens : 1) * sizeof(W64));
    for (g = 0; g < ngens; g++) {
      isometry_init(s);
      isometry_inv(s, genus->isoms[i]);
      isometry_mul_mateq_left(s, grp->gens[g], 1);
      isometry_muleq_left(s, genus->isoms[i]);
      assert(isometry_is_isom(s, genus->genus_reps->keys[0], genus->genus_reps->keys[0]));
      fmpq_mat_init(R[g], W->dim, W->dim);
      weight_rep_action(R[g], W, s);
      bits[g] = spinor_norm_isom(genus->spinor, s);
      isometry_clear(s);
    }
    for (c = 0; c < S->num_conductors; c++) {
      // the vectors w with chi_c(g) R(g) w = w for every generator g
      fmpq_mat_init(A, (ngens > 0 ? ngens : 1) * W->dim, W->dim);
      for (g = 0; g < ngens; g++)
	for (r = 0; r < W->dim; r++)
	  for (col = 0; col < W->dim; col++) {
	    fmpq_set(fmpq_mat_entry(A, g*W->dim + r, col), fmpq_mat_entry(R[g], r, col));
	    if (weight_char(bits[g], c) < 0)
	      fmpq_neg(fmpq_mat_entry(A, g*W->dim + r, col), fmpq_mat_entry(A, g*W->dim + r, col));
	    if (r == col)
	      fmpq_sub_si(fmpq_mat_entry(A, g*W->dim + r, col), fmpq_mat_entry(A, g*W->dim + r, col), 1);
	  }
      fmpz_mat_init(Z, fmpq_mat_nrows(A), W->dim);
      weight_clear_denominators(Z, A);
      weight_nullspace_rref(S->H[c] + i, &(S->pivots[c][i]), Z);
      S->offsets[c][i] = S->dims[c];
      S->dims[c] += fmpq_mat_nrows(S->H[c] + i);
      fmpz_mat_clear(Z);
      fmpq_mat_clear(A);
    }
    for (g = 0; g < ngens; g++)
      fmpq_mat_clear(R[g]);
    free(R);
    free(bits);
    aut_grp_clear(grp);
  }

  return;
}

void weight_space_clear(weight_space_t S)
{
  slong c, i;

  for (c = 0; c < S->num_conductors; c++) {
    for (i = 0; i < S->num_classes; i++) {
      fmpq_mat_clear(S->H[c] + i);
      free(S->pivots[c][i]);
    }
    free(S->H[c]);
    free(S->pivots[c]);
    free(S->offsets[c]);
  }
  free(S->H);
  free(S->pivots);
  free(S->offsets);
  free(S->dims);
  return;
}

// Adds the image of the basis vectors of H[c][i] under M (in the coordinates of W) to the columns of
// the block of rows j and columns i of T, in the coordinates of H[c][j]. Exits with an error if an
// image is not in H[c][j]: the sum over the neighbors has to be invariant under Aut(L_j).
static void weight_add_block(fmpq_mat_t T, const weight_space_t S, slong c, slong j, slong i, const fmpq_mat_t M)
{
  slong k, m, dW, ni, nj;
  fmpq_mat_t Ht, U, Up, Hjt, V;
  const fmpq_mat_struct* Hi = S->H[c] + i;
  const fmpq_mat_struct* Hj = S->H[c] + j;

  // the images of the basis vectors of H[c][i], the columns of U = M Hi^T
  dW = fmpq_mat_ncols(Hi);
  ni = fmpq_mat_nrows(Hi);
  nj = fmpq_mat_nrows(Hj);
  fmpq_mat_init(Ht, dW, ni);
  fmpq_mat_transpose(Ht, Hi);
  fmpq_mat_init(U, dW, ni);
  fmpq_mat_mul(U, M, Ht);
  // their coordinates in H[c][j] are their entries at the pivots, Up; check that they lie there,
  // U = Hj^T Up
  fmpq_mat_init(Up, nj, ni);
  for (m = 0; m < nj; m++)
    for (k = 0; k < ni; k++)
      fmpq_set(fmpq_mat_entry(Up, m, k), fmpq_mat_entry(U, S->pivots[c][j][m], k));
  fmpq_mat_init(Hjt, dW, nj);
  fmpq_mat_transpose(Hjt, Hj);
  fmpq_mat_init(V, dW, ni);
  fmpq_mat_mul(V, Hjt, Up);
  if (!fmpq_mat_equal(U, V)) {
    fprintf(stderr, "Error: at the conductor index %ld, the image of a form under a Hecke operator is not invariant under the automorphisms of the lattice\n", (long)c);
    exit(EXIT_FAILURE);
  }
  for (m = 0; m < nj; m++)
    for (k = 0; k < ni; k++)
      fmpq_add(fmpq_mat_entry(T, S->offsets[c][j] + m, S->offsets[c][i] + k),
	       fmpq_mat_entry(T, S->offsets[c][j] + m, S->offsets[c][i] + k), fmpq_mat_entry(Up, m, k));
  fmpq_mat_clear(Ht);
  fmpq_mat_clear(U);
  fmpq_mat_clear(Up);
  fmpq_mat_clear(Hjt);
  fmpq_mat_clear(V);
  return;
}

// The sums acc[i][b] for the lattice L_j, by orbits of the automorphisms of L_j on the isotropic
// lines, as process_isotropic_vector_all_conductors does: for an automorphism g of L_j the neighbor at
// g v is g N_v, whose isometry to the first lattice is t s, for t = isoms[j] g isoms[j]^(-1) and s that
// of N_v. So rho(t s) = rho(t) rho(s), and the orbit adds (sum of rho(t), by the spinor bits of t) times
// rho(s), with the spinor bits of t s those of t xor those of s. Only the representative of an orbit,
// the smallest vector in it, is tested for isometry.
static void weight_sums_orbits(fmpq_mat_t** acc, const genus_t genus, const weight_rep_t W, slong j, int p,
			       const fmpz_t pa)
{
  neighbor_manager_t nbr_man;
  fmpq_mat_t* Raut;
  W64* bits_aut;
  fmpq_mat_t* orbsum;
  fmpq_mat_t Rs, prod;
  slong* orbit_aut;
  vector_t* orbit;
  square_matrix_t* auts;
  vector_t g_vec;
  square_matrix_t nbr, aut_tr;
  isometry_t s_nbr, hash_isom, s_inv, t;
  slong k, e, b, nb, num_auts, orb_size, i;
  int line, vec_cmp;
  slong num_isom;
  double theta_time, isom_time;
  bool found, has_auts;
  W64 bits_s;

  nb = (slong)1 << genus->spinor->num_primes;
  theta_time = isom_time = 0;
  num_isom = 0;
  has_auts = false;
  Raut = NULL;
  bits_aut = NULL;
  num_auts = 0;
  orbsum = (fmpq_mat_t*)malloc(nb * sizeof(fmpq_mat_t));
  for (b = 0; b < nb; b++)
    fmpq_mat_init(orbsum[b], W->dim, W->dim);
  fmpq_mat_init(Rs, W->dim, W->dim);
  fmpq_mat_init(prod, W->dim, W->dim);

  for (line = 0; line < p; line++) {
    nbr_process_init(nbr_man, genus->genus_reps->keys[j], p, line, genus->isoms[j]);
    if (!has_auts) {
      // rho of every automorphism of L_j, taken to the first lattice
      num_auts = nbr_man->num_auts;
      Raut = (fmpq_mat_t*)malloc(num_auts * sizeof(fmpq_mat_t));
      bits_aut = (W64*)malloc(num_auts * sizeof(W64));
      isometry_inv(s_inv, genus->isoms[j]);
      for (k = 0; k < num_auts; k++) {
	square_matrix_transpose(aut_tr, nbr_man->auts[k]);
	isometry_init_set_square_matrix(t, aut_tr, 1);
	isometry_muleq_left(t, genus->isoms[j]);
	isometry_muleq_right(t, s_inv);
	assert(isometry_is_isom(t, genus->genus_reps->keys[0], genus->genus_reps->keys[0]));
	fmpq_mat_init(Raut[k], W->dim, W->dim);
	weight_rep_action(Raut[k], W, t);
	bits_aut[k] = spinor_norm_isom(genus->spinor, t);
	isometry_clear(t);
      }
      orbit = (vector_t*)malloc((num_auts + 1) * sizeof(vector_t));
      orbit_aut = (slong*)malloc((num_auts + 1) * sizeof(slong));
      // nbr_process_init finds the automorphisms again for every line, not always in the same
      // order, so the orbits are found with this list, whose k-th element is that of Raut[k]
      auts = (square_matrix_t*)malloc(num_auts * sizeof(square_matrix_t));
      for (k = 0; k < num_auts; k++)
	square_matrix_set(auts[k], nbr_man->auts[k]);
      has_auts = true;
    }
    assert(nbr_man->num_auts == num_auts);
    while (!nbr_process_has_ended(nbr_man)) {
      // the orbit of the isotropic vector, which is skipped unless it is the smallest in it
      vector_mod_p(nbr_man->iso_vec, nbr_man->p);
      normalize_mod_p(nbr_man->iso_vec, nbr_man->p);
      vector_set(orbit[0], nbr_man->iso_vec);
      orbit_aut[0] = -1;
      orb_size = 1;
      vec_cmp = 1;
      for (k = 0; k < num_auts; k++) {
	square_matrix_mul_vec_left(g_vec, nbr_man->iso_vec, auts[k]);
	vector_mod_p(g_vec, nbr_man->p);
	normalize_mod_p(g_vec, nbr_man->p);
	vec_cmp = vector_cmp(g_vec, nbr_man->iso_vec);
	if (vec_cmp < 0)
	  break;
	if (vec_cmp > 0) {
	  found = false;
	  for (e = 0; (e < orb_size) && !found; e++)
	    found = (vector_cmp(g_vec, orbit[e]) == 0);
	  if (!found) {
	    vector_set(orbit[orb_size], g_vec);
	    orbit_aut[orb_size] = k;
	    orb_size++;
	  }
	}
      }
      if ((vec_cmp >= 0) && nbr_process_build_nb_and_isom(nbr, s_nbr, nbr_man)) {
	i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, &theta_time, &isom_time, &num_isom);
	assert((0 <= i) && (i < genus->genus_reps->num_stored));
	isometry_muleq_left(s_nbr, genus->isoms[j]);
	isometry_muleq_right(s_nbr, hash_isom);
	isometry_inv(t, genus->isoms[i]);
	isometry_muleq_right(s_nbr, t);
	assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], genus->genus_reps->keys[0]));
	bits_s = spinor_norm_isom(genus->spinor, s_nbr);
	weight_rep_action(Rs, W, s_nbr);
	fmpq_mat_scalar_mul_fmpz(Rs, Rs, pa);
	for (b = 0; b < nb; b++)
	  fmpq_mat_zero(orbsum[b]);
	for (e = 0; e < orb_size; e++) {
	  if (orbit_aut[e] < 0) {
	    for (k = 0; k < W->dim; k++)
	      fmpq_add_si(fmpq_mat_entry(orbsum[0], k, k), fmpq_mat_entry(orbsum[0], k, k), 1);
	  }
	  else
	    fmpq_mat_add(orbsum[bits_aut[orbit_aut[e]]], orbsum[bits_aut[orbit_aut[e]]], Raut[orbit_aut[e]]);
	}
	for (b = 0; b < nb; b++) {
	  if (fmpq_mat_is_zero(orbsum[b]))
	    continue;
	  fmpq_mat_mul(prod, orbsum[b], Rs);
	  fmpq_mat_add(acc[i][b ^ bits_s], acc[i][b ^ bits_s], prod);
	}
	isometry_clear(s_nbr);
	isometry_clear(hash_isom);
	square_matrix_clear(nbr);
      }
      nbr_process_advance(nbr_man);
    }
    nbr_process_clear(nbr_man);
  }

  if (has_auts) {
    for (k = 0; k < num_auts; k++)
      fmpq_mat_clear(Raut[k]);
    free(Raut);
    free(bits_aut);
    free(orbit);
    free(orbit_aut);
    free(auts);
  }
  for (b = 0; b < nb; b++)
    fmpq_mat_clear(orbsum[b]);
  free(orbsum);
  fmpq_mat_clear(Rs);
  fmpq_mat_clear(prod);
  return;
}

// ---- The eigenvalues of T_p and T_{p^2} by evaluation at a point ----
//
// For a right eigenvector u of the Hecke operators (T u = lambda u), with the vector U_i of W of its block
// at the class i, sum_i M_ji U_i = lambda U_j, where M_ji = sum of chi(bits) p^a rho(s) over the
// neighbors of L_j of class i, as in hecke_matrices_weight_impl. rho(s) takes the polynomial F_w of w
// to sign(s) F_w(h x, h y), h = s^(-1), so the evaluation at the point (x0, y0) gives
//   lambda F_{U_j}(x0, y0) = p^a sum over the neighbors N of L_j of chi(bits) sign(s) F_{U_i}(h x0, h y0),
// a sum over the neighbors of one class j of the values of the monomials at (h x0, h y0), which are kept
// by class and spinor bits (weight_moments), and then paired with the coefficients of the forms.

static const Z64 weight_x0[QF_RANK] = {1, 3, -2, 5, 7};
static const Z64 weight_y0[QF_RANK] = {2, -1, 4, -3, 1};

// the sums of the values of the monomials, with a numerator vector for each denominator
typedef struct {
  slong n, cap, len;
  fmpz* dens;
  fmpz** nums;
} weight_mom_struct;

static void weight_mom_init(weight_mom_struct* A, slong len)
{
  A->n = 0;
  A->cap = 0;
  A->len = len;
  A->dens = NULL;
  A->nums = NULL;
}

static void weight_mom_clear(weight_mom_struct* A)
{
  slong e;
  for (e = 0; e < A->n; e++) {
    fmpz_clear(A->dens + e);
    _fmpz_vec_clear(A->nums[e], A->len);
  }
  free(A->dens);
  free(A->nums);
}

// adds vec / d
static void weight_mom_add(weight_mom_struct* A, const fmpz_t d, const fmpz* vec)
{
  slong e;
  for (e = 0; (e < A->n) && !fmpz_equal(A->dens + e, d); e++);
  if (e == A->n) {
    if (A->n == A->cap) {
      slong t;
      A->cap = (A->cap == 0) ? 4 : 2 * A->cap;
      A->dens = (fmpz*)realloc(A->dens, A->cap * sizeof(fmpz));
      A->nums = (fmpz**)realloc(A->nums, A->cap * sizeof(fmpz*));
      for (t = A->n; t < A->cap; t++)
	fmpz_init(A->dens + t);
    }
    fmpz_set(A->dens + e, d);
    A->nums[e] = _fmpz_vec_init(A->len);
    A->n++;
  }
  _fmpz_vec_add(A->nums[e], A->nums[e], vec, A->len);
}

// vec[kx * ny + ky] = X^ex[kx] Y^ey[ky], the bimonomials of bidegree (a, b) at the integral point (X, Y)
static void weight_mom_values(fmpz* vec, const weight_rep_t W, const fmpz* X, const fmpz* Y, fmpz* xm, fmpz* ym)
{
  slong kx, ky;
  int r;
  fmpz_t pw;

  fmpz_init(pw);
  for (kx = 0; kx < W->nx; kx++) {
    fmpz_one(xm + kx);
    for (r = 0; r < QF_RANK; r++)
      if (W->ex[QF_RANK * kx + r] > 0) {
	fmpz_pow_ui(pw, X + r, W->ex[QF_RANK * kx + r]);
	fmpz_mul(xm + kx, xm + kx, pw);
      }
  }
  for (ky = 0; ky < W->ny; ky++) {
    fmpz_one(ym + ky);
    for (r = 0; r < QF_RANK; r++)
      if (W->ey[QF_RANK * ky + r] > 0) {
	fmpz_pow_ui(pw, Y + r, W->ey[QF_RANK * ky + r]);
	fmpz_mul(ym + ky, ym + ky, pw);
      }
  }
  for (kx = 0; kx < W->nx; kx++)
    for (ky = 0; ky < W->ny; ky++)
      fmpz_mul(vec + kx * W->ny + ky, xm + kx, ym + ky);
  fmpz_clear(pw);
}

// sign(s) = det(s)^(a + b)
static int weight_isom_sign(const weight_rep_t W, const isometry_t s)
{
  fmpz_mat_t S;
  fmpz_t det;
  int sign;

  if (((W->a + W->b) % 2) == 0)
    return 1;
  fmpz_mat_init_set_square_matrix(S, s->s);
  fmpz_init(det);
  fmpz_mat_det(det, S);
  sign = ((fmpz_sgn(det) < 0) != (s->denom < 0)) ? -1 : 1;
  fmpz_clear(det);
  fmpz_mat_clear(S);
  return sign;
}

// the contribution of the isometry s to the neighbor class i, with the point (t^(-1) x0, t^(-1) y0) given as
// the integral vectors TX, TY over dt: adds the values at (h TX / dt, h TY / dt), h = s^(-1), to mom[i][bits]
static void weight_mom_point(weight_mom_struct** mom, slong i, W64 bits, int sign, const weight_rep_t W,
			     const isometry_t s, const fmpz* TX, const fmpz* TY, const fmpz_t dt, fmpz* vec,
			     fmpz* xm, fmpz* ym)
{
  fmpz X[QF_RANK], Y[QF_RANK];
  fmpz_t D, d;
  int r, k;

  for (r = 0; r < QF_RANK; r++) {
    fmpz_init(X + r);
    fmpz_init(Y + r);
    for (k = 0; k < QF_RANK; k++) {
      fmpz_addmul_si(X + r, TX + k, s->s_inv[r][k]);
      fmpz_addmul_si(Y + r, TY + k, s->s_inv[r][k]);
    }
  }
  fmpz_init(D);
  fmpz_init(d);
  fmpz_mul_si(D, dt, s->inv_denom);
  fmpz_pow_ui(d, D, (ulong)(W->a + W->b));
  if (sign < 0)
    fmpz_neg(d, d);
  weight_mom_values(vec, W, X, Y, xm, ym);
  weight_mom_add(&(mom[i][bits]), d, vec);
  fmpz_clear(D);
  fmpz_clear(d);
  for (r = 0; r < QF_RANK; r++) {
    fmpz_clear(X + r);
    fmpz_clear(Y + r);
  }
}

// TX = t->s_inv X0 and TY = t->s_inv Y0, over dt = t->inv_denom: the point t^(-1) (X0, Y0)
static void weight_inv_point(fmpz* TX, fmpz* TY, fmpz_t dt, const isometry_t t, const fmpz* X0, const fmpz* Y0)
{
  int r, k;
  for (r = 0; r < QF_RANK; r++) {
    fmpz_zero(TX + r);
    fmpz_zero(TY + r);
    for (k = 0; k < QF_RANK; k++) {
      fmpz_addmul_si(TX + r, X0 + k, t->s_inv[r][k]);
      fmpz_addmul_si(TY + r, Y0 + k, t->s_inv[r][k]);
    }
  }
  fmpz_set_si(dt, t->inv_denom);
}

// mom[t][i][b] for the p^k-neighbors of L_j, at the npts integral points (PX + QF_RANK t, PY + QF_RANK t): for
// k = 1 by the orbits of the automorphisms, as weight_sums_orbits, and for k = 2 every neighbor of nbr_data
static void weight_moments(weight_mom_struct*** mom, const genus_t genus, const weight_rep_t W, slong j, int p, int k,
			   slong npts, const fmpz* PX, const fmpz* PY)
{
  slong len = W->nx * W->ny, num_isom, kk, e, num_auts, orb_size, i, pt;
  fmpz *vec, *xm, *ym, *TX, *TY, *dts;
  fmpz_t one;
  double theta_time, isom_time;
  isometry_t s_nbr, hash_isom, s_inv, t;
  square_matrix_t nbr, aut_tr;

  vec = _fmpz_vec_init(len);
  xm = _fmpz_vec_init(W->nx);
  ym = _fmpz_vec_init(W->ny);
  fmpz_init_set_ui(one, 1);
  theta_time = isom_time = 0;
  num_isom = 0;

  if (k == 2) {
    nbr_data_t nd;
    fmpz_mat_t nbr_fmpz, nbr_isom;
    fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
    fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
    nbr_data_init(nd, genus->genus_reps->keys[j], p, 2);
    while (!nbr_data_has_ended(nd)) {
      nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nd);
      square_matrix_set_fmpz_mat(nbr, nbr_fmpz);
      isometry_init_set_fmpz_mat(s_nbr, nbr_isom, p);
      i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, &theta_time, &isom_time, &num_isom);
      isometry_muleq_left(s_nbr, genus->isoms[j]);
      isometry_muleq_right(s_nbr, hash_isom);
      isometry_inv(s_inv, genus->isoms[i]);
      isometry_muleq_right(s_nbr, s_inv);
      {
	W64 bits_s = spinor_norm_isom(genus->spinor, s_nbr);
	int sign_s = weight_isom_sign(W, s_nbr);
	for (pt = 0; pt < npts; pt++)
	  weight_mom_point(mom[pt], i, bits_s, sign_s, W, s_nbr, PX + QF_RANK * pt, PY + QF_RANK * pt, one, vec, xm, ym);
      }
      nbr_data_get_next_neighbor(nd);
    }
    nbr_data_clear(nd);
    fmpz_mat_clear(nbr_fmpz);
    fmpz_mat_clear(nbr_isom);
  }
  else {
    neighbor_manager_t nbr_man;
    W64* bits_aut = NULL;
    int* sign_aut = NULL;
    vector_t* orbit = NULL;
    slong* orbit_aut = NULL;
    square_matrix_t* auts = NULL;
    vector_t g_vec;
    bool has_auts = false, found;
    int line, vec_cmp;
    TX = TY = dts = NULL;
    num_auts = 0;
    for (line = 0; line < p; line++) {
      nbr_process_init(nbr_man, genus->genus_reps->keys[j], p, line, genus->isoms[j]);
      if (!has_auts) {
	// t = isoms[j] g isoms[j]^(-1) of every automorphism g, and t^(-1) of each point
	num_auts = nbr_man->num_auts;
	bits_aut = (W64*)malloc(num_auts * sizeof(W64));
	sign_aut = (int*)malloc(num_auts * sizeof(int));
	TX = _fmpz_vec_init(QF_RANK * num_auts * npts);
	TY = _fmpz_vec_init(QF_RANK * num_auts * npts);
	dts = _fmpz_vec_init(num_auts * npts);
	isometry_inv(s_inv, genus->isoms[j]);
	for (kk = 0; kk < num_auts; kk++) {
	  square_matrix_transpose(aut_tr, nbr_man->auts[kk]);
	  isometry_init_set_square_matrix(t, aut_tr, 1);
	  isometry_muleq_left(t, genus->isoms[j]);
	  isometry_muleq_right(t, s_inv);
	  bits_aut[kk] = spinor_norm_isom(genus->spinor, t);
	  sign_aut[kk] = weight_isom_sign(W, t);
	  for (pt = 0; pt < npts; pt++)
	    weight_inv_point(TX + QF_RANK * (pt * num_auts + kk), TY + QF_RANK * (pt * num_auts + kk),
			     dts + pt * num_auts + kk, t, PX + QF_RANK * pt, PY + QF_RANK * pt);
	  isometry_clear(t);
	}
	orbit = (vector_t*)malloc((num_auts + 1) * sizeof(vector_t));
	orbit_aut = (slong*)malloc((num_auts + 1) * sizeof(slong));
	auts = (square_matrix_t*)malloc(num_auts * sizeof(square_matrix_t));
	for (kk = 0; kk < num_auts; kk++)
	  square_matrix_set(auts[kk], nbr_man->auts[kk]);
	has_auts = true;
      }
      while (!nbr_process_has_ended(nbr_man)) {
	vector_mod_p(nbr_man->iso_vec, nbr_man->p);
	normalize_mod_p(nbr_man->iso_vec, nbr_man->p);
	vector_set(orbit[0], nbr_man->iso_vec);
	orbit_aut[0] = -1;
	orb_size = 1;
	vec_cmp = 1;
	for (kk = 0; kk < num_auts; kk++) {
	  square_matrix_mul_vec_left(g_vec, nbr_man->iso_vec, auts[kk]);
	  vector_mod_p(g_vec, nbr_man->p);
	  normalize_mod_p(g_vec, nbr_man->p);
	  vec_cmp = vector_cmp(g_vec, nbr_man->iso_vec);
	  if (vec_cmp < 0)
	    break;
	  if (vec_cmp > 0) {
	    found = false;
	    for (e = 0; (e < orb_size) && !found; e++)
	      found = (vector_cmp(g_vec, orbit[e]) == 0);
	    if (!found) {
	      vector_set(orbit[orb_size], g_vec);
	      orbit_aut[orb_size] = kk;
	      orb_size++;
	    }
	  }
	}
	if ((vec_cmp >= 0) && nbr_process_build_nb_and_isom(nbr, s_nbr, nbr_man)) {
	  W64 bits_s;
	  int sign_s;
	  i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, &theta_time, &isom_time, &num_isom);
	  isometry_muleq_left(s_nbr, genus->isoms[j]);
	  isometry_muleq_right(s_nbr, hash_isom);
	  isometry_inv(t, genus->isoms[i]);
	  isometry_muleq_right(s_nbr, t);
	  bits_s = spinor_norm_isom(genus->spinor, s_nbr);
	  sign_s = weight_isom_sign(W, s_nbr);
	  // the element t s of the orbit is at the point s^(-1) t^(-1) (X0, Y0)
	  for (pt = 0; pt < npts; pt++)
	    for (e = 0; e < orb_size; e++) {
	      if (orbit_aut[e] < 0)
		weight_mom_point(mom[pt], i, bits_s, sign_s, W, s_nbr, PX + QF_RANK * pt, PY + QF_RANK * pt, one, vec,
				 xm, ym);
	      else {
		kk = pt * num_auts + orbit_aut[e];
		weight_mom_point(mom[pt], i, bits_s ^ bits_aut[orbit_aut[e]], sign_s * sign_aut[orbit_aut[e]], W, s_nbr,
				 TX + QF_RANK * kk, TY + QF_RANK * kk, dts + kk, vec, xm, ym);
	      }
	    }
	  isometry_clear(s_nbr);
	  isometry_clear(hash_isom);
	  square_matrix_clear(nbr);
	}
	nbr_process_advance(nbr_man);
      }
      nbr_process_clear(nbr_man);
    }
    if (has_auts) {
      free(bits_aut);
      free(sign_aut);
      _fmpz_vec_clear(TX, QF_RANK * num_auts * npts);
      _fmpz_vec_clear(TY, QF_RANK * num_auts * npts);
      _fmpz_vec_clear(dts, num_auts * npts);
      free(orbit);
      free(orbit_aut);
      free(auts);
    }
  }

  _fmpz_vec_clear(vec, len);
  _fmpz_vec_clear(xm, W->nx);
  _fmpz_vec_clear(ym, W->ny);
  fmpz_clear(one);
}

// The polynomials of the basis vectors of H[c][i], as integral rows Bnum[m] over Bden[m].
typedef struct
{
  fmpz_mat_t num;
  fmpz* den;
  slong rows;
} weight_bpoly_struct;

static void weight_bpoly_init(weight_bpoly_struct* P, const fmpq_mat_t H, const weight_rep_t W)
{
  fmpq_mat_t Hs;
  fmpz_mat_t Hz, Bz;
  fmpz_t g;
  slong m, k, kx, ky, len = W->nx * W->ny;

  // B[m][kx ny + ky] = sum_k H[m][k] basis_z[kx][k ny + ky] / basis_den[k], the polynomial of the m-th basis
  // vector of H. So B = H D^(-1) Bz, with D = diag(basis_den) and Bz[k][kx ny + ky] = basis_z[kx][k ny + ky]:
  // the rows of H D^(-1) over their denominators, times Bz, and each row then in lowest terms, with the
  // smallest common denominator, as fmpq_mat_get_fmpz_mat_rowwise of B would give
  P->rows = fmpq_mat_nrows(H);
  fmpz_mat_init(P->num, (P->rows > 0) ? P->rows : 1, len);
  P->den = _fmpz_vec_init((P->rows > 0) ? P->rows : 1);
  if (P->rows == 0)
    return;
  fmpq_mat_init(Hs, P->rows, W->dim);
  for (m = 0; m < P->rows; m++)
    for (k = 0; k < W->dim; k++)
      fmpq_div_fmpz(fmpq_mat_entry(Hs, m, k), fmpq_mat_entry(H, m, k), W->basis_den + k);
  fmpz_mat_init(Hz, P->rows, W->dim);
  fmpq_mat_get_fmpz_mat_rowwise(Hz, P->den, Hs);
  fmpq_mat_clear(Hs);
  fmpz_mat_init(Bz, W->dim, len);
  for (k = 0; k < W->dim; k++)
    for (kx = 0; kx < W->nx; kx++)
      for (ky = 0; ky < W->ny; ky++)
	fmpz_set(fmpz_mat_entry(Bz, k, kx * W->ny + ky), fmpz_mat_entry(W->basis_z, kx, k * W->ny + ky));
  fmpz_mat_mul(P->num, Hz, Bz);
  fmpz_mat_clear(Hz);
  fmpz_mat_clear(Bz);
  fmpz_init(g);
  for (m = 0; m < P->rows; m++) {
    _fmpz_vec_content(g, fmpz_mat_entry(P->num, m, 0), len);
    fmpz_gcd(g, g, P->den + m);
    if (fmpz_is_zero(g) || fmpz_is_one(g))
      continue;
    _fmpz_vec_scalar_divexact_fmpz(fmpz_mat_entry(P->num, m, 0), fmpz_mat_entry(P->num, m, 0), len, g);
    fmpz_divexact(P->den + m, P->den + m, g);
  }
  fmpz_clear(g);
}

static void weight_bpoly_clear(weight_bpoly_struct* P)
{
  fmpz_mat_clear(P->num);
  _fmpz_vec_clear(P->den, (P->rows > 0) ? P->rows : 1);
}

// r[offsets[c][i] + m] = sum_b chi_c(b) <B_i[m], mom[i][b]>, the polynomial B_i[m] of the m-th basis vector
// of H[c][i] paired with the sums of the neighbors of class i: then lambda = pa <u, r> / phi for each form of
// the conductor index c whose pivot is the class of mom
static void weight_eval_pairing(fmpq* r, const weight_bpoly_struct* P, const weight_space_t S, slong c,
				weight_mom_struct** mom, slong h, slong nb, slong len)
{
  slong i, b, e, m, mono;
  fmpz_t dot;
  fmpq_t q;

  fmpz_init(dot);
  fmpq_init(q);
  for (i = 0; i < h; i++)
    for (m = 0; m < P[i].rows; m++) {
      fmpq* ri = r + S->offsets[c][i] + m;
      fmpq_zero(ri);
      for (b = 0; b < nb; b++)
	for (e = 0; e < mom[i][b].n; e++) {
	  fmpz_zero(dot);
	  for (mono = 0; mono < len; mono++)
	    fmpz_addmul(dot, fmpz_mat_entry(P[i].num, m, mono), mom[i][b].nums[e] + mono);
	  if (fmpz_is_zero(dot))
	    continue;
	  fmpz_mul(fmpq_denref(q), mom[i][b].dens + e, P[i].den + m);
	  fmpz_set(fmpq_numref(q), dot);
	  fmpq_canonicalise(q);
	  if (weight_char((W64)b, c) > 0)
	    fmpq_add(ri, ri, q);
	  else
	    fmpq_sub(ri, ri, q);
	}
    }
  fmpz_clear(dot);
  fmpq_clear(q);
}

// ---- The Hecke matrices by evaluation ----
//
// The block (j, i) of T at the conductor index c takes the k-th basis vector w of H[c][i] to the vector U of
// H[c][j] with F_U(x, y) = p^a sum over the neighbors N of L_j of class i of chi_c(bits) sign(s) F_w(h x, h y),
// h = s^(-1) (see weight_eval_pairing). The coordinates u of U in H[c][j] are found from the values of F_U at
// dim H[c][j] points at which the polynomials of the basis of H[c][j] are independent, E u = (F_U(z_r))_r with
// E[r][m] the value at z_r of the m-th of them. The value at one more point checks that U is in H[c][j], as
// weight_add_block does. So no matrix of rho(s) on W is needed.

// The point of index t: (x0, y0) for t = 0, and otherwise with entries in [-w, w] from a fixed sequence,
// w = 7 + t / 64, so that the candidates are not confined to a grid on which a polynomial of the weight vanishes.
// The points of the check, seed 1 and 2 of another sequence, have entries in [-60, 60].
static void weight_eval_point_seq(fmpz* X, fmpz* Y, ulong seed, slong w)
{
  ulong st;
  int r;

  st = seed * UWORD(0x9E3779B97F4A7C15);
  for (r = 0; r < 2 * QF_RANK; r++) {
    st = st * UWORD(6364136223846793005) + UWORD(1442695040888963407);
    if (r < QF_RANK)
      fmpz_set_si(X + r, (slong)((st >> 33) % (ulong)(2 * w + 1)) - w);
    else
      fmpz_set_si(Y + r - QF_RANK, (slong)((st >> 33) % (ulong)(2 * w + 1)) - w);
  }
}

static void weight_eval_candidate(fmpz* X, fmpz* Y, slong t)
{
  int r;

  if (t == 0) {
    for (r = 0; r < QF_RANK; r++) {
      fmpz_set_si(X + r, weight_x0[r]);
      fmpz_set_si(Y + r, weight_y0[r]);
    }
    return;
  }
  weight_eval_point_seq(X, Y, (ulong)t, 7 + t / 64);
}

static void weight_eval_check_point(fmpz* X, fmpz* Y, int q)
{
  weight_eval_point_seq(X, Y, UWORD(0xC0FFEE) + (ulong)q, 60);
}

// v[m] = the value of the m-th polynomial of P at the point (X, Y)
static void weight_bpoly_values(fmpq* v, const weight_bpoly_struct* P, const weight_rep_t W, const fmpz* X,
				const fmpz* Y, fmpz* vec, fmpz* xm, fmpz* ym)
{
  slong m, mono, len = W->nx * W->ny;
  fmpz_t dot;

  fmpz_init(dot);
  weight_mom_values(vec, W, X, Y, xm, ym);
  for (m = 0; m < P->rows; m++) {
    fmpz_zero(dot);
    for (mono = 0; mono < len; mono++)
      fmpz_addmul(dot, fmpz_mat_entry(P->num, m, mono), vec + mono);
    fmpq_set_fmpz_frac(v + m, dot, P->den + m);
  }
  fmpz_clear(dot);
}

#define WEIGHT_EVAL_CHECKS 2
#ifndef WEIGHT_EVAL_FACTOR
#define WEIGHT_EVAL_FACTOR 4.0
#endif

// The points of the class j: PX, PY hold npts points, the last WEIGHT_EVAL_CHECKS of them only for the check.
// For each conductor index c, sel[c][r], r < dim H[c][j], is the index of the r-th point of E[c], Einv[c] the
// inverse of E[c], and chk[c] (WEIGHT_EVAL_CHECKS x dim H[c][j]) the values of the polynomials of H[c][j] at the
// points of the check.
typedef struct
{
  slong npts;
  fmpz *PX, *PY;
  slong** sel;
  fmpq_mat_struct* Einv;
  fmpq_mat_struct* chk;
  slong nc, cap;
} weight_eval_points_struct;

static void weight_eval_points_init(weight_eval_points_struct* E, const weight_bpoly_struct* P, slong nc,
				    const weight_rep_t W)
{
  slong c, t, cap, n, done, m, r, q, len = W->nx * W->ny;
  slong *cnt, **piv;
  fmpq_mat_struct* rows;
  ulong ***red, *row, prime, f;
  nmod_t mod;
  fmpz *X, *Y, *vec, *xm, *ym;
  fmpq* v;
  bool used, ok;

  E->nc = nc;
  E->sel = (slong**)malloc(nc * sizeof(slong*));
  E->Einv = (fmpq_mat_struct*)malloc(nc * sizeof(fmpq_mat_struct));
  E->chk = (fmpq_mat_struct*)malloc(nc * sizeof(fmpq_mat_struct));
  cnt = (slong*)calloc(nc, sizeof(slong));
  rows = (fmpq_mat_struct*)malloc(nc * sizeof(fmpq_mat_struct));
  piv = (slong**)malloc(nc * sizeof(slong*));
  red = (ulong***)malloc(nc * sizeof(ulong**));
  // each point kept raises the rank for some c, so there are at most sum_c dim H[c][j] of them, and the checks
  n = 0;
  cap = WEIGHT_EVAL_CHECKS;
  for (c = 0; c < nc; c++) {
    slong d = (P[c].rows > 0) ? P[c].rows : 1;
    E->sel[c] = (slong*)malloc(d * sizeof(slong));
    fmpq_mat_init(rows + c, d, d);
    piv[c] = (slong*)malloc(d * sizeof(slong));
    red[c] = (ulong**)malloc(d * sizeof(ulong*));
    for (r = 0; r < d; r++)
      red[c][r] = (ulong*)calloc(d, sizeof(ulong));
    n = (P[c].rows > n) ? P[c].rows : n;
    cap += P[c].rows;
  }
  E->cap = cap;
  E->PX = _fmpz_vec_init(QF_RANK * cap);
  E->PY = _fmpz_vec_init(QF_RANK * cap);
  X = _fmpz_vec_init(QF_RANK);
  Y = _fmpz_vec_init(QF_RANK);
  vec = _fmpz_vec_init(len);
  xm = _fmpz_vec_init(W->nx);
  ym = _fmpz_vec_init(W->ny);
  v = _fmpq_vec_init(n + 1);
  row = (ulong*)malloc((n + 1) * sizeof(ulong));
  prime = n_nextprime(UWORD(1) << 62, 1);
  nmod_init(&mod, prime);

  // a candidate is kept for c when its row of values of the polynomials of H[c][j] is independent of the rows
  // kept for c modulo a prime, which is reduced against them (the rows in echelon form red[c], with the pivots
  // piv[c]); then it is independent over Q too, so E[c] is invertible
  E->npts = 0;
  done = 0;
  for (c = 0; c < nc; c++)
    done += (P[c].rows == 0);
  for (t = 0; done < nc; t++) {
    // checked in every build: the polynomials of a basis are independent, so some points will do
    if (t > 100000) {
      fprintf(stderr, "Error: no points at which the forms of a weight are independent were found\n");
      abort();
    }
    weight_eval_candidate(X, Y, t);
    used = false;
    for (c = 0; c < nc; c++) {
      if (cnt[c] == P[c].rows)
	continue;
      weight_bpoly_values(v, P + c, W, X, Y, vec, xm, ym);
      ok = true;
      for (m = 0; (m < P[c].rows) && ok; m++) {
	f = fmpz_fdiv_ui(fmpq_denref(v + m), prime);
	ok = (f != 0);
	if (ok)
	  row[m] = nmod_mul(fmpz_fdiv_ui(fmpq_numref(v + m), prime), n_invmod(f, prime), mod);
      }
      if (!ok)
	continue;
      for (r = 0; r < cnt[c]; r++) {
	f = row[piv[c][r]];
	if (f != 0)
	  for (m = 0; m < P[c].rows; m++)
	    row[m] = nmod_sub(row[m], nmod_mul(f, red[c][r][m], mod), mod);
      }
      for (m = 0; (m < P[c].rows) && (row[m] == 0); m++);
      if (m == P[c].rows)
	continue;
      piv[c][cnt[c]] = m;
      f = n_invmod(row[m], prime);
      for (q = 0; q < P[c].rows; q++)
	red[c][cnt[c]][q] = nmod_mul(row[q], f, mod);
      for (q = 0; q < P[c].rows; q++)
	fmpq_set(fmpq_mat_entry(rows + c, cnt[c], q), v + q);
      E->sel[c][cnt[c]++] = E->npts;
      used = true;
      if (cnt[c] == P[c].rows)
	done++;
    }
    if (used) {
      _fmpz_vec_set(E->PX + QF_RANK * E->npts, X, QF_RANK);
      _fmpz_vec_set(E->PY + QF_RANK * E->npts, Y, QF_RANK);
      E->npts++;
    }
  }
  for (c = 0; c < nc; c++) {
    fmpq_mat_init(E->Einv + c, P[c].rows, P[c].rows);
    if (P[c].rows > 0)
      fmpq_mat_inv(E->Einv + c, rows + c);
    fmpq_mat_init(E->chk + c, WEIGHT_EVAL_CHECKS, P[c].rows);
  }
  // the points of the check
  for (q = 0; q < WEIGHT_EVAL_CHECKS; q++) {
    weight_eval_check_point(X, Y, q);
    _fmpz_vec_set(E->PX + QF_RANK * E->npts, X, QF_RANK);
    _fmpz_vec_set(E->PY + QF_RANK * E->npts, Y, QF_RANK);
    E->npts++;
    for (c = 0; c < nc; c++) {
      if (P[c].rows == 0)
	continue;
      weight_bpoly_values(v, P + c, W, X, Y, vec, xm, ym);
      for (m = 0; m < P[c].rows; m++)
	fmpq_set(fmpq_mat_entry(E->chk + c, q, m), v + m);
    }
  }

  for (c = 0; c < nc; c++) {
    slong d = (P[c].rows > 0) ? P[c].rows : 1;
    fmpq_mat_clear(rows + c);
    for (r = 0; r < d; r++)
      free(red[c][r]);
    free(red[c]);
    free(piv[c]);
  }
  free(red);
  free(piv);
  free(row);
  free(rows);
  free(cnt);
  _fmpz_vec_clear(X, QF_RANK);
  _fmpz_vec_clear(Y, QF_RANK);
  _fmpz_vec_clear(vec, len);
  _fmpz_vec_clear(xm, W->nx);
  _fmpz_vec_clear(ym, W->ny);
  _fmpq_vec_clear(v, n + 1);
}

static void weight_eval_points_clear(weight_eval_points_struct* E)
{
  slong c;
  for (c = 0; c < E->nc; c++) {
    free(E->sel[c]);
    fmpq_mat_clear(E->Einv + c);
    fmpq_mat_clear(E->chk + c);
  }
  free(E->sel);
  free(E->Einv);
  free(E->chk);
  _fmpz_vec_clear(E->PX, QF_RANK * E->cap);
  _fmpz_vec_clear(E->PY, QF_RANK * E->cap);
}

// Rv[pt][c][offsets[c][i] + m] = sum_b chi_c(b) <F_m, mom[pt][i][b]>, the values at the point pt of the images
// (without p^a) of the basis vectors of conductor c, with F_m the polynomial of the m-th basis vector of H[c][i].
// F_m = sum_k H[c][i][m][k] F_{e_k}, F_{e_k} the k-th basis polynomial of W, basis_z / basis_den, so the moments
// are paired with these first, all of them in one product with BzT[mono][k] = basis_z[kx][k ny + ky]; this needs
// no polynomial of the forms of the other classes.
static void weight_eval_values(fmpq*** Rv, weight_mom_struct*** mom, slong npts, const fmpz_mat_t BzT,
			       const weight_rep_t W, const weight_space_t S, slong h, slong nb, slong nc)
{
  fmpz_mat_t Mu, Q;
  fmpq_t t;
  fmpz_t den;
  fmpq **sv;
  slong n, row, pt, i, b, e, c, k, m, len = W->nx * W->ny;

  n = 0;
  for (pt = 0; pt < npts; pt++)
    for (i = 0; i < h; i++)
      for (b = 0; b < nb; b++)
	n += mom[pt][i][b].n;
  fmpz_mat_init(Mu, (n > 0) ? n : 1, len);
  row = 0;
  for (pt = 0; pt < npts; pt++)
    for (i = 0; i < h; i++)
      for (b = 0; b < nb; b++)
	for (e = 0; e < mom[pt][i][b].n; e++, row++)
	  _fmpz_vec_set(fmpz_mat_entry(Mu, row, 0), mom[pt][i][b].nums[e], len);
  fmpz_mat_init(Q, (n > 0) ? n : 1, W->dim);
  if (n > 0)
    fmpz_mat_mul(Q, Mu, BzT);
  fmpz_mat_clear(Mu);

  fmpq_init(t);
  fmpz_init(den);
  sv = (fmpq**)malloc(nc * sizeof(fmpq*));
  for (c = 0; c < nc; c++)
    sv[c] = _fmpq_vec_init(W->dim);
  for (pt = 0; pt < npts; pt++)
    for (c = 0; c < nc; c++)
      Rv[pt][c] = _fmpq_vec_init((S->dims[c] > 0) ? S->dims[c] : 1);
  row = 0;
  for (pt = 0; pt < npts; pt++)
    for (i = 0; i < h; i++) {
      // sv[c][k] = sum_b chi_c(b) <F_{e_k}, mom[pt][i][b]>
      for (c = 0; c < nc; c++)
	for (k = 0; k < W->dim; k++)
	  fmpq_zero(sv[c] + k);
      for (b = 0; b < nb; b++)
	for (e = 0; e < mom[pt][i][b].n; e++, row++)
	  for (k = 0; k < W->dim; k++) {
	    if (fmpz_is_zero(fmpz_mat_entry(Q, row, k)))
	      continue;
	    fmpz_mul(den, W->basis_den + k, mom[pt][i][b].dens + e);
	    fmpq_set_fmpz_frac(t, fmpz_mat_entry(Q, row, k), den);
	    for (c = 0; c < nc; c++) {
	      if (fmpq_mat_nrows(S->H[c] + i) == 0)
		continue;
	      if (weight_char((W64)b, c) > 0)
		fmpq_add(sv[c] + k, sv[c] + k, t);
	      else
		fmpq_sub(sv[c] + k, sv[c] + k, t);
	    }
	  }
      for (c = 0; c < nc; c++)
	for (m = 0; m < fmpq_mat_nrows(S->H[c] + i); m++) {
	  fmpq* rr = Rv[pt][c] + S->offsets[c][i] + m;
	  for (k = 0; k < W->dim; k++)
	    if (!fmpq_is_zero(fmpq_mat_entry(S->H[c] + i, m, k)))
	      fmpq_addmul(rr, fmpq_mat_entry(S->H[c] + i, m, k), sv[c] + k);
	}
    }
  for (c = 0; c < nc; c++)
    _fmpq_vec_clear(sv[c], W->dim);
  free(sv);
  fmpq_clear(t);
  fmpz_clear(den);
  fmpz_mat_clear(Q);
}

static void hecke_weight_rho_classes(fmpq_mat_struct* T, const genus_t genus, const weight_rep_t W,
				     const weight_space_t S, int p, int k, const bool* classes, bool use_orbits);

// The matrices of T_{p^k} on the spaces of each conductor, as hecke_matrices_weight_impl. The rows of a class j
// are found by evaluation when that is cheaper, when the npts points times the nmons bimonomials of each moved
// point are at most WEIGHT_EVAL_FACTOR dim W^2, dim W^2 being the size of a matrix of rho(s) (4 was the fastest
// of 1/4, 1, 4, 16 and infinity, on the weights (1, 0) to (6, 3)); otherwise, when L_j has few automorphisms and
// so H[c][j] is almost as large as W, from the matrices of rho(s), by the orbits of the automorphisms for T_p.
static fmpq_mat_struct* hecke_matrices_weight_eval(const genus_t genus, const weight_rep_t W, const weight_space_t S,
						   int p, int k)
{
  fmpq_mat_struct* T;
  weight_bpoly_struct* Pj;
  weight_eval_points_struct E;
  weight_mom_struct*** mom;
  fmpq*** Rv;
  fmpz_mat_t BzT;
  fmpz_t pa;
  slong c, i, j, h, nb, b, nc, len, pt, m, r, col, npts, k2, kx, ky;
  bool needed, any_rho, any_eval, *by_rho;

  h = genus->genus_reps->num_stored;
  nb = (slong)1 << genus->spinor->num_primes;
  nc = genus->num_conductors;
  len = W->nx * W->ny;
  T = (fmpq_mat_struct*)malloc(nc * sizeof(fmpq_mat_struct));
  for (c = 0; c < nc; c++)
    fmpq_mat_init(T + c, S->dims[c], S->dims[c]);
  fmpz_init(pa);
  fmpz_set_ui(pa, p);
  fmpz_pow_ui(pa, pa, (k == 1) ? W->a : 2 * W->a);

  // the classes whose rows are found from the matrices of rho(s)
  by_rho = (bool*)calloc((h > 0) ? h : 1, sizeof(bool));
  any_rho = any_eval = false;
  for (j = 0; j < h; j++) {
    needed = false;
    npts = WEIGHT_EVAL_CHECKS;
    for (c = 0; c < nc; c++) {
      needed = needed || (fmpq_mat_nrows(S->H[c] + j) > 0);
      npts += fmpq_mat_nrows(S->H[c] + j);
    }
    if (!needed)
      continue;
    if ((double)npts * (double)len > WEIGHT_EVAL_FACTOR * (double)W->dim * (double)W->dim)
      by_rho[j] = any_rho = true;
    else
      any_eval = true;
  }

  // BzT[mono][k] = basis_z[kx][k ny + ky], when some class is found by evaluation
  fmpz_mat_init(BzT, any_eval ? len : 1, W->dim);
  for (k2 = 0; (k2 < W->dim) && any_eval; k2++)
    for (kx = 0; kx < W->nx; kx++)
      for (ky = 0; ky < W->ny; ky++)
	fmpz_set(fmpz_mat_entry(BzT, kx * W->ny + ky, k2), fmpz_mat_entry(W->basis_z, kx, k2 * W->ny + ky));
  Pj = (weight_bpoly_struct*)malloc(nc * sizeof(weight_bpoly_struct));

  for (j = 0; j < h; j++) {
    needed = false;
    for (c = 0; c < nc; c++)
      needed = needed || (fmpq_mat_nrows(S->H[c] + j) > 0);
    if (!needed || by_rho[j])
      continue;
    // the polynomials of the basis of H[c][j], for the points
    for (c = 0; c < nc; c++)
      weight_bpoly_init(Pj + c, S->H[c] + j, W);
    weight_eval_points_init(&E, Pj, nc, W);

    mom = (weight_mom_struct***)malloc(E.npts * sizeof(weight_mom_struct**));
    for (pt = 0; pt < E.npts; pt++) {
      mom[pt] = (weight_mom_struct**)malloc(h * sizeof(weight_mom_struct*));
      for (i = 0; i < h; i++) {
	mom[pt][i] = (weight_mom_struct*)malloc(nb * sizeof(weight_mom_struct));
	for (b = 0; b < nb; b++)
	  weight_mom_init(&(mom[pt][i][b]), len);
      }
    }
    weight_moments(mom, genus, W, j, p, k, E.npts, E.PX, E.PY);
    Rv = (fmpq***)malloc(E.npts * sizeof(fmpq**));
    for (pt = 0; pt < E.npts; pt++)
      Rv[pt] = (fmpq**)malloc(nc * sizeof(fmpq*));
    weight_eval_values(Rv, mom, E.npts, BzT, W, S, h, nb, nc);

    for (c = 0; c < nc; c++) {
      if (Pj[c].rows == 0)
	continue;
      // Rsel[r][col], the value at the r-th point of E[c] of the image of the basis vector col, without p^a, and
      // Rchk the same at the points of the check; then the block is p^a Einv Rsel
      {
	fmpq_mat_t Rsel, Rchk, B, C;
	fmpq_mat_init(Rsel, Pj[c].rows, S->dims[c]);
	fmpq_mat_init(Rchk, WEIGHT_EVAL_CHECKS, S->dims[c]);
	for (r = 0; r < Pj[c].rows + WEIGHT_EVAL_CHECKS; r++) {
	  pt = (r < Pj[c].rows) ? E.sel[c][r] : E.npts - WEIGHT_EVAL_CHECKS + (r - Pj[c].rows);
	  for (col = 0; col < S->dims[c]; col++)
	    fmpq_set((r < Pj[c].rows) ? fmpq_mat_entry(Rsel, r, col) : fmpq_mat_entry(Rchk, r - Pj[c].rows, col),
		     Rv[pt][c] + col);
	}
	fmpq_mat_init(B, Pj[c].rows, S->dims[c]);
	fmpq_mat_mul(B, E.Einv + c, Rsel);
	// the check: F_U at the points of the check is the combination of the values there of the basis of
	// H[c][j]
	fmpq_mat_init(C, WEIGHT_EVAL_CHECKS, S->dims[c]);
	fmpq_mat_mul(C, E.chk + c, B);
	if (!fmpq_mat_equal(C, Rchk)) {
	  fprintf(stderr, "Error: at the conductor index %ld, the image of a form under a Hecke operator, found by evaluation, does not agree with its values at the points of the check, so it is not in the space of the forms of the lattice\n", (long)c);
	  exit(EXIT_FAILURE);
	}
	for (m = 0; m < Pj[c].rows; m++)
	  for (col = 0; col < S->dims[c]; col++)
	    fmpq_mul_fmpz(fmpq_mat_entry(T + c, S->offsets[c][j] + m, col), fmpq_mat_entry(B, m, col), pa);
	fmpq_mat_clear(Rsel);
	fmpq_mat_clear(Rchk);
	fmpq_mat_clear(B);
	fmpq_mat_clear(C);
      }
    }

    for (pt = 0; pt < E.npts; pt++) {
      for (c = 0; c < nc; c++)
	_fmpq_vec_clear(Rv[pt][c], (S->dims[c] > 0) ? S->dims[c] : 1);
      free(Rv[pt]);
      for (i = 0; i < h; i++) {
	for (b = 0; b < nb; b++)
	  weight_mom_clear(&(mom[pt][i][b]));
	free(mom[pt][i]);
      }
      free(mom[pt]);
    }
    free(Rv);
    free(mom);
    weight_eval_points_clear(&E);
    for (c = 0; c < nc; c++)
      weight_bpoly_clear(Pj + c);
  }

  if (any_rho)
    hecke_weight_rho_classes(T, genus, W, S, p, k, by_rho, true);
  free(by_rho);

  fmpz_mat_clear(BzT);
  free(Pj);
  fmpz_clear(pa);
  return T;
}

// Adds to the matrices T of T_{p^k} the rows of the classes j with classes[j] true, from the matrices of rho(s)
// on W: for k = 1 with use_orbits by the orbits of the automorphisms of L_j (weight_sums_orbits), and otherwise
// summed over every neighbor of nbr_data, each tested for isometry.
static void hecke_weight_rho_classes(fmpq_mat_struct* T, const genus_t genus, const weight_rep_t W,
				     const weight_space_t S, int p, int k, const bool* classes, bool use_orbits)
{
  fmpq_mat_t** acc;
  fmpq_mat_t R, M;
  nbr_data_t nbr_man;
  fmpz_mat_t nbr_fmpz, nbr_isom;
  square_matrix_t nbr;
  isometry_t s_nbr, hash_isom, s_inv;
  fmpz_t pa;
  slong c, i, j, h, nb, b;
  slong num_isom;
  double theta_time, isom_time;
  W64 bits;
  bool needed;

  h = genus->genus_reps->num_stored;
  nb = (slong)1 << genus->spinor->num_primes;

  fmpz_init(pa);
  // fmpz_ui_pow_ui is not in flint 2.9
  fmpz_set_ui(pa, p);
  // T_p is normalized by p^a and T_{p^2} by p^(2a). The eigenvalues are then algebraic integers (p^a
  // and p^(a + b) are the smallest powers which make them so), and the factor at p of the spin
  // L-function of a form is, with w = 2a + 3,
  //   1 - T_p x + p (T_{p^2} + p^(2a) (p^2 + 1)) x^2 - p^w T_p x^3 + p^(2w) x^4,
  // as for the weight 0
  fmpz_pow_ui(pa, pa, (k == 1) ? W->a : 2 * W->a);
  fmpq_mat_init(R, W->dim, W->dim);
  fmpq_mat_init(M, W->dim, W->dim);
  fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
  fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
  theta_time = isom_time = 0;
  num_isom = 0;

  // acc[i][b], the sum of p^a R(s) (p^(2a) R(s) for k = 2) over the p^k-neighbors of L_j isometric
  // to L_i by s with spinor bits b
  acc = (fmpq_mat_t**)malloc(h * sizeof(fmpq_mat_t*));
  for (i = 0; i < h; i++) {
    acc[i] = (fmpq_mat_t*)malloc(nb * sizeof(fmpq_mat_t));
    for (b = 0; b < nb; b++)
      fmpq_mat_init(acc[i][b], W->dim, W->dim);
  }

  for (j = 0; j < h; j++) {
    needed = false;
    for (c = 0; (c < genus->num_conductors) && classes[j]; c++)
      needed = needed || (fmpq_mat_nrows(S->H[c] + j) > 0);
    if (!needed)
      continue;
    for (i = 0; i < h; i++)
      for (b = 0; b < nb; b++)
	fmpq_mat_zero(acc[i][b]);

    // the orbits are taken for T_p, also at the primes dividing the discriminant, as for the Hecke
    // matrices of the weight 0 (hecke_col_all_conductors)
    if (use_orbits && (k == 1))
      weight_sums_orbits(acc, genus, W, j, p, pa);
    else {
    nbr_data_init(nbr_man, genus->genus_reps->keys[j], p, k);
    while (!nbr_data_has_ended(nbr_man)) {
      // the isometry from the first lattice to itself through this neighbor, as in
      // process_isotropic_vector_nbr_data_all_conductors
      nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
      square_matrix_set_fmpz_mat(nbr, nbr_fmpz);
      isometry_init_set_fmpz_mat(s_nbr, nbr_isom, p);
      i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, &theta_time, &isom_time, &num_isom);
      assert((0 <= i) && (i < h));
      isometry_muleq_left(s_nbr, genus->isoms[j]);
      isometry_muleq_right(s_nbr, hash_isom);
      isometry_inv(s_inv, genus->isoms[i]);
      isometry_muleq_right(s_nbr, s_inv);
      assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], genus->genus_reps->keys[0]));
      bits = spinor_norm_isom(genus->spinor, s_nbr);
      weight_rep_action(R, W, s_nbr);
      fmpq_mat_scalar_mul_fmpz(R, R, pa);
      fmpq_mat_add(acc[i][bits], acc[i][bits], R);
      nbr_data_get_next_neighbor(nbr_man);
    }
    nbr_data_clear(nbr_man);
    }

    for (c = 0; c < genus->num_conductors; c++) {
      if (fmpq_mat_nrows(S->H[c] + j) == 0)
	continue;
      for (i = 0; i < h; i++) {
	if (fmpq_mat_nrows(S->H[c] + i) == 0)
	  continue;
	fmpq_mat_zero(M);
	for (b = 0; b < nb; b++) {
	  if (weight_char((W64)b, c) > 0)
	    fmpq_mat_add(M, M, acc[i][b]);
	  else
	    fmpq_mat_sub(M, M, acc[i][b]);
	}
	weight_add_block(T + c, S, c, j, i, M);
      }
    }
  }

  for (i = 0; i < h; i++) {
    for (b = 0; b < nb; b++)
      fmpq_mat_clear(acc[i][b]);
    free(acc[i]);
  }
  free(acc);
  fmpq_mat_clear(R);
  fmpq_mat_clear(M);
  fmpz_mat_clear(nbr_fmpz);
  fmpz_mat_clear(nbr_isom);
  fmpz_clear(pa);
}

// The matrices of T_{p^k} by the matrices of rho(s) on W, with every neighbor of nbr_data tested for
// isometry: the definition, for the tests of hecke_matrices_weight_eval.
static fmpq_mat_struct* hecke_matrices_weight_impl(const genus_t genus, const weight_rep_t W,
						   const weight_space_t S, int p, int k)
{
  fmpq_mat_struct* T;
  bool* classes;
  slong c, j, h = genus->genus_reps->num_stored;

  T = (fmpq_mat_struct*)malloc(genus->num_conductors * sizeof(fmpq_mat_struct));
  for (c = 0; c < genus->num_conductors; c++)
    fmpq_mat_init(T + c, S->dims[c], S->dims[c]);
  classes = (bool*)malloc((h > 0 ? h : 1) * sizeof(bool));
  for (j = 0; j < h; j++)
    classes[j] = true;
  hecke_weight_rho_classes(T, genus, W, S, p, k, classes, false);
  free(classes);
  return T;
}

fmpq_mat_struct* hecke_matrices_weight(const genus_t genus, const weight_rep_t W, const weight_space_t S, int p)
{
  return hecke_matrices_weight_eval(genus, W, S, p, 1);
}

fmpq_mat_struct* hecke_matrices_weight_by_rho(const genus_t genus, const weight_rep_t W, const weight_space_t S, int p,
					      int k)
{
  return hecke_matrices_weight_impl(genus, W, S, p, k);
}

fmpq_mat_struct* hecke_matrices_weight_square(const genus_t genus, const weight_rep_t W, const weight_space_t S, int p)
{
  return hecke_matrices_weight_eval(genus, W, S, p, 2);
}

int print_hecke_matrices_weight(const genus_t genus, int a, int b, int p, int k)
{
  weight_rep_t W;
  weight_space_t S;
  fmpq_mat_struct* T;
  slong c, r, col;
  char* str;

  // a genus file may give no lattice: there are no forms
  if (genus->genus_reps->num_stored == 0) {
    printf("{");
    for (c = 0; c < genus->num_conductors; c++)
      printf("%ld : []%s", (long)genus->conductors[c], (c + 1 < genus->num_conductors) ? "," : "");
    printf("}");
    return 0;
  }
  weight_rep_init(W, a, b, genus->genus_reps->keys[0]);
  weight_space_init(S, genus, W);
  T = (k == 1) ? hecke_matrices_weight(genus, W, S, p) : hecke_matrices_weight_square(genus, W, S, p);

  if (b == 0)
    fprintf(stderr, "With the weight %d, the corresponding dimensions are: ", a);
  else
    fprintf(stderr, "With the weight (%d, %d), the corresponding dimensions are: ", a, b);
  for (c = 0; c < genus->num_conductors; c++)
    fprintf(stderr, "%ld ", (long)S->dims[c]);
  fprintf(stderr, "\n");

  printf("{");
  for (c = 0; c < genus->num_conductors; c++) {
    printf("%ld : [", (long)genus->conductors[c]);
    for (r = 0; r < fmpq_mat_nrows(T + c); r++) {
      printf("[");
      for (col = 0; col < fmpq_mat_ncols(T + c); col++) {
	str = fmpq_get_str(NULL, 10, fmpq_mat_entry(T + c, r, col));
	printf("%s%s", str, (col + 1 < fmpq_mat_ncols(T + c)) ? "," : "");
	flint_free(str);
      }
      printf("]%s", (r + 1 < fmpq_mat_nrows(T + c)) ? "," : "");
    }
    printf("]%s", (c + 1 < genus->num_conductors) ? "," : "");
  }
  // as compute_hecke_matrix_all_conds prints it, with no newline at the end
  printf("}");

  for (c = 0; c < genus->num_conductors; c++)
    fmpq_mat_clear(T + c);
  free(T);
  weight_space_clear(S);
  weight_rep_clear(W);
  return 0;
}

/* Eigenvalues of the forms of weight (a, b). */

#include "arith.h"
#include "decomposition.h"
#include "eigenvalues.h"
#include "hecke.h"

typedef struct
{
  const weight_rep_struct* W;
  const weight_space_struct* S;
  // the primes of the matrices the decomposition took, in the order of its list D->hecke
  int* ps;
  slong num_ps;
} weight_source_struct;

// The source of Hecke matrices of the decomposition for the forms of a weight: see decomposition.h.
// It records the prime of each matrix, so that the eigenvalues reuse the matrices of the decomposition.
static void weight_hecke_source(fmpq_mat_t* hecke, const genus_t genus, int p, void* data)
{
  weight_source_struct* d = (weight_source_struct*)data;
  fmpq_mat_struct* T;
  slong c;

  T = hecke_matrices_weight(genus, d->W, d->S, p);
  for (c = 0; c < genus->num_conductors; c++) {
    fmpq_mat_init_set(hecke[c], T + c);
    fmpq_mat_clear(T + c);
  }
  free(T);
  d->ps = (int*)realloc(d->ps, (d->num_ps + 1) * sizeof(int));
  d->ps[d->num_ps++] = p;
  return;
}

// Sets ev to the eigenvalue of T at the row vector v, v T = ev v, and returns true, if v is an
// eigenvector of T; returns false otherwise.
static bool weight_eigenvalue(nf_elem_t ev, const nf_elem_t* v, const fmpq_mat_t T, const nf_t nf)
{
  slong dim, j, k, k0;
  nf_elem_t* w;
  nf_elem_t tmp;
  bool is_eigen;

  dim = fmpq_mat_nrows(T);
  w = (nf_elem_t*)malloc((dim > 0 ? dim : 1) * sizeof(nf_elem_t));
  nf_elem_init(tmp, nf);
  for (k = 0; k < dim; k++) {
    nf_elem_init(w[k], nf);
    nf_elem_zero(w[k], nf);
    for (j = 0; j < dim; j++)
      if (!fmpq_is_zero(fmpq_mat_entry(T, j, k))) {
	nf_elem_scalar_mul_fmpq(tmp, v[j], fmpq_mat_entry(T, j, k), nf);
	nf_elem_add(w[k], w[k], tmp, nf);
      }
  }
  for (k0 = 0; (k0 < dim) && nf_elem_is_zero(v[k0], nf); k0++);
  is_eigen = (k0 < dim);
  if (is_eigen) {
    nf_elem_div(ev, w[k0], v[k0], nf);
    for (k = 0; is_eigen && (k < dim); k++) {
      nf_elem_mul(tmp, ev, v[k], nf);
      is_eigen = nf_elem_equal(tmp, w[k], nf);
    }
  }
  for (k = 0; k < dim; k++)
    nf_elem_clear(w[k], nf);
  free(w);
  nf_elem_clear(tmp, nf);
  return is_eigen;
}

// Sets ev (initialized in the field of the form i of evs) to the eigenvalue of T, the matrix of a
// Hecke operator on the whole space of the conductor, at the form, and returns true; returns false
// when it is not defined: when the vector of an old form is not an eigenvector of T, or when T is
// not an element of the field of the form on its space, so that the number would depend on the
// vector which was kept.
static bool weight_form_eigenvalue(nf_elem_t ev, const eigenvalues_t evs, slong i, const fmpq_mat_t T)
{
  if (evs->is_old[i] && !ev_is_scalar_on_old_space(evs, i, T))
    return false;
  return weight_eigenvalue(ev, (const nf_elem_t*)evs->eigenvecs[i], T, evs->nfs[i]);
}

// At a prime p which divides the discriminant exactly once, the eigenvalue t of T_p (normalized by
// p^a) of a form of conductor c is reported as lambda = t + p^a, negated when p divides c, as T_p + 1 is
// for the weight 0. lambda is the number for which the factor at p of the spin L-function of a newform
// is, with w = 2a + 3 and a sign eps,
//   (1 - (lambda + eps p^(a+1)) x + p^w x^2) (1 + eps p^(a+1) x),
// the shape of Roberts-Schmidt for the paramodular level p.
static void weight_bad_prime_eigenvalue(nf_elem_t ev, const nf_t nf, int p, int a, slong conductor)
{
  fmpz_t pa;

  fmpz_init_set_ui(pa, p);
  fmpz_pow_ui(pa, pa, a);
  nf_elem_add_fmpz(ev, ev, pa, nf);
  if (conductor % p == 0)
    nf_elem_neg(ev, ev, nf);
  fmpz_clear(pa);
}

// the source of the transposes of the Hecke matrices, for the right eigenvectors: from the matrices of the
// first decomposition when it has the prime
typedef struct
{
  const weight_rep_struct* W;
  const weight_space_struct* S;
  const decomposition_struct* D;
  const int* ps;
  slong num_ps;
  // the primes of the transposed matrices, in the order of the second decomposition
  int* tps;
  slong num_tps;
} weight_tsource_struct;

static void weight_hecke_tsource(fmpq_mat_t* hecke, const genus_t genus, int p, void* data)
{
  weight_tsource_struct* d = (weight_tsource_struct*)data;
  fmpq_mat_struct* T = NULL;
  slong c, k;

  d->tps = (int*)realloc(d->tps, (d->num_tps + 1) * sizeof(int));
  d->tps[d->num_tps++] = p;
  for (k = 0; (k < d->num_ps) && (d->ps[k] != p); k++);
  if (k == d->num_ps)
    T = hecke_matrices_weight(genus, d->W, d->S, p);
  for (c = 0; c < genus->num_conductors; c++) {
    const fmpq_mat_struct* M = (T == NULL) ? &(d->D->hecke[k][c][0]) : T + c;
    fmpq_mat_init(hecke[c], fmpq_mat_ncols(M), fmpq_mat_nrows(M));
    fmpq_mat_transpose(hecke[c], M);
  }
  if (T != NULL) {
    for (c = 0; c < genus->num_conductors; c++)
      fmpq_mat_clear(T + c);
    free(T);
  }
}

// The data of the evaluation for a form: its right eigenvector u (of dim entries), its pivot class, and
// phi = F_{U_pivot}(x0, y0), which is not zero.
//
// An old form (m > 0) has no eigenvector, and its eigenvalue needs the operator T_W on its space W (the rows
// of old_basis, W T = T_W W; see ev_is_scalar_on_old_space). The spaces of all the forms of the conductor
// are preserved by the Hecke operators, which commute, and their sum is the whole space; so the columns of
// the solution U of P U = E, P the rows of all of them and E the columns of the unit matrix at the rows of
// W, satisfy T U = U T_W, and W U = 1. A row b^t of the values at a point z of the polynomials of the basis
// of H[c][j] (zero at the other classes) gives b^t T, pa times the sums of weight_eval_pairing over the
// neighbors of L_j at z, so that b^t T U = (b^t U) T_W. With m points z_r, at the classes pcls[r] and of the
// indices pid[r] of weight_eval_candidate, at which the rows b_r^t U are independent, T_W = Q^(-1) R, Q[r] =
// b_r^t U and R[r] = b_r^t T U. WEIGHT_EVAL_CHECKS more points, of pid -1, -2, ... (weight_eval_check_point)
// at the first and the last of these classes, check that T U = U T_W there. coords = v U are the coordinates of the vector v of
// the form in the basis W.
typedef struct
{
  bool ok;
  slong pivot;
  nf_elem_t phi;
  nf_elem_t* u;
  slong dim;
  slong m;
  fmpq_mat_t U, Qinv, chkQ, R;
  nf_elem_t* coords;
  slong *pcls, *pid;
  // at each prime: 1 when the eigenvalue was found by the evaluation, 2 when T is not a scalar on the space
  // of the old form, 0 when it is left to the matrices
  int st;
} weight_eval_struct;

// the point of index t of the evaluation: weight_eval_candidate for t >= 0, weight_eval_check_point(-t) else
static void weight_eval_point_of(fmpz* X, fmpz* Y, slong t)
{
  if (t >= 0)
    weight_eval_candidate(X, Y, t);
  else
    weight_eval_check_point(X, Y, (int)(-t));
}

// q[l] = sum_m v[m] U[off + m][l], for the values v of the polynomials of a class at a point
static void weight_old_row(fmpq* q, const fmpq* v, slong rows, const fmpq_mat_t U, slong off)
{
  slong l, m;
  fmpq_t t;

  fmpq_init(t);
  for (l = 0; l < fmpq_mat_ncols(U); l++) {
    fmpq_zero(q + l);
    for (m = 0; m < rows; m++) {
      if (fmpq_is_zero(v + m) || fmpq_is_zero(fmpq_mat_entry(U, off + m, l)))
	continue;
      fmpq_mul(t, v + m, fmpq_mat_entry(U, off + m, l));
      fmpq_add(q + l, q + l, t);
    }
  }
  fmpq_clear(t);
}

// Sets the data of the evaluation of the old forms of the conductor index c (see weight_eval_struct), whose
// E[i].m is 0 and E[i].ok false on entry. The classes of the points are those of chosen first, the pivots of
// the new forms, at which the moments are summed anyway, and then by the order of their automorphism groups.
// A form whose points are not found stays without the data, and its eigenvalues come from the matrices.
static void weight_old_eval_init(weight_eval_struct* E, const eigenvalues_t evs, const decomposition_t D,
				 slong c, const weight_space_t S, const weight_rep_t W, const weight_bpoly_struct* P,
				 const bool* chosen, const slong* aut_order, slong h)
{
  slong i, j, jj, l, r, t, num_old, tot, dim, row, col;
  slong *ord, *row_off, *col_off;
  fmpq_mat_t PM, B, X;
  fmpz *vec, *xm, *ym;
  fmpz PX[QF_RANK], PY[QF_RANK];
  fmpq* vals;

  dim = S->dims[c];
  if ((dim == 0) || (D->num[c] != evs->num))
    return;
  num_old = 0;
  tot = 0;
  for (i = 0; i < evs->num; i++) {
    tot += fmpq_mat_nrows(D->bases[c][i]);
    if (evs->is_old[i])
      num_old += fmpq_mat_nrows(evs->old_basis[i]);
  }
  if ((num_old == 0) || (tot != dim))
    return;

  // P U = E for all the old forms at once
  row_off = (slong*)malloc(evs->num * sizeof(slong));
  col_off = (slong*)malloc(evs->num * sizeof(slong));
  fmpq_mat_init(PM, dim, dim);
  fmpq_mat_init(B, dim, num_old);
  fmpq_mat_init(X, dim, num_old);
  row = col = 0;
  for (i = 0; i < evs->num; i++) {
    slong nr = fmpq_mat_nrows(D->bases[c][i]);
    row_off[i] = row;
    col_off[i] = col;
    for (r = 0; r < nr; r++)
      for (l = 0; l < dim; l++)
	fmpq_set(fmpq_mat_entry(PM, row + r, l), fmpq_mat_entry(D->bases[c][i], r, l));
    if (evs->is_old[i]) {
      for (r = 0; r < nr; r++)
	fmpq_one(fmpq_mat_entry(B, row + r, col + r));
      col += nr;
    }
    row += nr;
  }
  if (!fmpq_mat_solve_fraction_free(X, PM, B)) {
    fmpq_mat_clear(PM);
    fmpq_mat_clear(B);
    fmpq_mat_clear(X);
    free(row_off);
    free(col_off);
    return;
  }
  fmpq_mat_clear(PM);
  fmpq_mat_clear(B);

  // the order of the classes
  ord = (slong*)malloc(h * sizeof(slong));
  for (j = 0; j < h; j++)
    ord[j] = j;
  for (j = 1; j < h; j++)
    for (jj = j; (jj > 0) && ((chosen[ord[jj]] && !chosen[ord[jj - 1]]) ||
			      ((chosen[ord[jj]] == chosen[ord[jj - 1]]) && (aut_order[ord[jj]] > aut_order[ord[jj - 1]])));
	 jj--) {
      t = ord[jj];
      ord[jj] = ord[jj - 1];
      ord[jj - 1] = t;
    }

  vec = _fmpz_vec_init(W->nx * W->ny);
  xm = _fmpz_vec_init(W->nx);
  ym = _fmpz_vec_init(W->ny);
  for (r = 0; r < QF_RANK; r++) {
    fmpz_init(PX + r);
    fmpz_init(PY + r);
  }
  vals = _fmpq_vec_init(dim);
  for (i = 0; i < evs->num; i++) {
    weight_eval_struct* F = &(E[i]);
    slong m, got;
    fmpq_mat_t Q, ech, tmp;
    if (!evs->is_old[i])
      continue;
    m = fmpq_mat_nrows(evs->old_basis[i]);
    fmpq_mat_init(F->U, dim, m);
    for (r = 0; r < dim; r++)
      for (l = 0; l < m; l++)
	fmpq_set(fmpq_mat_entry(F->U, r, l), fmpq_mat_entry(X, r, col_off[i] + l));
    F->pcls = (slong*)malloc((m + WEIGHT_EVAL_CHECKS) * sizeof(slong));
    F->pid = (slong*)malloc((m + WEIGHT_EVAL_CHECKS) * sizeof(slong));
    fmpq_mat_init(Q, m, m);
    got = 0;
    for (j = 0; (j < h) && (got < m); j++) {
      slong rows;
      jj = ord[j];
      rows = P[jj].rows;
      for (t = 0; (t < 2 * rows + 16) && (got < m) && (rows > 0); t++) {
	slong rank;
	weight_eval_point_of(PX, PY, t);
	weight_bpoly_values(vals, &(P[jj]), W, PX, PY, vec, xm, ym);
	weight_old_row(fmpq_mat_entry(Q, got, 0), vals, rows, F->U, S->offsets[c][jj]);
	fmpq_mat_init(ech, got + 1, m);
	fmpq_mat_init(tmp, got + 1, m);
	for (r = 0; r <= got; r++)
	  for (l = 0; l < m; l++)
	    fmpq_set(fmpq_mat_entry(ech, r, l), fmpq_mat_entry(Q, r, l));
	rank = fmpq_mat_rref(tmp, ech);
	fmpq_mat_clear(ech);
	fmpq_mat_clear(tmp);
	if (rank == got + 1) {
	  F->pcls[got] = jj;
	  F->pid[got] = t;
	  got++;
	}
      }
    }
    if (got < m) {
      fmpq_mat_clear(Q);
      fmpq_mat_clear(F->U);
      free(F->pcls);
      free(F->pid);
      continue;
    }
    fmpq_mat_init(F->Qinv, m, m);
    fmpq_mat_inv(F->Qinv, Q);
    fmpq_mat_clear(Q);
    fmpq_mat_init(F->chkQ, WEIGHT_EVAL_CHECKS, m);
    for (r = 0; r < WEIGHT_EVAL_CHECKS; r++) {
      // at the first and the last of the classes of the points, when they differ
      jj = F->pcls[(WEIGHT_EVAL_CHECKS > 1) ? (r * (m - 1)) / (WEIGHT_EVAL_CHECKS - 1) : 0];
      F->pcls[m + r] = jj;
      F->pid[m + r] = -(r + 1);
      weight_eval_point_of(PX, PY, -(r + 1));
      weight_bpoly_values(vals, &(P[jj]), W, PX, PY, vec, xm, ym);
      weight_old_row(fmpq_mat_entry(F->chkQ, r, 0), vals, P[jj].rows, F->U, S->offsets[c][jj]);
    }
    fmpq_mat_init(F->R, m + WEIGHT_EVAL_CHECKS, m);
    F->coords = (nf_elem_t*)malloc(m * sizeof(nf_elem_t));
    for (l = 0; l < m; l++)
      nf_elem_init(F->coords[l], evs->nfs[i]);
    nf_elem_vec_mul_fmpq_mat(F->coords, (const nf_elem_t*)evs->eigenvecs[i], F->U, evs->nfs[i]);
    F->m = m;
    F->pivot = -1;
    F->ok = true;
  }
  _fmpq_vec_clear(vals, dim);
  for (r = 0; r < QF_RANK; r++) {
    fmpz_clear(PX + r);
    fmpz_clear(PY + r);
  }
  _fmpz_vec_clear(vec, W->nx * W->ny);
  _fmpz_vec_clear(xm, W->nx);
  _fmpz_vec_clear(ym, W->ny);
  free(ord);
  fmpq_mat_clear(X);
  free(row_off);
  free(col_off);
}

static void weight_eval_clear(weight_eval_struct* E, const nf_t nf)
{
  slong m;

  if (E->u != NULL) {
    for (m = 0; m < E->dim; m++)
      nf_elem_clear(E->u[m], nf);
    free(E->u);
  }
  nf_elem_clear(E->phi, nf);
  if (E->m > 0) {
    for (m = 0; m < E->m; m++)
      nf_elem_clear(E->coords[m], nf);
    free(E->coords);
    fmpq_mat_clear(E->U);
    fmpq_mat_clear(E->Qinv);
    fmpq_mat_clear(E->chkQ);
    fmpq_mat_clear(E->R);
    free(E->pcls);
    free(E->pid);
  }
}

int print_eigenvalues_weight(const genus_t genus, int a, int b, const int* ps, int num_ps, const int* ps2,
			     int num_ps2, bool set_types)
{
  weight_rep_t W;
  weight_space_t S;
  weight_source_struct data;
  decomposition_t D;
  eigenvalues_t* evs;
  fmpq_mat_struct* Tp;
  // evals[c][i * nq + n]: the eigenvalue of the form i of conductor c of T_q for n < xnum_ps, and of
  // T_{q^2} otherwise, where q = xps[n] or xps2[n - xnum_ps], if has_ev[c][...]
  nf_elem_t** evals;
  bool** has_ev;
  slong c, i, k, nq;
  int n, m, q, num_type_ps, xnum_ps, xnum_ps2;
  int *xps, *xps2;
  bool once, found;
  aut_type t;
  // the evaluation: ev_data[c][i] for the form i of conductor c, and the orders of the automorphism groups
  weight_eval_struct** ev_data;
  weight_bpoly_struct** bpoly;
  slong* aut_order;
  slong h, nb, len, jj, i2, kk, k2;

  // The primes at which the eigenvalues are computed: those of ps and ps2, which are printed, and, for
  // the types, the first primes not dividing the discriminant after those of ps2, until there are
  // LIFT_TYPE_MIN_PRIMES of them at which T_p and T_{p^2} are known. The eigenvalues at these extra
  // primes are not printed. As for the weight 0, which computes them in eigenvalues_set_lifts.
  xps = (int*)malloc((num_ps + LIFT_TYPE_MIN_PRIMES) * sizeof(int));
  xps2 = (int*)malloc((num_ps2 + LIFT_TYPE_MIN_PRIMES) * sizeof(int));
  for (n = 0; n < num_ps; n++)
    xps[n] = ps[n];
  for (n = 0; n < num_ps2; n++)
    xps2[n] = ps2[n];
  xnum_ps = num_ps;
  xnum_ps2 = num_ps2;
  num_type_ps = 0;
  if (set_types) {
    for (n = 0; n < num_ps2; n++)
      if (fmpz_fdiv_ui(genus->disc, ps2[n]) != 0)
	num_type_ps++;
    for (q = 2; num_type_ps < LIFT_TYPE_MIN_PRIMES; q = (int)n_nextprime(q, 1)) {
      if (fmpz_fdiv_ui(genus->disc, q) == 0)
	continue;
      for (found = false, n = 0; n < xnum_ps2; n++)
	found = found || (xps2[n] == q);
      if (found)
	continue;
      xps2[xnum_ps2++] = q;
      for (found = false, n = 0; n < xnum_ps; n++)
	found = found || (xps[n] == q);
      if (!found)
	xps[xnum_ps++] = q;
      num_type_ps++;
    }
  }

  nq = xnum_ps + xnum_ps2;
  if (genus->genus_reps->num_stored == 0) {
    printf("{");
    for (c = 0; c < genus->num_conductors; c++)
      printf("%ld : [], ", (long)genus->conductors[c]);
    printf("}");
    free(xps);
    free(xps2);
    return 0;
  }

  weight_rep_init(W, a, b, genus->genus_reps->keys[0]);
  weight_space_init(S, genus, W);
  if (b == 0)
    fprintf(stderr, "With the weight %d, the corresponding dimensions are: ", a);
  else
    fprintf(stderr, "With the weight (%d, %d), the corresponding dimensions are: ", a, b);
  for (c = 0; c < genus->num_conductors; c++)
    fprintf(stderr, "%ld ", (long)S->dims[c]);
  fprintf(stderr, "\n");

  // the eigenforms, from the decomposition of the spaces of the weight
  data.W = W;
  data.S = S;
  data.ps = NULL;
  data.num_ps = 0;
  decomposition_init(D, genus->num_conductors);
  D->source = weight_hecke_source;
  D->source_data = &data;
  D->dims = S->dims;
  evs = (eigenvalues_t*)malloc(genus->num_conductors * sizeof(eigenvalues_t));
  for (c = 0; c < genus->num_conductors; c++) {
    decompose(D, genus, c);
    hecke_eigenforms(evs[c], D, genus, c);
    // the types are set below, from the eigenvalues
    for (i = 0; i < evs[c]->num; i++)
      if (evs[c]->lift_type[i] != O)
	evs[c]->lift_type[i] = U;
  }
  assert(data.num_ps == D->num_primes);

  // The right eigenvectors of the forms, from the decomposition of the transposed matrices, and the data
  // of the evaluation of each form whose right eigenvector is found (see weight_moments).
  h = genus->genus_reps->num_stored;
  nb = (slong)1 << genus->spinor->num_primes;
  len = W->nx * W->ny;
  aut_order = (slong*)malloc(h * sizeof(slong));
  for (jj = 0; jj < h; jj++) {
    aut_grp_t grp;
    aut_grp_init_square_matrix(grp, genus->genus_reps->keys[jj]);
    aut_order[jj] = grp->order;
    aut_grp_clear(grp);
  }
  ev_data = (weight_eval_struct**)malloc(genus->num_conductors * sizeof(weight_eval_struct*));
  bpoly = (weight_bpoly_struct**)malloc(genus->num_conductors * sizeof(weight_bpoly_struct*));
  {
    weight_tsource_struct tdata;
    decomposition_t D2;
    eigenvalues_t* evs2;
    bool* chosen;
    // the transposes of the matrices of D at the primes which D2 did not use, for the conductor c
    fmpq_mat_struct* tr;
    bool* has_tr;
    fmpz* mon0;
    fmpz X0[QF_RANK], Y0[QF_RANK];
    fmpz *xm, *ym;
    fmpz_t dot;
    fmpq_t q0;
    nf_elem_t tmp;
    int r;

    tdata.W = W;
    tdata.S = S;
    tdata.D = D;
    tdata.ps = data.ps;
    tdata.num_ps = data.num_ps;
    tdata.tps = NULL;
    tdata.num_tps = 0;
    decomposition_init(D2, genus->num_conductors);
    D2->source = weight_hecke_tsource;
    D2->source_data = &tdata;
    D2->dims = S->dims;
    evs2 = (eigenvalues_t*)malloc(genus->num_conductors * sizeof(eigenvalues_t));
    // the bimonomials at (x0, y0)
    mon0 = _fmpz_vec_init(len);
    xm = _fmpz_vec_init(W->nx);
    ym = _fmpz_vec_init(W->ny);
    for (r = 0; r < QF_RANK; r++) {
      fmpz_init_set_si(X0 + r, weight_x0[r]);
      fmpz_init_set_si(Y0 + r, weight_y0[r]);
    }
    weight_mom_values(mon0, W, X0, Y0, xm, ym);
    fmpz_init(dot);
    fmpq_init(q0);
    chosen = (bool*)calloc(h, sizeof(bool));
    // the eigenvectors of the transposes are the same work again: their polynomials are not printed
    matrix_tools_quiet = true;
    for (c = 0; c < genus->num_conductors; c++) {
      fmpq* bphi;
      decompose(D2, genus, c);
      hecke_eigenforms(evs2[c], D2, genus, c);
      bpoly[c] = (weight_bpoly_struct*)malloc(h * sizeof(weight_bpoly_struct));
      for (jj = 0; jj < h; jj++)
	weight_bpoly_init(&(bpoly[c][jj]), S->H[c] + jj, W);
      // bphi[offsets[c][j] + m] = B_j[m](x0, y0)
      bphi = _fmpq_vec_init((S->dims[c] > 0) ? S->dims[c] : 1);
      for (jj = 0; jj < h; jj++) {
	slong m;
	for (m = 0; m < bpoly[c][jj].rows; m++) {
	  slong mono;
	  fmpz_zero(dot);
	  for (mono = 0; mono < len; mono++)
	    fmpz_addmul(dot, fmpz_mat_entry(bpoly[c][jj].num, m, mono), mon0 + mono);
	  fmpq_set_fmpz_frac(bphi + S->offsets[c][jj] + m, dot, bpoly[c][jj].den + m);
	}
      }
      ev_data[c] = (weight_eval_struct*)malloc((evs[c]->num > 0 ? evs[c]->num : 1) * sizeof(weight_eval_struct));
      tr = (fmpq_mat_struct*)malloc((D->num_primes > 0 ? D->num_primes : 1) * sizeof(fmpq_mat_struct));
      has_tr = (bool*)calloc((D->num_primes > 0 ? D->num_primes : 1), sizeof(bool));
      for (i = 0; i < evs[c]->num; i++) {
	weight_eval_struct* E = &(ev_data[c][i]);
	const nf_struct* nf = evs[c]->nfs[i];
	slong match = -1;
	nf_elem_t lam_l, lam_r;
	E->ok = false;
	E->u = NULL;
	E->dim = S->dims[c];
	E->m = 0;
	nf_elem_init(E->phi, nf);
	if (evs[c]->is_old[i])
	  continue;
	nf_elem_init(lam_l, nf);
	nf_elem_init(lam_r, nf);
	// the form of the transposes with the same field and the same eigenvalues at the primes of D
	for (i2 = 0; (i2 < evs2[c]->num) && (match < 0); i2++) {
	  bool same;
	  if (evs2[c]->is_old[i2] || !fmpq_poly_equal(evs2[c]->nfs[i2]->pol, nf->pol))
	    continue;
	  same = true;
	  for (kk = 0; same && (kk < D->num_primes); kk++) {
	    const fmpq_mat_struct* Tr;
	    for (k2 = 0; (k2 < tdata.num_tps) && (tdata.tps[k2] != data.ps[kk]); k2++);
	    if (k2 < tdata.num_tps)
	      Tr = &(D2->hecke[k2][c][0]);
	    else {
	      // a prime which the transposes of this conductor did not need (D took it for a later conductor,
	      // or the transposes were split sooner): the transpose of the matrix of D
	      if (!has_tr[kk]) {
		fmpq_mat_init(tr + kk, fmpq_mat_ncols(&(D->hecke[kk][c][0])), fmpq_mat_nrows(&(D->hecke[kk][c][0])));
		fmpq_mat_transpose(tr + kk, &(D->hecke[kk][c][0]));
		has_tr[kk] = true;
	      }
	      Tr = tr + kk;
	    }
	    same = weight_eigenvalue(lam_l, (const nf_elem_t*)evs[c]->eigenvecs[i], &(D->hecke[kk][c][0]), nf) &&
	      weight_eigenvalue(lam_r, (const nf_elem_t*)evs2[c]->eigenvecs[i2], Tr, nf) &&
	      nf_elem_equal(lam_l, lam_r, nf);
	  }
	  if (same)
	    match = i2;
	}
	nf_elem_clear(lam_l, nf);
	nf_elem_clear(lam_r, nf);
	if (match < 0)
	  continue;
	// the pivot: a class at which phi != 0, one chosen already if there is one, with the most automorphisms
	{
	  const nf_elem_t* u = (const nf_elem_t*)evs2[c]->eigenvecs[match];
	  slong best = -1, m;
	  bool best_chosen = false;
	  nf_elem_t ph;
	  nf_elem_init(tmp, nf);
	  nf_elem_init(ph, nf);
	  for (jj = 0; jj < h; jj++) {
	    slong off = S->offsets[c][jj];
	    nf_elem_zero(ph, nf);
	    for (m = 0; m < bpoly[c][jj].rows; m++) {
	      if (nf_elem_is_zero(u[off + m], nf) || fmpq_is_zero(bphi + off + m))
		continue;
	      nf_elem_scalar_mul_fmpq(tmp, u[off + m], bphi + off + m, nf);
	      nf_elem_add(ph, ph, tmp, nf);
	    }
	    if (nf_elem_is_zero(ph, nf))
	      continue;
	    if ((best < 0) || (chosen[jj] && !best_chosen) ||
		((chosen[jj] == best_chosen) && (aut_order[jj] > aut_order[best]))) {
	      best = jj;
	      best_chosen = chosen[jj];
	      nf_elem_set(E->phi, ph, nf);
	    }
	  }
	  nf_elem_clear(ph, nf);
	  nf_elem_clear(tmp, nf);
	  if (best >= 0) {
	    E->ok = true;
	    E->pivot = best;
	    chosen[best] = true;
	    E->u = (nf_elem_t*)malloc((E->dim > 0 ? E->dim : 1) * sizeof(nf_elem_t));
	    for (m = 0; m < E->dim; m++) {
	      nf_elem_init(E->u[m], nf);
	      nf_elem_set(E->u[m], u[m], nf);
	    }
	  }
	}
      }
      _fmpq_vec_clear(bphi, (S->dims[c] > 0) ? S->dims[c] : 1);
      for (kk = 0; kk < D->num_primes; kk++)
	if (has_tr[kk])
	  fmpq_mat_clear(tr + kk);
      free(tr);
      free(has_tr);
    }
    // the old forms, once the pivots of the new forms of every conductor are chosen
    for (c = 0; c < genus->num_conductors; c++)
      weight_old_eval_init(ev_data[c], evs[c], D, c, S, W, bpoly[c], chosen, aut_order, h);
    matrix_tools_quiet = false;
    for (c = 0; c < genus->num_conductors; c++)
      eigenvalues_clear(evs2[c]);
    free(evs2);
    decomposition_clear(D2);
    free(tdata.tps);
    free(chosen);
    fmpz_clear(dot);
    fmpq_clear(q0);
    _fmpz_vec_clear(mon0, len);
    _fmpz_vec_clear(xm, W->nx);
    _fmpz_vec_clear(ym, W->ny);
    for (r = 0; r < QF_RANK; r++) {
      fmpz_clear(X0 + r);
      fmpz_clear(Y0 + r);
    }
  }

  // the eigenvalues, one prime at a time: the matrices of the decomposition are used again, and any
  // other is computed and cleared once its eigenvalues are found
  evals = (nf_elem_t**)malloc(genus->num_conductors * sizeof(nf_elem_t*));
  has_ev = (bool**)malloc(genus->num_conductors * sizeof(bool*));
  for (c = 0; c < genus->num_conductors; c++) {
    k = (evs[c]->num * nq > 0) ? evs[c]->num * nq : 1;
    evals[c] = (nf_elem_t*)malloc(k * sizeof(nf_elem_t));
    has_ev[c] = (bool*)malloc(k * sizeof(bool));
    for (i = 0; i < evs[c]->num * nq; i++)
      nf_elem_init(evals[c][i], evs[c]->nfs[i / nq]);
  }
  for (n = 0; n < nq; n++) {
    bool by_eval, need_full, old_eval;
    q = (n < xnum_ps) ? xps[n] : xps2[n - xnum_ps];
    if (n < xnum_ps)
      for (k = 0; (k < data.num_ps) && (data.ps[k] != q); k++);
    else
      k = data.num_ps;
    // the matrices of the decomposition are used at its primes; elsewhere the forms with the data of the
    // evaluation take their eigenvalues from it, and the others from the matrices
    by_eval = (k == data.num_ps);
    need_full = !by_eval;
    for (c = 0; (c < genus->num_conductors) && !need_full; c++)
      for (i = 0; (i < evs[c]->num) && !need_full; i++)
	need_full = !ev_data[c][i].ok;
    if (n >= xnum_ps)
      // the decomposition has only matrices of T_p, so k is not used then
      k = 0;
    if (!need_full)
      Tp = NULL;
    else if (n < xnum_ps)
      Tp = (!by_eval) ? NULL : hecke_matrices_weight(genus, W, S, q);
    else
      Tp = hecke_matrices_weight_square(genus, W, S, q);
    // the old forms take their eigenvalues from the evaluation when the matrices are not computed
    old_eval = (Tp == NULL);
    for (c = 0; c < genus->num_conductors; c++)
      for (i = 0; i < evs[c]->num; i++)
	ev_data[c][i].st = 0;
    if (by_eval) {
      fmpz_t pa;
      fmpz_init(pa);
      fmpz_set_ui(pa, q);
      fmpz_pow_ui(pa, pa, (n < xnum_ps) ? a : 2 * a);
      for (jj = 0; jj < h; jj++) {
	weight_mom_struct*** mom;
	fmpq** rrs;
	fmpz *PX, *PY;
	slong *ids, nids, ii, bb, pt, r, l;
	// the indices of the points at the class jj: 0, that is (x0, y0), for the new forms of pivot jj, and
	// those of the old forms
	ids = NULL;
	nids = 0;
	for (c = 0; c < genus->num_conductors; c++)
	  for (i = 0; i < evs[c]->num; i++) {
	    weight_eval_struct* E = &(ev_data[c][i]);
	    slong np = (E->m > 0) ? E->m + WEIGHT_EVAL_CHECKS : 1;
	    if (!E->ok || ((E->m > 0) && !old_eval))
	      continue;
	    for (r = 0; r < np; r++) {
	      slong id;
	      if ((E->m > 0) ? (E->pcls[r] != jj) : (E->pivot != jj))
		continue;
	      id = (E->m > 0) ? E->pid[r] : 0;
	      for (pt = 0; (pt < nids) && (ids[pt] != id); pt++);
	      if (pt == nids) {
		ids = (slong*)realloc(ids, (nids + 1) * sizeof(slong));
		ids[nids++] = id;
	      }
	    }
	  }
	if (nids == 0)
	  continue;
	PX = _fmpz_vec_init(QF_RANK * nids);
	PY = _fmpz_vec_init(QF_RANK * nids);
	for (pt = 0; pt < nids; pt++)
	  weight_eval_point_of(PX + QF_RANK * pt, PY + QF_RANK * pt, ids[pt]);
	mom = (weight_mom_struct***)malloc(nids * sizeof(weight_mom_struct**));
	for (pt = 0; pt < nids; pt++) {
	  mom[pt] = (weight_mom_struct**)malloc(h * sizeof(weight_mom_struct*));
	  for (ii = 0; ii < h; ii++) {
	    mom[pt][ii] = (weight_mom_struct*)malloc(nb * sizeof(weight_mom_struct));
	    for (bb = 0; bb < nb; bb++)
	      weight_mom_init(&(mom[pt][ii][bb]), len);
	  }
	}
	weight_moments(mom, genus, W, jj, q, (n < xnum_ps) ? 1 : 2, nids, PX, PY);
	rrs = (fmpq**)malloc(nids * sizeof(fmpq*));
	for (c = 0; c < genus->num_conductors; c++) {
	  slong dc = (S->dims[c] > 0) ? S->dims[c] : 1;
	  for (pt = 0; pt < nids; pt++)
	    rrs[pt] = NULL;
	  for (i = 0; i < evs[c]->num; i++) {
	    weight_eval_struct* E = &(ev_data[c][i]);
	    const nf_struct* nf = evs[c]->nfs[i];
	    nf_elem_t tmp;
	    slong m;
	    if (!E->ok || ((E->m > 0) && !old_eval))
	      continue;
	    if (E->m > 0) {
	      // R[r] = pa (b_r^t T / pa) U at the points of the class jj
	      for (r = 0; r < E->m + WEIGHT_EVAL_CHECKS; r++) {
		fmpq_t t;
		if (E->pcls[r] != jj)
		  continue;
		for (pt = 0; ids[pt] != E->pid[r]; pt++);
		if (rrs[pt] == NULL) {
		  rrs[pt] = _fmpq_vec_init(dc);
		  weight_eval_pairing(rrs[pt], bpoly[c], S, c, mom[pt], h, nb, len);
		}
		fmpq_init(t);
		for (l = 0; l < E->m; l++) {
		  fmpq* e = fmpq_mat_entry(E->R, r, l);
		  fmpq_zero(e);
		  for (m = 0; m < E->dim; m++) {
		    if (fmpq_is_zero(rrs[pt] + m) || fmpq_is_zero(fmpq_mat_entry(E->U, m, l)))
		      continue;
		    fmpq_mul(t, rrs[pt] + m, fmpq_mat_entry(E->U, m, l));
		    fmpq_add(e, e, t);
		  }
		  fmpq_mul_fmpz(e, e, pa);
		}
		fmpq_clear(t);
	      }
	      continue;
	    }
	    if (E->pivot != jj)
	      continue;
	    for (pt = 0; ids[pt] != 0; pt++);
	    if (rrs[pt] == NULL) {
	      rrs[pt] = _fmpq_vec_init(dc);
	      weight_eval_pairing(rrs[pt], bpoly[c], S, c, mom[pt], h, nb, len);
	    }
	    // lambda = pa <u, r> / phi
	    nf_elem_init(tmp, nf);
	    nf_elem_zero(evals[c][i * nq + n], nf);
	    for (m = 0; m < E->dim; m++) {
	      if (fmpq_is_zero(rrs[pt] + m) || nf_elem_is_zero(E->u[m], nf))
		continue;
	      nf_elem_scalar_mul_fmpq(tmp, E->u[m], rrs[pt] + m, nf);
	      nf_elem_add(evals[c][i * nq + n], evals[c][i * nq + n], tmp, nf);
	    }
	    nf_elem_scalar_mul_fmpz(evals[c][i * nq + n], evals[c][i * nq + n], pa, nf);
	    nf_elem_div(evals[c][i * nq + n], evals[c][i * nq + n], E->phi, nf);
	    nf_elem_clear(tmp, nf);
	    E->st = 1;
	  }
	  for (pt = 0; pt < nids; pt++)
	    if (rrs[pt] != NULL)
	      _fmpq_vec_clear(rrs[pt], dc);
	}
	free(rrs);
	for (pt = 0; pt < nids; pt++) {
	  for (ii = 0; ii < h; ii++) {
	    for (bb = 0; bb < nb; bb++)
	      weight_mom_clear(&(mom[pt][ii][bb]));
	    free(mom[pt][ii]);
	  }
	  free(mom[pt]);
	}
	free(mom);
	_fmpz_vec_clear(PX, QF_RANK * nids);
	_fmpz_vec_clear(PY, QF_RANK * nids);
	free(ids);
      }
      fmpz_clear(pa);
      // the old forms: T_W = Q^(-1) R at the first m points, checked at the others
      if (old_eval)
	for (c = 0; c < genus->num_conductors; c++)
	  for (i = 0; i < evs[c]->num; i++) {
	    weight_eval_struct* E = &(ev_data[c][i]);
	    fmpq_mat_t R0, TW, RC, QT;
	    slong r, l;
	    if (!E->ok || (E->m == 0))
	      continue;
	    fmpq_mat_init(R0, E->m, E->m);
	    fmpq_mat_init(RC, WEIGHT_EVAL_CHECKS, E->m);
	    for (l = 0; l < E->m; l++) {
	      for (r = 0; r < E->m; r++)
		fmpq_set(fmpq_mat_entry(R0, r, l), fmpq_mat_entry(E->R, r, l));
	      for (r = 0; r < WEIGHT_EVAL_CHECKS; r++)
		fmpq_set(fmpq_mat_entry(RC, r, l), fmpq_mat_entry(E->R, E->m + r, l));
	    }
	    fmpq_mat_init(TW, E->m, E->m);
	    fmpq_mat_mul(TW, E->Qinv, R0);
	    fmpq_mat_init(QT, WEIGHT_EVAL_CHECKS, E->m);
	    fmpq_mat_mul(QT, E->chkQ, TW);
	    // a check which fails leaves the form to the matrices
	    if (fmpq_mat_equal(QT, RC))
	      E->st = ev_old_space_eigenvalue(evals[c][i * nq + n], evs[c], i, TW,
					      (const nf_elem_t*)E->coords) ? 1 : 2;
	    fmpq_mat_clear(QT);
	    fmpq_mat_clear(TW);
	    fmpq_mat_clear(RC);
	    fmpq_mat_clear(R0);
	  }
    }
    // the matrices, for a check of an old form which failed
    for (c = 0; (c < genus->num_conductors) && by_eval && (Tp == NULL); c++)
      for (i = 0; (i < evs[c]->num) && (Tp == NULL); i++)
	if (ev_data[c][i].ok && (ev_data[c][i].st == 0))
	  Tp = (n < xnum_ps) ? hecke_matrices_weight(genus, W, S, q) : hecke_matrices_weight_square(genus, W, S, q);
    // whether q divides the discriminant exactly once (q^2 fits in a ulong, since q is an int)
    once = (fmpz_fdiv_ui(genus->disc, q) == 0) && (fmpz_fdiv_ui(genus->disc, (ulong)q * q) != 0);
    for (c = 0; c < genus->num_conductors; c++) {
      // D->hecke[k] is an array of fmpq_mat_t, i.e. of one-element arrays of fmpq_mat_struct
      const fmpq_mat_struct* Tc = (Tp == NULL) ? ((!by_eval) ? &(D->hecke[k][c][0]) : NULL) : Tp + c;
      for (i = 0; i < evs[c]->num; i++) {
	if (by_eval && (ev_data[c][i].st > 0))
	  has_ev[c][i * nq + n] = (ev_data[c][i].st == 1);
	else
	  has_ev[c][i * nq + n] = weight_form_eigenvalue(evals[c][i * nq + n], evs[c], i, Tc);
	if (has_ev[c][i * nq + n] && (n < xnum_ps) && once)
	  weight_bad_prime_eigenvalue(evals[c][i * nq + n], evs[c]->nfs[i], q, a, genus->conductors[c]);
	// (at the extra primes for the types, which are not printed, an old form is not typed)
	if (!has_ev[c][i * nq + n] && ((n < num_ps) || ((n >= xnum_ps) && (n - xnum_ps < num_ps2))))
	  fprintf(stderr, "Warning: T_{%d^%d} is not a scalar on the space of the old form %ld of conductor %ld\n",
		  q, (n < xnum_ps) ? 1 : 2, (long)i, (long)genus->conductors[c]);
      }
    }
    if (Tp != NULL) {
      for (c = 0; c < genus->num_conductors; c++)
	fmpq_mat_clear(Tp + c);
      free(Tp);
    }
  }
  decomposition_clear(D);
  free(data.ps);

  // The types, as for the weight 0: from the factors at the primes q of xps2 which do not divide the
  // discriminant (and are in xps), LIFT_TYPE_MIN_PRIMES of them at least. A form is a lift when it has the
  // same type at all of them, and G otherwise.
  if (set_types)
    for (c = 0; c < genus->num_conductors; c++)
      for (i = 0; i < evs[c]->num; i++) {
	if (evs[c]->lift_type[i] == O)
	  continue;
	for (n = 0; (n < xnum_ps2) && (evs[c]->lift_type[i] != G); n++) {
	  q = xps2[n];
	  if (fmpz_fdiv_ui(genus->disc, q) == 0)
	    continue;
	  for (m = 0; (m < xnum_ps) && (xps[m] != q); m++);
	  assert(m < xnum_ps);
	  assert(has_ev[c][i * nq + m] && has_ev[c][i * nq + xnum_ps + n]);
	  t = lpoly_type(evals[c][i * nq + m], evals[c][i * nq + xnum_ps + n], q, a, b, evs[c]->nfs[i]);
	  evs[c]->lift_type[i] = (evs[c]->lift_type[i] == U) ? t : ((evs[c]->lift_type[i] == t) ? t : G);
	}
      }

  printf("{");
  for (c = 0; c < genus->num_conductors; c++) {
    printf("%ld : [", (long)genus->conductors[c]);
    for (i = 0; i < evs[c]->num; i++) {
      printf("{ 'aut_rep_type' : '%c', ", "GYPFOU"[evs[c]->lift_type[i]]);
      printf("'field_poly' : '");
      fmpq_poly_print(evs[c]->nfs[i]->pol);
      printf("', 'lambda_p' : [");
      // only the eigenvalues at the primes of ps and ps2, which come first in xps and xps2
      for (n = 0; n < num_ps + num_ps2; n++) {
	k = (n < num_ps) ? n : xnum_ps + (n - num_ps);
	if (n == num_ps)
	  printf("], 'lambda_p_square' : [");
	if (has_ev[c][i * nq + k]) {
	  printf("'");
	  nf_elem_print_pretty(evals[c][i * nq + k], evs[c]->nfs[i], "a");
	  printf("',");
	}
	else
	  printf("'NULL',");
      }
      if (num_ps2 == 0)
	printf("], 'lambda_p_square' : [");
      printf("], }, ");
    }
    printf("], ");
  }
  printf("}");

  for (c = 0; c < genus->num_conductors; c++) {
    for (i = 0; i < evs[c]->num; i++)
      weight_eval_clear(&(ev_data[c][i]), evs[c]->nfs[i]);
    free(ev_data[c]);
    for (jj = 0; jj < h; jj++)
      weight_bpoly_clear(&(bpoly[c][jj]));
    free(bpoly[c]);
  }
  free(ev_data);
  free(bpoly);
  free(aut_order);
  for (c = 0; c < genus->num_conductors; c++) {
    for (i = 0; i < evs[c]->num * nq; i++)
      nf_elem_clear(evals[c][i], evs[c]->nfs[i / nq]);
    free(evals[c]);
    free(has_ev[c]);
    eigenvalues_clear(evs[c]);
  }
  free(evals);
  free(has_ev);
  free(evs);
  weight_space_clear(S);
  weight_rep_clear(W);
  free(xps);
  free(xps2);
  return 0;
}
