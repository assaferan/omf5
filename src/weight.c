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
					double* total_time, int* num_isom, int gen_idx)
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
  int lc, i;
  int num_isom;
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
  printf("theta_time = %f, isom_time = %f, total_time = %f, num_isom = %d / %d \n",
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
#include <flint/ulong_extras.h>
#include <flint/fmpq_mat.h>
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
  int num_isom, line, vec_cmp;
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

static fmpq_mat_struct* hecke_matrices_weight_impl(const genus_t genus, const weight_rep_t W,
						   const weight_space_t S, int p, int k, bool use_orbits)
{
  fmpq_mat_struct* T;
  fmpq_mat_t** acc;
  fmpq_mat_t R, M;
  nbr_data_t nbr_man;
  fmpz_mat_t nbr_fmpz, nbr_isom;
  square_matrix_t nbr;
  isometry_t s_nbr, hash_isom, s_inv;
  fmpz_t pa;
  slong c, i, j, h, nb, b;
  int num_isom;
  double theta_time, isom_time;
  W64 bits;
  bool needed;

  h = genus->genus_reps->num_stored;
  nb = (slong)1 << genus->spinor->num_primes;
  T = (fmpq_mat_struct*)malloc(genus->num_conductors * sizeof(fmpq_mat_struct));
  for (c = 0; c < genus->num_conductors; c++)
    fmpq_mat_init(T + c, S->dims[c], S->dims[c]);

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
    for (c = 0; c < genus->num_conductors; c++)
      needed = needed || (fmpq_mat_nrows(S->H[c] + j) > 0);
    if (!needed)
      continue;
    for (i = 0; i < h; i++)
      for (b = 0; b < nb; b++)
	fmpq_mat_zero(acc[i][b]);

    // the orbits are taken for T_p, also at the primes dividing the discriminant, as for the Hecke
    // matrices of the weight 0 (hecke_col_all_conductors). For T_{p^2} every neighbor of nbr_data is
    // tested for isometry
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
  return T;
}

fmpq_mat_struct* hecke_matrices_weight(const genus_t genus, const weight_rep_t W, const weight_space_t S, int p)
{
  return hecke_matrices_weight_impl(genus, W, S, p, 1, true);
}

fmpq_mat_struct* hecke_matrices_weight_by_neighbor(const genus_t genus, const weight_rep_t W, const weight_space_t S, int p)
{
  return hecke_matrices_weight_impl(genus, W, S, p, 1, false);
}

fmpq_mat_struct* hecke_matrices_weight_square(const genus_t genus, const weight_rep_t W, const weight_space_t S, int p)
{
  return hecke_matrices_weight_impl(genus, W, S, p, 2, true);
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

int print_eigenvalues_weight(const genus_t genus, int a, int b, const int* ps, int num_ps, const int* ps2,
			     int num_ps2, bool set_types)
{
  weight_rep_t W;
  weight_space_t S;
  weight_source_struct data;
  decomposition_t D;
  eigenvalues_t* evs;
  fmpq_mat_struct* Tp;
  // evals[c][i * nq + n]: the eigenvalue of the form i of conductor c of T_q for n < num_ps, and of
  // T_{q^2} otherwise, where q = ps[n] or ps2[n - num_ps], if has_ev[c][...]
  nf_elem_t** evals;
  bool** has_ev;
  slong c, i, k, nq;
  int n, m, q, num_type_ps;
  bool once;
  aut_type t;

  nq = num_ps + num_ps2;
  if (genus->genus_reps->num_stored == 0) {
    printf("{");
    for (c = 0; c < genus->num_conductors; c++)
      printf("%ld : [], ", (long)genus->conductors[c]);
    printf("}");
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
    if (n < num_ps) {
      q = ps[n];
      for (k = 0; (k < data.num_ps) && (data.ps[k] != q); k++);
      Tp = (k < data.num_ps) ? NULL : hecke_matrices_weight(genus, W, S, q);
    }
    else {
      // the decomposition has only matrices of T_p, so k is not used
      k = 0;
      q = ps2[n - num_ps];
      Tp = hecke_matrices_weight_square(genus, W, S, q);
    }
    // whether q divides the discriminant exactly once (q^2 fits in a ulong, since q is an int)
    once = (fmpz_fdiv_ui(genus->disc, q) == 0) && (fmpz_fdiv_ui(genus->disc, (ulong)q * q) != 0);
    for (c = 0; c < genus->num_conductors; c++) {
      // D->hecke[k] is an array of fmpq_mat_t, i.e. of one-element arrays of fmpq_mat_struct
      const fmpq_mat_struct* Tc = (Tp == NULL) ? &(D->hecke[k][c][0]) : Tp + c;
      for (i = 0; i < evs[c]->num; i++) {
	has_ev[c][i * nq + n] = weight_form_eigenvalue(evals[c][i * nq + n], evs[c], i, Tc);
	if (has_ev[c][i * nq + n] && (n < num_ps) && once)
	  weight_bad_prime_eigenvalue(evals[c][i * nq + n], evs[c]->nfs[i], q, a, genus->conductors[c]);
	if (!has_ev[c][i * nq + n])
	  fprintf(stderr, "Warning: T_{%d^%d} is not a scalar on the space of the old form %ld of conductor %ld\n",
		  q, (n < num_ps) ? 1 : 2, (long)i, (long)genus->conductors[c]);
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

  // The types, as for the weight 0: from the factors at the primes q of ps2 which do not divide the
  // discriminant and are in ps, at two of them at least. A form is a lift when it has the same type at
  // all of them, and G otherwise.
  num_type_ps = 0;
  if (set_types)
    for (n = 0; n < num_ps2; n++)
      if (fmpz_fdiv_ui(genus->disc, ps2[n]) != 0)
	num_type_ps++;
  if (set_types && (num_type_ps < 2))
    fprintf(stderr, "Warning: with a weight the types are determined at two primes p with p^2 <= L at least, which do not divide the discriminant: they are printed as U, and as O for the old forms. Use a larger L to determine them.\n");
  if (num_type_ps >= 2)
    for (c = 0; c < genus->num_conductors; c++)
      for (i = 0; i < evs[c]->num; i++) {
	if (evs[c]->lift_type[i] == O)
	  continue;
	for (n = 0; (n < num_ps2) && (evs[c]->lift_type[i] != G); n++) {
	  q = ps2[n];
	  if (fmpz_fdiv_ui(genus->disc, q) == 0)
	    continue;
	  for (m = 0; (m < num_ps) && (ps[m] != q); m++);
	  assert(m < num_ps);
	  assert(has_ev[c][i * nq + m] && has_ev[c][i * nq + num_ps + n]);
	  t = lpoly_type(evals[c][i * nq + m], evals[c][i * nq + num_ps + n], q, a, b, evs[c]->nfs[i]);
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
      for (n = 0; n < nq; n++) {
	if (n == num_ps)
	  printf("], 'lambda_p_square' : [");
	if (has_ev[c][i * nq + n]) {
	  printf("'");
	  nf_elem_print_pretty(evals[c][i * nq + n], evs[c]->nfs[i], "a");
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
  return 0;
}
