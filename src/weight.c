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

/* Algebraic modular forms of weight (a, 0): see weight.h. */

#include <flint/fmpz_mat.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_mpoly.h>

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
  fmpz_mat_t X;
  slong nullity, i, j, n;

  n = fmpz_mat_ncols(A);
  fmpz_mat_init(X, n, n);
  nullity = fmpz_mat_nullspace(X, A);
  fmpq_mat_init(B, nullity, n);
  for (i = 0; i < nullity; i++)
    for (j = 0; j < n; j++)
      fmpq_set_fmpz(fmpq_mat_entry(B, i, j), fmpz_mat_entry(X, j, i));
  fmpq_mat_rref(B, B);
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

void weight_rep_init(weight_rep_t W, int a, const square_matrix_t q0)
{
  slong i, j, k, r, nlow;
  ulong* low;
  ulong e[QF_RANK];
  fmpz_mat_t q, qinv, lap;
  fmpz_t den, c;
  fmpq_mat_t B;
  ulong* mon;

  W->a = a;
  W->nmons = weight_monomials(NULL, a);
  W->exps = (ulong*)malloc(QF_RANK * W->nmons * sizeof(ulong));
  weight_monomials(W->exps, a);

  if (a < 2) {
    // every polynomial of degree at most 1 is harmonic
    W->dim = W->nmons;
    fmpq_mat_init(W->basis, W->dim, W->nmons);
    fmpq_mat_one(W->basis);
    W->pivots = (slong*)malloc(W->dim * sizeof(slong));
    for (i = 0; i < W->dim; i++)
      W->pivots[i] = i;
  }
  else {
    // the laplacian sum_{i,j} qinv_ij d_i d_j, from degree a to degree a - 2, with qinv the adjugate
    // of q0, a multiple of its inverse, which has the same kernel
    nlow = weight_monomials(NULL, a - 2);
    low = (ulong*)malloc(QF_RANK * nlow * sizeof(ulong));
    weight_monomials(low, a - 2);
    fmpz_mat_init_set_square_matrix(q, q0);
    fmpz_mat_init(qinv, QF_RANK, QF_RANK);
    fmpz_init(den);
    fmpz_init(c);
    fmpz_mat_inv(qinv, den, q);
    fmpz_mat_init(lap, nlow, W->nmons);
    for (k = 0; k < W->nmons; k++) {
      mon = W->exps + QF_RANK*k;
      for (i = 0; i < QF_RANK; i++)
	for (j = i; j < QF_RANK; j++) {
	  for (r = 0; r < QF_RANK; r++)
	    e[r] = mon[r];
	  if (i == j) {
	    if (mon[i] < 2) continue;
	    e[i] -= 2;
	    fmpz_mul_ui(c, fmpz_mat_entry(qinv, i, i), mon[i] * (mon[i] - 1));
	  }
	  else {
	    if ((mon[i] == 0) || (mon[j] == 0)) continue;
	    e[i]--;
	    e[j]--;
	    fmpz_mul_ui(c, fmpz_mat_entry(qinv, i, j), 2 * mon[i] * mon[j]);
	  }
	  r = weight_monomial_index(low, nlow, e);
	  assert(r >= 0);
	  fmpz_add(fmpz_mat_entry(lap, r, k), fmpz_mat_entry(lap, r, k), c);
	}
    }
    W->dim = weight_nullspace_rref(B, &(W->pivots), lap);
    fmpq_mat_init_set(W->basis, B);
    fmpq_mat_clear(B);
    fmpz_mat_clear(lap);
    fmpz_mat_clear(q);
    fmpz_mat_clear(qinv);
    fmpz_clear(den);
    fmpz_clear(c);
    free(low);
  }

  // the basis as polynomials
  fmpq_mpoly_ctx_init(W->ctx, QF_RANK, ORD_LEX);
  W->polys = (fmpq_mpoly_struct*)malloc(W->dim * sizeof(fmpq_mpoly_struct));
  for (i = 0; i < W->dim; i++) {
    fmpq_mpoly_init(W->polys + i, W->ctx);
    for (k = 0; k < W->nmons; k++)
      if (!fmpq_is_zero(fmpq_mat_entry(W->basis, i, k)))
	fmpq_mpoly_set_coeff_fmpq_ui(W->polys + i, fmpq_mat_entry(W->basis, i, k), W->exps + QF_RANK*k, W->ctx);
  }

  return;
}

void weight_rep_clear(weight_rep_t W)
{
  slong i;

  for (i = 0; i < W->dim; i++)
    fmpq_mpoly_clear(W->polys + i, W->ctx);
  free(W->polys);
  fmpq_mpoly_ctx_clear(W->ctx);
  fmpq_mat_clear(W->basis);
  free(W->pivots);
  free(W->exps);
  return;
}

void weight_rep_action(fmpq_mat_t R, const weight_rep_t W, const isometry_t s)
{
  fmpq_mpoly_struct lin[QF_RANK];
  fmpq_mpoly_struct* lin_ptr[QF_RANK];
  fmpq_mpoly_t img, var;
  fmpq_t c;
  fmpz_mat_t S;
  fmpz_t det;
  slong i, j, k, m;
  int sign;

  // the linear forms x_i -> sum_j h_ij x_j, for h = s^(-1) = s_inv / inv_denom
  fmpq_init(c);
  fmpq_mpoly_init(var, W->ctx);
  for (i = 0; i < QF_RANK; i++) {
    fmpq_mpoly_init(lin + i, W->ctx);
    for (j = 0; j < QF_RANK; j++) {
      fmpq_set_si(c, s->s_inv[i][j], s->inv_denom);
      fmpq_mpoly_gen(var, j, W->ctx);
      fmpq_mpoly_scalar_mul_fmpq(var, var, c, W->ctx);
      fmpq_mpoly_add(lin + i, lin + i, var, W->ctx);
    }
    lin_ptr[i] = lin + i;
  }

  // det(s) = det(S) / denom^QF_RANK, which is 1 or -1
  fmpz_mat_init_set_square_matrix(S, s->s);
  fmpz_init(det);
  fmpz_mat_det(det, S);
  sign = ((fmpz_sgn(det) < 0) != (s->denom < 0)) ? -1 : 1;
  if ((W->a % 2) == 0)
    sign = 1;

  fmpq_mpoly_init(img, W->ctx);
  for (k = 0; k < W->dim; k++) {
    fmpq_mpoly_compose_fmpq_mpoly(img, W->polys + k, lin_ptr, W->ctx, W->ctx);
    for (m = 0; m < W->dim; m++) {
      fmpq_mpoly_get_coeff_fmpq_ui(c, img, W->exps + QF_RANK*(W->pivots[m]), W->ctx);
      if (sign < 0)
	fmpq_neg(c, c);
      fmpq_set(fmpq_mat_entry(R, m, k), c);
    }
  }

  fmpq_mpoly_clear(img, W->ctx);
  for (i = 0; i < QF_RANK; i++)
    fmpq_mpoly_clear(lin + i, W->ctx);
  fmpq_mpoly_clear(var, W->ctx);
  fmpz_mat_clear(S);
  fmpz_clear(det);
  fmpq_clear(c);
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
  slong k, m, r, dW;
  fmpq_mat_t b, u, v;
  const fmpq_mat_struct* Hi = S->H[c] + i;
  const fmpq_mat_struct* Hj = S->H[c] + j;

  dW = fmpq_mat_ncols(Hi);
  fmpq_mat_init(b, dW, 1);
  fmpq_mat_init(u, dW, 1);
  fmpq_mat_init(v, dW, 1);
  for (k = 0; k < fmpq_mat_nrows(Hi); k++) {
    for (r = 0; r < dW; r++)
      fmpq_set(fmpq_mat_entry(b, r, 0), fmpq_mat_entry(Hi, k, r));
    fmpq_mat_mul(u, M, b);
    // its coordinates in H[c][j] are its entries at the pivots; check that it lies there
    fmpq_mat_zero(v);
    for (m = 0; m < fmpq_mat_nrows(Hj); m++) {
      fmpq_add(fmpq_mat_entry(T, S->offsets[c][j] + m, S->offsets[c][i] + k),
	       fmpq_mat_entry(T, S->offsets[c][j] + m, S->offsets[c][i] + k),
	       fmpq_mat_entry(u, S->pivots[c][j][m], 0));
      for (r = 0; r < dW; r++)
	fmpq_addmul(fmpq_mat_entry(v, r, 0), fmpq_mat_entry(u, S->pivots[c][j][m], 0), fmpq_mat_entry(Hj, m, r));
    }
    if (!fmpq_mat_equal(u, v)) {
      fprintf(stderr, "Error: at the conductor index %ld, the image of a form under a Hecke operator is not invariant under the automorphisms of the lattice\n", (long)c);
      exit(EXIT_FAILURE);
    }
  }
  fmpq_mat_clear(b);
  fmpq_mat_clear(u);
  fmpq_mat_clear(v);
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
						   const weight_space_t S, int p, bool use_orbits)
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
  fmpz_pow_ui(pa, pa, W->a);
  fmpq_mat_init(R, W->dim, W->dim);
  fmpq_mat_init(M, W->dim, W->dim);
  fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
  fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
  theta_time = isom_time = 0;
  num_isom = 0;

  // acc[i][b], the sum of p^a R(s) over the neighbors of L_j isometric to L_i by s with spinor bits b
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

    if (use_orbits)
      weight_sums_orbits(acc, genus, W, j, p, pa);
    else {
    nbr_data_init(nbr_man, genus->genus_reps->keys[j], p, 1);
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
  return hecke_matrices_weight_impl(genus, W, S, p, true);
}

fmpq_mat_struct* hecke_matrices_weight_by_neighbor(const genus_t genus, const weight_rep_t W, const weight_space_t S, int p)
{
  return hecke_matrices_weight_impl(genus, W, S, p, false);
}

int print_hecke_matrices_weight(const genus_t genus, int a, int p)
{
  weight_rep_t W;
  weight_space_t S;
  fmpq_mat_struct* T;
  slong c, r, col;
  char* str;

  if (fmpz_fdiv_ui(genus->disc, p) == 0) {
    fprintf(stderr, "Error: with -weight, -p has to be a prime which does not divide the discriminant\n");
    return 1;
  }
  // a genus file may give no lattice: there are no forms
  if (genus->genus_reps->num_stored == 0) {
    printf("{");
    for (c = 0; c < genus->num_conductors; c++)
      printf("%ld : []%s", (long)genus->conductors[c], (c + 1 < genus->num_conductors) ? "," : "");
    printf("}");
    return 0;
  }
  weight_rep_init(W, a, genus->genus_reps->keys[0]);
  weight_space_init(S, genus, W);
  T = hecke_matrices_weight(genus, W, S, p);

  fprintf(stderr, "With the weight %d, the corresponding dimensions are: ", a);
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

/* Eigenvalues of the forms of weight (a, 0). */

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

int print_eigenvalues_weight(const genus_t genus, int a, const int* ps, int num_ps)
{
  weight_rep_t W;
  weight_space_t S;
  weight_source_struct data;
  decomposition_t D;
  eigenvalues_t* evs;
  fmpq_mat_struct* Tp;
  // evals[c][i * num_ps + n]: the eigenvalue of the form i of conductor c at ps[n], if has_ev[c][...]
  nf_elem_t** evals;
  bool** has_ev;
  slong c, i, k;
  int n;

  for (n = 0; n < num_ps; n++)
    if (fmpz_fdiv_ui(genus->disc, ps[n]) == 0) {
      fprintf(stderr, "Error: with -weight, the primes have to be prime to the discriminant\n");
      return 1;
    }
  if (genus->genus_reps->num_stored == 0) {
    printf("{");
    for (c = 0; c < genus->num_conductors; c++)
      printf("%ld : [], ", (long)genus->conductors[c]);
    printf("}");
    return 0;
  }

  weight_rep_init(W, a, genus->genus_reps->keys[0]);
  weight_space_init(S, genus, W);
  fprintf(stderr, "With the weight %d, the corresponding dimensions are: ", a);
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
    // the types are not determined for a weight
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
    k = (evs[c]->num * num_ps > 0) ? evs[c]->num * num_ps : 1;
    evals[c] = (nf_elem_t*)malloc(k * sizeof(nf_elem_t));
    has_ev[c] = (bool*)malloc(k * sizeof(bool));
    for (i = 0; i < evs[c]->num * num_ps; i++)
      nf_elem_init(evals[c][i], evs[c]->nfs[i / num_ps]);
  }
  for (n = 0; n < num_ps; n++) {
    for (k = 0; (k < data.num_ps) && (data.ps[k] != ps[n]); k++);
    Tp = (k < data.num_ps) ? NULL : hecke_matrices_weight(genus, W, S, ps[n]);
    for (c = 0; c < genus->num_conductors; c++) {
      // D->hecke[k] is an array of fmpq_mat_t, i.e. of one-element arrays of fmpq_mat_struct
      const fmpq_mat_struct* Tc = (Tp == NULL) ? &(D->hecke[k][c][0]) : Tp + c;
      for (i = 0; i < evs[c]->num; i++) {
	has_ev[c][i * num_ps + n] = weight_form_eigenvalue(evals[c][i * num_ps + n], evs[c], i, Tc);
	if (!has_ev[c][i * num_ps + n])
	  fprintf(stderr, "Warning: T_%d is not a scalar on the space of the old form %ld of conductor %ld\n",
		  ps[n], (long)i, (long)genus->conductors[c]);
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

  printf("{");
  for (c = 0; c < genus->num_conductors; c++) {
    printf("%ld : [", (long)genus->conductors[c]);
    for (i = 0; i < evs[c]->num; i++) {
      printf("{ 'aut_rep_type' : '%c', ", (evs[c]->lift_type[i] == O) ? 'O' : 'U');
      printf("'field_poly' : '");
      fmpq_poly_print(evs[c]->nfs[i]->pol);
      printf("', 'lambda_p' : [");
      for (n = 0; n < num_ps; n++) {
	if (has_ev[c][i * num_ps + n]) {
	  printf("'");
	  nf_elem_print_pretty(evals[c][i * num_ps + n], evs[c]->nfs[i], "a");
	  printf("',");
	}
	else
	  printf("'NULL',");
      }
      printf("], 'lambda_p_square' : [], }, ");
    }
    printf("], ");
  }
  printf("}");

  for (c = 0; c < genus->num_conductors; c++) {
    for (i = 0; i < evs[c]->num * num_ps; i++)
      nf_elem_clear(evals[c][i], evs[c]->nfs[i / num_ps]);
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
