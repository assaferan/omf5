#include <assert.h>

#include "carat/matrix.h"

#include "flint/fq_nmod.h"
#include "flint/fq_nmod_mat.h"
#include "flint/fq_nmod_mpoly.h"

#include "fq_nmod_mat.h"
#include "fq_nmod_mpoly.h"
#include "fq_nmod_quad.h"
#include "matrix_tools.h"

#include "nbr_data.h"
#include "typedefs.h"

/*****************************************************************************
 * This file implements p^2-neighbors
 * Since we care a lot about performance of p-neighbors,
 * and slightly less about p^2-neighbors, we separate the implementations
 * This implementation is more amenable to generalization, but less efficient
 *****************************************************************************/

void nbr_data_init(nbr_data_t nbr_man, const square_matrix_t q, slong p_int, slong k)
{
  slong idx;
  fmpz_t p, tmp;
#ifdef DEBUG_LEVEL_FULL
  fq_nmod_t value;
#endif // DEBUG_LEVEL_FULL

  assert(k > 0); // no perestroika operators in odd dimension
  
  fmpz_init_set_si(p, p_int);
  assert(fmpz_is_prime(p));
  fmpz_mat_init_set_square_matrix(nbr_man->q, q);
  fmpz_init(nbr_man->disc);
  fmpz_mat_det(nbr_man->disc, nbr_man->q);
  // (half) discriminant = det/2 when QF_RANK = 5 
  fmpz_divexact_si(nbr_man->disc, nbr_man->disc, 2);

  fq_nmod_ctx_init(nbr_man->GF, p, 1, "1");
  fq_nmod_mat_init_set_fmpz_mat(nbr_man->b, nbr_man->q, nbr_man->GF);
  if (p_int == 2) {
    fmpz_init(tmp);
    for (idx = 0; idx < QF_RANK; idx++) {
      fmpz_divexact_si(tmp, fmpz_mat_entry(nbr_man->q, idx, idx), 2); 
      fq_nmod_set_fmpz(fq_nmod_mat_entry(nbr_man->b,idx,idx), tmp, nbr_man->GF);
    }
    fmpz_clear(tmp);
  }
  nmod_mat_init_set_fmpz_mat(nbr_man->quot_gram, nbr_man->q, p_int*p_int);

  fq_nmod_mat_init(nbr_man->vec, 1, QF_RANK, nbr_man->GF);
#ifdef DEBUG
  fq_nmod_quad_isotropic_vector(nbr_man->vec, nbr_man->b, nbr_man->GF, 0, true);
#else
  fq_nmod_quad_isotropic_vector(nbr_man->vec, nbr_man->b, nbr_man->GF, 0, false);
#endif // DEBUG

#ifdef DEBUG_LEVEL_FULL
  fq_nmod_init(value, nbr_man->GF);
  fq_nmod_quad_evaluate(value, nbr_man->b, nbr_man->vec, nbr_man->GF);
  printf("b = \n");
  fq_nmod_mat_print_pretty(nbr_man->b, nbr_man->GF);
  printf("\n");
  printf("vec = \n");
  fq_nmod_mat_print_pretty(nbr_man->vec, nbr_man->GF);
  printf("\n");
  printf("value = ");
  fq_nmod_print(value, nbr_man->GF);
  printf("\n");
  assert(fq_nmod_is_zero(value, nbr_man->GF));
#endif // DEBUG_LEVEL_FULL

  fq_nmod_mat_init(nbr_man->p_std_gram, QF_RANK, QF_RANK, nbr_man->GF);
  fq_nmod_mat_init(nbr_man->p_basis, QF_RANK, QF_RANK, nbr_man->GF);
  fq_nmod_mat_init(nbr_man->p_skew, k, k, nbr_man->GF);

  nbr_man->is_skew_init = true;

#ifdef DEBUG
  fq_nmod_quad_decompose(nbr_man->p_std_gram, nbr_man->p_basis, nbr_man->b, nbr_man->GF, true);
#else
  fq_nmod_quad_decompose(nbr_man->p_std_gram, nbr_man->p_basis, nbr_man->b, nbr_man->GF, false);  
#endif // DEBUG

  fq_nmod_mpoly_ctx_init(nbr_man->p_q_std_ctx, QF_RANK, ORD_DEGREVLEX, nbr_man->GF);
  fq_nmod_mpoly_init(nbr_man->p_q_std,nbr_man->p_q_std_ctx);
  fq_nmod_poly_set_fq_nmod_quad(nbr_man->p_q_std, nbr_man->p_std_gram, nbr_man->GF, nbr_man->p_q_std_ctx);

#ifdef DEBUG_LEVEL_FULL
  printf("Performed Witt Decomposition on\n");
  fq_nmod_mat_print_pretty(nbr_man->b,nbr_man->GF);
  printf("\n");
  printf("Resulting gram matrix is \n");
  fq_nmod_mat_print_pretty(nbr_man->p_std_gram,nbr_man->GF);
  printf("\n");
  printf("Resulting basis is \n");
  fq_nmod_mat_print_pretty(nbr_man->p_basis,nbr_man->GF);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  // Count the rows at the end of the matrix which are exactly zero.
  idx = QF_RANK;
  while ((idx >= 1) && fq_nmod_mat_is_zero_row(nbr_man->p_std_gram,idx-1,nbr_man->GF)) idx--;

  // The dimension of the radical.
  nbr_man->rad_dim = QF_RANK - idx;

  // Determine the dimension of the totally hyperbolic subspace.
  idx = 1;
  while ((idx <= QF_RANK - nbr_man->rad_dim) && fq_nmod_is_zero(fq_nmod_mat_entry(nbr_man->p_std_gram,idx-1,idx-1),nbr_man->GF)) idx++;

  // Dimension of the anistotropic subspace.
  nbr_man->aniso_dim = QF_RANK - nbr_man->rad_dim - idx + 1;

  // The number of hyperbolic planes in the Witt decomposition.
  nbr_man->witt_index = (idx - 1) / 2;

  pivot_data_init(nbr_man->pivots, QF_RANK - nbr_man->rad_dim, nbr_man->aniso_dim, k);
  nbr_man->pivots->pivot_ptr = 0;
  nbr_man->k = k;
  nbr_man->skew_dim = k*(k-1)/2;
  nbr_man->is_done = false;
  nbr_man->pivots->num_params = 0;
  nbr_man->is_iso_subspace_init = false;

  nbr_data_next_isotropic_subspace(nbr_man);

  nmod_mat_init(nbr_man->X, nbr_man->k, QF_RANK, p_int*p_int);
  nmod_mat_init(nbr_man->Z, nbr_man->k, QF_RANK, p_int*p_int);
  nmod_mat_init(nbr_man->U, QF_RANK - 2*nbr_man->k, QF_RANK, p_int*p_int);
  nbr_data_lift_subspace(nbr_man);

  nmod_mat_init(nbr_man->X_skew, nbr_man->k, QF_RANK, p_int*p_int);
  if (!(nbr_man->is_done))
    nmod_mat_set(nbr_man->X_skew, nbr_man->X);
  
#ifdef DEBUG_LEVEL_FULL
  fq_nmod_clear(value, nbr_man->GF);
#endif // DEBUG_LEVEL_FULL
  fmpz_clear(p);
  
  return;
}

void nbr_data_clear(nbr_data_t nbr_man)
{
  nmod_mat_clear(nbr_man->X_skew);
  nmod_mat_clear(nbr_man->X);
  nmod_mat_clear(nbr_man->Z);
  nmod_mat_clear(nbr_man->U);
  if (nbr_man->is_iso_subspace_init)
    fq_nmod_mat_clear(nbr_man->iso_subspace, nbr_man->GF);
  if (nbr_man->pivots->is_params_init)
    pivot_data_params_clear(nbr_man->pivots);
  pivot_data_clear(nbr_man->pivots);
  fq_nmod_mpoly_clear(nbr_man->p_q_std,nbr_man->p_q_std_ctx);
  fq_nmod_mpoly_ctx_clear(nbr_man->p_q_std_ctx);
  fq_nmod_mat_clear(nbr_man->p_std_gram, nbr_man->GF);
  fq_nmod_mat_clear(nbr_man->p_basis, nbr_man->GF);
  fq_nmod_mat_clear(nbr_man->p_skew, nbr_man->GF);
  fq_nmod_mat_clear(nbr_man->vec, nbr_man->GF);
  nmod_mat_clear(nbr_man->quot_gram);
  fq_nmod_mat_clear(nbr_man->b,nbr_man->GF);
  fq_nmod_ctx_clear(nbr_man->GF);
  fmpz_clear(nbr_man->disc);
  fmpz_mat_clear(nbr_man->q);
}

void nbr_data_next_isotropic_subspace(nbr_data_t nbr_man)
{
  // fq_nmod_mat_t eval_list;
  fq_nmod_struct** eval_list;
  slong i, j, pos, n;
  bool all_zero;
  fq_nmod_t one;
#ifdef DEBUG_LEVEL_FULL
  const char* var_names[10] = {"x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10"};
#endif // DEBUG_LEVEL_FULL

  n = fmpz_mat_nrows(nbr_man->q);
  assert(n == fmpz_mat_ncols(nbr_man->q));
  
  if (nbr_man->pivots->num_params == 0) {
    // Move to the next pivot.
    nbr_man->pivots->pivot_ptr++;

    // If we've exceeded the list of pivots, we're done.
    if (nbr_man->pivots->pivot_ptr > nbr_man->pivots->total_len) {
      // Reset the pivot pointer so that we can loop through
      //  the isotropic subspaces again if needed.
      nbr_man->pivots->pivot_ptr = 0;
      if (nbr_man->is_iso_subspace_init) {
	fq_nmod_mat_clear(nbr_man->iso_subspace, nbr_man->GF);
	nbr_man->is_iso_subspace_init = false;
      }
      nbr_man->is_done = true;
      return;
    }
    // Initialize the new pivot.
    if (nbr_man->pivots->is_params_init)
      pivot_data_params_clear(nbr_man->pivots);
    nbr_data_params_init(nbr_man->pivots, nbr_man);
  }

  // The list of evaluation values.
  // fq_nmod_mat_init(eval_list, 1, n * nbr_man->k, nbr_man->GF);
  // fq_nmod_mat_zero(eval_list, nbr_man->GF);
  eval_list = (fq_nmod_struct**)malloc(n*nbr_man->k*sizeof(fq_nmod_struct*));
  for (i = 0; i < n * nbr_man->k; i++) {
    eval_list[i] = (fq_nmod_struct*)malloc(sizeof(fq_nmod_struct));
    fq_nmod_init(eval_list[i], nbr_man->GF);
    fq_nmod_zero(eval_list[i], nbr_man->GF);
  }

  assert(nbr_man->pivots->is_params_init);
  
  // Produce the isotropic subspace corresponding to the current
  //  parameters.

  for (i = 0; i < nbr_man->pivots->num_params; i++)
    // fq_nmod_set(fq_nmod_mat_entry(eval_list,0,i), nbr_man->pivots->params[i], nbr_man->GF);
    fq_nmod_set(eval_list[nbr_man->pivots->free_vars[i]], nbr_man->pivots->params[i], nbr_man->GF);

  if (nbr_man->is_iso_subspace_init)
    fq_nmod_mat_clear(nbr_man->iso_subspace, nbr_man->GF);
  // The basis for the current isotropic subspace.
  fq_nmod_mat_init(nbr_man->iso_subspace, nbr_man->k, n, nbr_man->GF);
  nbr_man->is_iso_subspace_init = true;
  for (i = 0; i < nbr_man->k; i++) {
    for (j = 0; j < n; j++) {
      fq_nmod_mpoly_evaluate_all_fq_nmod(fq_nmod_mat_entry(nbr_man->iso_subspace,i,j),
					 fq_nmod_mpoly_mat_entry(nbr_man->pivots->p_isotropic_param,i,j),
					 eval_list, nbr_man->pivots->R);
    }
  }

#ifdef DEBUG_LEVEL_FULL
  printf("params = ");
  for (i = 0; i < nbr_man->pivots->num_params; i++)
    fq_nmod_print_pretty(nbr_man->pivots->params[i], nbr_man->GF);
  printf("\n");
  printf("free vars = ");
  for (i = 0; i < nbr_man->pivots->num_free_vars; i++)
    printf("%ld ", nbr_man->pivots->free_vars[i]);
  printf("\n");
  printf("isotropic_param = ");
  fq_nmod_mpoly_mat_print(nbr_man->pivots->p_isotropic_param, var_names, nbr_man->pivots->R);
  printf("\n");
  printf("isotropic_subspace = ");
  fq_nmod_mat_print_pretty(nbr_man->iso_subspace, nbr_man->GF);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  if (nbr_man->pivots->num_free_vars != 0) {
    fq_nmod_init(one, nbr_man->GF);
    fq_nmod_one(one, nbr_man->GF);
    // The current position in the parameterization.
    pos = 0;
    // Terminate loop once we found the next new subspace, or we
    //  hit the end of the list.
    do {
      // Increment position.
      pos++;
      // Manually move to the next element.
      assert(pos <= nbr_man->pivots->num_params);
      fq_nmod_add(nbr_man->pivots->params[pos-1], nbr_man->pivots->params[pos-1], one, nbr_man->GF);
    } while((pos != nbr_man->pivots->num_free_vars) && (fq_nmod_is_zero(nbr_man->pivots->params[pos-1], nbr_man->GF)) );
    fq_nmod_clear(one, nbr_man->GF);
  }

  // If we've hit the end of the list, indicate we need to move on to the
  //  next pivot.
  all_zero = true;
  for (i = 0; i < nbr_man->pivots->num_params; i++)
    if (!fq_nmod_is_zero(nbr_man->pivots->params[i], nbr_man->GF)) {
      all_zero = false;
      break;
    }

  if (all_zero) {
    pivot_data_params_clear(nbr_man->pivots);
  }
  
  // fq_nmod_mat_clear(eval_list, nbr_man->GF);
  for (i = 0; i < n * nbr_man->k; i++) {
    fq_nmod_clear(eval_list[i], nbr_man->GF);
    free(eval_list[i]);
  }
  free(eval_list);
  
  return;
}

void nbr_data_params_init(pivot_data_t pivots, const nbr_data_t nbr_man)
{
  slong* pivot;
  slong rank, row, col, remove_idx, remove_size, rows, i, j, c, data_idx, vec_idx, n;
  // Keep a list of non-free variables from which we will determine the
  //  free variables when we are done.
  slong* remove;
  fq_nmod_mpoly_mat_t data;
  // fq_nmod_mpoly_mat_t eval_list;
  fq_nmod_mpoly_struct** eval_list;
  // fq_nmod_mpoly_mat_t vec;
  fq_nmod_mpoly_struct** vec;
  
  fq_nmod_mat_t mat, trans, l;
  fq_nmod_mpoly_t f, g;
  bool appears;

#ifdef DEBUG_LEVEL_FULL
  // TODO - find a better way to do that
  const char* var_names[10] = {"x1","x2","x3","x4","x5","x6","x7","x8","x9","x10"};
#endif // DEBUG_LEVEL_FULL
  
  assert(pivots->pivot_ptr >= 1);

  n = fmpz_mat_nrows(nbr_man->q);
  assert(n == fmpz_mat_ncols(nbr_man->q));
  
  pivot = pivots->pivots[pivots->pivot_ptr-1];
  rank = nbr_man->k * n;

  // Initialize matrix that will determine parameterization.
  fq_nmod_mpoly_ctx_init(pivots->R, rank, ORD_DEGREVLEX, nbr_man->GF);
  /* fq_nmod_mpoly_mat_init(data, 1, rank, pivots->R); */
  /* for (i = 0; i < rank; i++) */
  /*   fq_nmod_mpoly_gen(fq_nmod_mpoly_mat_entry(data,0,i),i,pivots->R); */

  /* fq_nmod_mpoly_mat_init_set_fq_nmod_mpoly_vec(pivots->p_isotropic_param, data, nbr_man->k, n, pivots->R); */

  fq_nmod_mpoly_mat_init(pivots->p_isotropic_param, nbr_man->k, n, pivots->R);
  i = 0;
  for (row = 0; row < nbr_man->k; row++)
    for (col = 0; col < n; col++)
      fq_nmod_mpoly_gen(fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,row,col),i++,pivots->R);

  // Setup the columns corresponding to the pivots.
  remove_size = (nbr_man->k)*(nbr_man->k);
  for (row = 0; row < nbr_man->k; row++) {
    remove_size += pivot[row];
    if (pivot[row] >= nbr_man->witt_index)
      remove_size += nbr_man->aniso_dim;
  }

  // Determine the number of rows of the matrix that we'll echelonize.
  rows = nbr_man->k * (nbr_man->k + 1) / 2;

  remove_size += rows;
  remove_size += rank; // We fill in only some of these
  
  remove = (slong*)malloc(remove_size*sizeof(slong));
  remove_idx = 0;
  for (row = 0; row < nbr_man->k; row++)
    for (col = 0; col < nbr_man->k; col++) {
      if (row == col)
	fq_nmod_mpoly_one(fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param, row, pivot[col]), pivots->R);
      else
	fq_nmod_mpoly_zero(fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param, row, pivot[col]), pivots->R);
      remove[remove_idx++] = row*n+pivot[col];
    }

  // Clear the rows prior to the pivot positions (but not the radical).
  for (row = 0; row < nbr_man->k; row++)
    for (col = 0; col < pivot[row]; col++) {
      fq_nmod_mpoly_zero(fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param, row, col), pivots->R);
      remove[remove_idx++] = row*n+col;
    }

  // Check if one or more of the anisotropic coordinates need to be zero.
  for (row = 0; row < nbr_man->k; row++) {
    if (pivot[row] >= nbr_man->witt_index) {
      for (col = 0; col < nbr_man->aniso_dim; col++) {
	fq_nmod_mpoly_zero(fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param, row, n-1-nbr_man->rad_dim-col), pivots->R);
	remove[remove_idx++] = (row+1)*n - 1 - nbr_man->rad_dim-col;
      }
    }
  }

  // Here we will save the quadratic polynomials
  fq_nmod_mpoly_mat_init(data, 1, rows, pivots->R);

  // The matrix that we're going to echelonize.
  fq_nmod_mat_init(mat, rows, rank, nbr_man->GF);

  // fq_nmod_mpoly_mat_init(vec, 1, rows*n, pivots->R);
  vec = (fq_nmod_mpoly_struct**)malloc(n*sizeof(fq_nmod_mpoly_struct*));
  for (vec_idx = 0; vec_idx < n; vec_idx++) {
    vec[vec_idx] = (fq_nmod_mpoly_struct*)malloc(sizeof(fq_nmod_mpoly_struct));
    fq_nmod_mpoly_init(vec[vec_idx], pivots->R);
  }
  
  fq_nmod_mat_init(l, 1, rank, nbr_man->GF);
  fq_nmod_mpoly_init(f, pivots->R);

#ifdef DEBUG_LEVEL_FULL
  printf("isotropic_param = \n");
  fq_nmod_mpoly_mat_print(pivots->p_isotropic_param, var_names, pivots->R);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  data_idx = 0;
  row = 0;
  for (i = 0; i < nbr_man->k; i++)
    for (j = i; j < nbr_man->k; j++) {
      // The appropriate vector that we want to be isotropic.
      for (vec_idx = 0; vec_idx < n; vec_idx++) {
	fq_nmod_mpoly_set(vec[vec_idx], fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,i,vec_idx), pivots->R);
	if (i != j)
	  fq_nmod_mpoly_add(vec[vec_idx],
			    vec[vec_idx], fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,j,vec_idx), pivots->R);
      }
#ifdef DEBUG_LEVEL_FULL
      printf("Substituting vec = ");
      for (vec_idx = 0; vec_idx < n; vec_idx++) {
	fq_nmod_mpoly_print_pretty(vec[vec_idx], var_names, pivots->R);
	printf(" ");
      }
      printf(" in q_std = ");
      fq_nmod_mpoly_print_pretty(nbr_man->p_q_std, var_names, nbr_man->p_q_std_ctx);
#endif // DEBUG_LEVEL_FULL
      fq_nmod_mpoly_compose_fq_nmod_mpoly(f, nbr_man->p_q_std, vec, nbr_man->p_q_std_ctx, pivots->R);
#ifdef DEBUG_LEVEL_FULL
      printf(" yields f = ");
      fq_nmod_mpoly_print_pretty(f, var_names, pivots->R);
      printf("\n");
#endif // DEBUG_LEVEL_FULL
      // Degree 2 terms are inhomogenous.
      fq_nmod_mpoly_quadratic_part(fq_nmod_mpoly_mat_entry(data,0,data_idx), f, pivots->R);
#ifdef DEBUG_LEVEL_FULL
      printf("Quadratic part = ");
      fq_nmod_mpoly_print_pretty(fq_nmod_mpoly_mat_entry(data,0,data_idx), var_names, pivots->R);
      printf("\n");
#endif // DEBUG_LEVEL_FULL
      fq_nmod_mpoly_neg(fq_nmod_mpoly_mat_entry(data,0,data_idx), fq_nmod_mpoly_mat_entry(data,0,data_idx), pivots->R);
      data_idx++;
      // The other terms are linear
      // so fill in the matrix accordingly.
      fq_nmod_mpoly_linear_part(l, f, pivots->R);
      for (vec_idx = 0; vec_idx < rank; vec_idx++)
	fq_nmod_set(fq_nmod_mat_entry(mat, row, vec_idx), fq_nmod_mat_entry(l, 0, vec_idx), nbr_man->GF);
      // Move on to the next row.
      row++;
      // Silly question. Isn't row == i?
    }

#ifdef DEBUG_LEVEL_FULL
  printf("The matrix before echelon is mat = \n");
  fq_nmod_mat_print_pretty(mat, nbr_man->GF);
  printf("\n");
  printf("The last entry is the quadratic data = \n");
  for (i = 0; i < rows; i++) {
    fq_nmod_mpoly_print_pretty(fq_nmod_mpoly_mat_entry(data,0,i), var_names, pivots->R);
    printf(" ");
  }
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  // Compute the Echelon form of the matrix.
  fq_nmod_mat_init(trans, rows, rows, nbr_man->GF);
  // fq_nmod_mat_rref(mat, trans, nbr_man->GF);
  fq_nmod_mat_rref_trans(mat, trans, nbr_man->GF);

#ifdef DEBUG_LEVEL_FULL
  printf("The matrix after echelon is mat = \n");
  fq_nmod_mat_print_pretty(mat, nbr_man->GF);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  // The evaluation list for replacing variables with their dependence
  //  relations.

  // fq_nmod_mpoly_mat_init(eval_list, 1, rank, pivots->R);
  eval_list = (fq_nmod_mpoly_struct**)malloc(rank*sizeof(fq_nmod_mpoly_struct*));
  for (i = 0; i < rank; i++) {
    eval_list[i] = (fq_nmod_mpoly_struct*)malloc(sizeof(fq_nmod_mpoly_struct));
    fq_nmod_mpoly_init(eval_list[i], pivots->R);
    fq_nmod_mpoly_gen(eval_list[i],i,pivots->R);
  }

  fq_nmod_mpoly_init(g, pivots->R);
  
  for (i = 0; i < rows; i++) {
    // Find the pivot in the i-th row.
    c = 0;
    while ((c < rank) && (!fq_nmod_is_one(fq_nmod_mat_entry(mat,i,c), nbr_man->GF))) c++;
    // Add this pivot to the list of non-free variables.
    remove[remove_idx++] = c;

    // If the row is entirely zero, skip it.
    if (c >= rank) continue;

    // Build the multinomial for which x_c is dependent.
    fq_nmod_mpoly_zero(f, pivots->R);
    for (j = 0; j < rows; j++) {
      fq_nmod_mpoly_scalar_mul_fq_nmod(g, fq_nmod_mpoly_mat_entry(data,0,j), fq_nmod_mat_entry(trans,i,j), pivots->R);
      fq_nmod_mpoly_add(f, f, g, pivots->R);
    }
    for (j = 0; j < rank; j++) {
      if (j != c) {
	fq_nmod_mpoly_gen(g, j, pivots->R);
	fq_nmod_mpoly_scalar_mul_fq_nmod(g, g, fq_nmod_mat_entry(mat,i,j), pivots->R);
	fq_nmod_mpoly_sub(f, f, g, pivots->R);
      }
    }
    fq_nmod_mpoly_set(eval_list[c], f, pivots->R);
  }

  // The matrix whose rows parameterize all isotropic subspaces.
  for (row = 0; row < fq_nmod_mpoly_mat_nrows(pivots->p_isotropic_param); row++)
    for (col = 0; col < fq_nmod_mpoly_mat_ncols(pivots->p_isotropic_param); col++)
      fq_nmod_mpoly_compose_fq_nmod_mpoly(fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param, row, col),
					  fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param, row, col), eval_list, pivots->R, pivots->R);

#ifdef DEBUG
  // Verify that we didn't screw up somewhere along the line.
#ifdef DEBUG_LEVEL_FULL
  printf("testing parametrization\n");
  printf("p_isotropic_param = \n");
  fq_nmod_mpoly_mat_print(pivots->p_isotropic_param, var_names, pivots->R);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  for (i = 0; i < nbr_man->k; i++)
    for (j = i; j < nbr_man->k; j++) {
      for (vec_idx = 0; vec_idx < n; vec_idx++) {
	fq_nmod_mpoly_set(vec[vec_idx], fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,i,vec_idx), pivots->R);
	if (i != j)
	  fq_nmod_mpoly_add(vec[vec_idx],
			    vec[vec_idx], fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,j,vec_idx), pivots->R);
      }
  
#ifdef DEBUG_LEVEL_FULL
      printf("Substituting vec = ");
      for (vec_idx = 0; vec_idx < n; vec_idx++) {
	fq_nmod_mpoly_print_pretty(vec[vec_idx], var_names, pivots->R);
	printf(" ");
      }
      printf(" in q_std = ");
      fq_nmod_mpoly_print_pretty(nbr_man->p_q_std, var_names, nbr_man->p_q_std_ctx);
#endif // DEBUG_LEVEL_FULL
      
      fq_nmod_mpoly_compose_fq_nmod_mpoly(f, nbr_man->p_q_std, vec, nbr_man->p_q_std_ctx, pivots->R);
#ifdef DEBUG_LEVEL_FULL
      printf(" yields f = ");
      fq_nmod_mpoly_print_pretty(f, var_names, pivots->R);
      printf("\n");
#endif // DEBUG_LEVEL_FULL
      
      assert(fq_nmod_mpoly_is_zero(f,pivots->R));
    }
#endif // DEBUG

  // !! TODO - needs to enlarge remove!!!
  // Determine the free variables.
  for (i = 0; i < rank; i++) {
    appears = false;
    for (row = 0; row < fq_nmod_mpoly_mat_nrows(pivots->p_isotropic_param); row++)
      for (col = 0; col < fq_nmod_mpoly_mat_ncols(pivots->p_isotropic_param); col++)
	if (fq_nmod_mpoly_degree_si(fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,row,col),i,pivots->R) > 0) {
	  appears = true;
	  break;
	}
    if (!appears)
      remove[remove_idx++] = i;
  }

  if (pivots->is_params_init)
    free(pivots->free_vars);
  pivots->free_vars = (slong*)malloc(rank*sizeof(slong)); // not using all of them
  pivots->params = (fq_nmod_t*)malloc(rank*sizeof(fq_nmod_t));
  pivots->num_free_vars = 0;
  pivots->num_params = 0;
  for (i = 0; i < rank; i++) {
    appears = false;
    for (j = 0; j < remove_idx; j++) {
      if (remove[j] == i) {
	appears = true;
	break;
      }
    }
    if (!appears) {
      pivots->free_vars[pivots->num_free_vars++] = i;
      fq_nmod_init(pivots->params[pivots->num_params],nbr_man->GF);
      fq_nmod_zero(pivots->params[pivots->num_params++], nbr_man->GF);
    }
  }

  pivots->is_params_init = true;
  
  for (i = 0; i < rank; i++) {
    fq_nmod_mpoly_clear(eval_list[i], pivots->R);
    free(eval_list[i]);
  }
  
  free(eval_list);
  
  fq_nmod_mpoly_clear(g, pivots->R);

  fq_nmod_mat_clear(trans, nbr_man->GF);
  
  fq_nmod_mat_clear(l, nbr_man->GF);
  fq_nmod_mpoly_clear(f, pivots->R);

  for (i = 0; i < n; i++) {
    fq_nmod_mpoly_clear(vec[i], pivots->R);
    free(vec[i]);
  }
  free(vec);
  
  fq_nmod_mat_clear(mat,nbr_man->GF);
    
  free(remove);

  fq_nmod_mpoly_mat_clear(data, pivots->R);
  
  return;
}

void nmod_mat_gram(nmod_mat_t gram, const nmod_mat_t B, const nmod_mat_t q, slong pp)
{
  nmod_mat_t Bq, B_t;
  slong i,j,n;

  n = nmod_mat_nrows(B);
  assert(n == nmod_mat_ncols(B));
  
  nmod_mat_init(Bq, n, n, pp);
  nmod_mat_init(B_t, n, n, pp);

  nmod_mat_mul(Bq, B, q);
  nmod_mat_transpose(B_t, B);
  nmod_mat_mul(gram, Bq, B_t);

  // !! TODO - this isn't necessary only for debugging versus magma
  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(gram,i,j) %= pp;
  
  nmod_mat_clear(B_t);
  nmod_mat_clear(Bq);
  return;
}

void nmod_mat_gram_fmpz_mat(nmod_mat_t gram, const nmod_mat_t B, const fmpz_mat_t q)
{
  fmpz_mat_t B_fmpz, Bq, B_t, gram_fmpz;
  slong i,j,n;

  n = nmod_mat_nrows(B);
  assert(n == nmod_mat_ncols(B));

  fmpz_mat_init(B_fmpz, n, n);
  fmpz_mat_init(gram_fmpz, n, n);
  fmpz_mat_init(Bq, n, n);
  fmpz_mat_init(B_t, n, n);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      fmpz_set_si(fmpz_mat_entry(B_fmpz,i,j), nmod_mat_entry(B,i,j));

  fmpz_mat_mul(Bq, B_fmpz, q);
  fmpz_mat_transpose(B_t, B_fmpz);
  fmpz_mat_mul(gram_fmpz, Bq, B_t);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(gram, i, j) = fmpz_get_si(fmpz_mat_entry(gram_fmpz,i,j));
  
  fmpz_mat_clear(B_t);
  fmpz_mat_clear(Bq);
  fmpz_mat_clear(B_fmpz);
  fmpz_mat_clear(gram_fmpz);
  return;
}

void nbr_data_lift_subspace(nbr_data_t nbr_man)
{
  slong p;
  slong* pivots;
  slong* paired;
  slong j = 0;
  slong i, l, u_idx, n;
  slong h_dim, delta, a, half, gram2, scalar;
  fq_nmod_mat_t basis, basis_t, x, z, u;
  bool* excluded;
  nmod_mat_t B, X_new, Z_new;
  // fmpz_mat_t B;
  nmod_mat_t gram;

#ifdef DEBUG
  nmod_mat_t temp;
#endif // DEBUG

#ifdef DEBUG_LEVEL_FULL
  slong num_pivots;
#endif // DEBUG_LEVEL_FULL
  
  if (nbr_man->is_done) return;

  n = fmpz_mat_nrows(nbr_man->q);
  assert(n == fmpz_mat_ncols(nbr_man->q));
  
  p = fmpz_get_si(fq_nmod_ctx_prime(nbr_man->GF));

  assert(nbr_man->pivots->pivot_ptr >= 1);

  // Get the pivots for the bases of the isotropic subspaces.
  pivots = nbr_man->pivots->pivots[nbr_man->pivots->pivot_ptr-1];

#ifdef DEBUG_LEVEL_FULL
  num_pivots = nbr_man->pivots->pivot_lens[nbr_man->pivots->pivot_ptr-1];
  printf("before lifting, p_basis is \n");
  fq_nmod_mat_print_pretty(nbr_man->p_basis, nbr_man->GF);
  printf("\n");
  printf("iso_subspace is \n");
  fq_nmod_mat_print_pretty(nbr_man->iso_subspace, nbr_man->GF);
  printf("\n");
  printf("pivots = ");
  for (i = 0; i < num_pivots; i++)
    printf("%ld", pivots[i]);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_init_set(basis, nbr_man->p_basis, nbr_man->GF);
  // Set up the correct basis vectors.
  for (i = 0; i < nbr_man->k; i++)
    for (j = pivots[i] + 1; j < n; j++)

#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
      {
	printf("iso_subspace is \n");
	fq_nmod_mat_print_pretty(nbr_man->iso_subspace, nbr_man->GF);
	printf("\n");
	printf("basis is \n");
	fq_nmod_mat_print_pretty(basis, nbr_man->GF);
	printf("\n");
	printf("pivots = ");
	num_pivots = nbr_man->pivots->pivot_lens[nbr_man->pivots->pivot_ptr-1];
	for (l = 0; l < num_pivots; l++)
	  printf("%ld", pivots[l]);
	printf("\n");
	printf("j = %ld\n", j);
	fq_nmod_mat_add_col(basis, pivots[i], j, fq_nmod_mat_entry(nbr_man->iso_subspace,i,j), nbr_man->GF);
      }
  
#ifdef DEBUG_LEVEL_FULL
  printf("the correct basis vectors are \n");
  fq_nmod_mat_print_pretty(basis, nbr_man->GF);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_init(basis_t, fq_nmod_mat_ncols(basis, nbr_man->GF), fq_nmod_mat_nrows(basis, nbr_man->GF), nbr_man->GF);
#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_transpose(basis_t, basis, nbr_man->GF);
  
  // Extract our target isotropic subspace modulo p
#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_init(x, nbr_man->k, n, nbr_man->GF);
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      fq_nmod_set(fq_nmod_mat_entry(x,i,j), fq_nmod_mat_entry(basis_t, pivots[i],j), nbr_man->GF);

#ifdef DEBUG_LEVEL_FULL
  printf("x = ");
  fq_nmod_mat_print_pretty(x, nbr_man->GF);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  // Extract the hyperbolic complement modulo p.
  //  printf("initializing z and paired...\n");

#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_init(z, nbr_man->k, n, nbr_man->GF);
  paired = (slong*)malloc((nbr_man->k)*sizeof(slong));
  h_dim = 2 * nbr_man->witt_index;
  // printf("h_dim = %ld\n", h_dim);
  for (i = 0 ; i < nbr_man->k; i++) {
    // printf("paired[%ld] = %ld - 1 - pivots[%ld]\n", i, h_dim, nbr_man->k-1-i);
    paired[i] = h_dim - 1 - pivots[nbr_man->k-1-i];
  }
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++) {
      // printf("z(%ld,%ld) = basis_t(%ld,%ld)\n", i, j, paired[i],j);
      fq_nmod_set(fq_nmod_mat_entry(z,i,j), fq_nmod_mat_entry(basis_t, paired[i],j), nbr_man->GF);
    }

#ifdef DEBUG_LEVEL_FULL
  printf("z = ");
  fq_nmod_mat_print_pretty(z, nbr_man->GF);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  // Extract the remaining basis vectors.
  excluded = (bool*)malloc(n * sizeof(bool));
  for (i = 0; i < n; i++)
    excluded[i] = true;
  for (i = 0; i < nbr_man->k; i++) {
    excluded[pivots[i]] = false;
    excluded[paired[i]] = false;
  }

#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_init(u, n-2*nbr_man->k, n, nbr_man->GF);
  u_idx = 0;
  for (i = 0; i < n; i++)
    if (excluded[i]) {
       for (j = 0; j < n; j++)
	 fq_nmod_set(fq_nmod_mat_entry(u,u_idx,j), fq_nmod_mat_entry(basis_t, i,j), nbr_man->GF);
       u_idx++;
    }

#ifdef DEBUG_LEVEL_FULL
  printf("u = ");
  fq_nmod_mat_print_pretty(u, nbr_man->GF);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  // Convert to coordinates modulo p^2.
  nmod_mat_init(B, n, n, p*p);
  
  // Build the coordinate matrix.
  // !! TODO - the mod p is not necessary, good for debugging
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++) {
      nmod_mat_entry(B,i,j) = nmod_mat_entry(nbr_man->X,i,j) = nmod_poly_get_coeff_ui(fq_nmod_mat_entry(x,i,j),0) % p;
    }

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
       nmod_mat_entry(B,nbr_man->k+i,j) =
	 nmod_mat_entry(nbr_man->Z,i,j) = nmod_poly_get_coeff_ui(fq_nmod_mat_entry(z,i,j),0) % p;

  for (i = 0; i < n-2*nbr_man->k; i++)
    for (j = 0; j < n; j++)
       nmod_mat_entry(B,2*nbr_man->k+i,j) =
	 nmod_mat_entry(nbr_man->U,i,j) = nmod_poly_get_coeff_ui(fq_nmod_mat_entry(u,i,j),0) % p;

#ifdef DEBUG_LEVEL_FULL
  printf("X = ");
  nmod_mat_print(nbr_man->X);
  printf("Z = ");
  nmod_mat_print(nbr_man->Z);
  printf("U = ");
  nmod_mat_print(nbr_man->U);
#endif // DEBUG_LEVEL_FULL

  // Compute the Gram matrix of the subspace with respect to the spaces
  //  we will perform the following computations upon.

  nmod_mat_init(gram, n, n, p*p);
  nmod_mat_gram(gram, B, nbr_man->quot_gram, p*p);

  // Lift Z so that it is in a hyperbolic pair with X modulo p^2.
  nmod_mat_init(Z_new, nbr_man->k, n, p*p);
  for (i = 0; i < nbr_man->k; i++) {
    for (j = 0; j < n; j++)
      nmod_mat_entry(Z_new,i,j) = nmod_mat_entry(nbr_man->Z,i,j);
    for (j = 0; j < nbr_man->k; j++) {
      delta = (i == j) ? 1 : 0;
      // we split the operation due to signed type issues
      a = nmod_mat_entry(gram, nbr_man->k-j-1, i + nbr_man->k);
      // a nonnegative value with the same residue mod p*p
      // !!! TODO - we might be able to get rid of that
      a = (a / (p*p) + 1)*p*p-a+delta;
      if (a >= (p*p))
	a -= p*p;
      for (l = 0; l < n; l++)
        nmod_mat_entry(Z_new,i,l) += a * nmod_mat_entry(nbr_man->Z,j,l);
    }
  }

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(nbr_man->Z,i,j) = nmod_mat_entry(Z_new,i,j) % (p*p);

#ifdef DEBUG_LEVEL_FULL
  printf("after setting <X,Z> = 1\n");
  printf("Z = ");
  nmod_mat_print(nbr_man->Z);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  // Verify that X and Z form a hyperbolic pair.
  // Compute the Gram matrix thusfar.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(B,nbr_man->k+i,j) = nmod_mat_entry(nbr_man->Z,i,j);

  nmod_mat_init(temp, n, n, p*p);
  nmod_mat_gram(temp, B, nbr_man->quot_gram, p*p);
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < nbr_man->k; j++)
      // This is beacuse negative % is negative
      assert((nmod_mat_entry(temp, i, nbr_man->k+j) - ((i+j == nbr_man->k-1) ? 1 : 0)) % (p*p) == 0);	
#endif // DEBUG

  if (p == 2) {
    for (i = 0; i < nbr_man->k; i++)
      for (j = 0; j < n; j++)
	nmod_mat_entry(B, nbr_man->k+i,j) = nmod_mat_entry(nbr_man->Z,i,j);
    nmod_mat_clear(gram);
    nmod_mat_init(gram, n, n, p*p*p);
    nmod_mat_gram_fmpz_mat(gram, B, nbr_man->q);
  }
  // Lift X so that it is isotropic modulo p^2.
  nmod_mat_init(X_new, nbr_man->k, n, p*p);
  half = (p*p + 1)/2;
  for (i = 0; i < nbr_man->k; i++) {
    for (l = 0; l < n; l++)
      nmod_mat_entry(X_new,i,l) = nmod_mat_entry(nbr_man->X,i,l);
    gram2 = nmod_mat_entry(gram,i,i);
    gram2 = gram2/2 + (((nmod_mat_entry(gram,i,i) % 2) == 0) ? 0 : half);
    for (j = nbr_man->k-1-i; j < nbr_man->k; j++) {
      scalar = (i+j == nbr_man->k-1) ? gram2 : nmod_mat_entry(gram,i, nbr_man->k-1-j);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      for (l = 0; l < n; l++)
	nmod_mat_entry(X_new,i,l) += scalar * nmod_mat_entry(nbr_man->Z, j, l);
    }
  }

  // we are reducing to keep the sizes small
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(nbr_man->X,i,j) = nmod_mat_entry(X_new,i,j) % (p*p);
  
#ifdef DEBUG_LEVEL_FULL
  printf("after setting <X,X> = 0\n");
  printf("X = ");
  nmod_mat_print(nbr_man->X);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  // Verify that X is isotropic modulo p^2.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(B,i,j) = nmod_mat_entry(nbr_man->X,i,j);

  // The Gram matrix on this basis.
  nmod_mat_gram(temp,B,nbr_man->quot_gram,p*p);

  // Verify all is well.
  for ( i = 0; i < nbr_man->k; i++)
    for ( j = 0; j < nbr_man->k; j++)
      assert(nmod_mat_entry(temp,i,j) % (p*p) == 0);

#endif // DEBUG

  // Lift Z so that it is isotropic modulo p^2.
  for (i = 0; i < nbr_man->k; i++) {
    for (l = 0; l < n; l++)
      nmod_mat_entry(Z_new, i, l) = nmod_mat_entry(nbr_man->Z, i, l);
    for (j = nbr_man->k-1-i; j < nbr_man->k; j++) {
      scalar = nmod_mat_entry(gram, nbr_man->k+i, 2*nbr_man->k-1-j);
      if (i + j == nbr_man->k-1)
	scalar = (scalar/2) + (((scalar % 2) == 0) ? 0 : half);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      for (l = 0; l < n; l++)
	nmod_mat_entry(Z_new,i,l) += scalar * nmod_mat_entry(nbr_man->X,j,l);
    }
  }

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(nbr_man->Z,i,j) = nmod_mat_entry(Z_new,i,j) % (p*p);

#ifdef DEBUG_LEVEL_FULL
  printf("after setting <Z,Z> = 0\n");
  printf("Z = ");
  nmod_mat_print(nbr_man->Z);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  // Verify that Z is isotropic modulo p^2.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(B,nbr_man->k+i,j) = nmod_mat_entry(nbr_man->Z,i,j);

  // The Gram matrix on this basis.
  nmod_mat_gram(temp,B, nbr_man->quot_gram,p*p);

  // Verify all is well.
  for ( i = 0; i < nbr_man->k; i++)
    for ( j = 0; j < nbr_man->k; j++)
      assert(nmod_mat_entry(temp,nbr_man->k+i,nbr_man->k+j) % (p*p) == 0);

#endif // DEBUG

  // The Gram matrix thusfar.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(B, i, j) = nmod_mat_entry(nbr_man->X,i,j);
  
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(B, nbr_man->k + i, j) = nmod_mat_entry(nbr_man->Z,i,j);

  nmod_mat_gram(gram,B, nbr_man->quot_gram,p*p);

  // Make U orthogonal to X+Z.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n - 2*nbr_man->k; j++) {
      // Clear components corresponding to X.
      scalar = nmod_mat_entry(gram, 2*nbr_man->k-1-i, 2*nbr_man->k+j);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      for (l = 0; l < n; l++)
	nmod_mat_entry(nbr_man->U,j,l) += scalar * nmod_mat_entry(nbr_man->X,i,l);

      // Clear components corresponding to Z.
      scalar = nmod_mat_entry(gram, nbr_man->k-1-i, 2*nbr_man->k+j);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      for (l = 0; l < n; l++)
	nmod_mat_entry(nbr_man->U,j,l) += scalar * nmod_mat_entry(nbr_man->Z,i,l);
    }

  // make sure that all the entries of U are between 0 and p^2
  for (i = 0; i < n - 2*nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(nbr_man->U,i,j) = nmod_mat_entry(nbr_man->U,i,j) % (p*p);

#ifdef DEBUG_LEVEL_FULL
  printf("after setting <U,X+Z> = 0\n");
  printf("U = ");
  nmod_mat_print(nbr_man->U);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  // Verify that U is now orthogonal to X+Z.
  for (i = 0; i < n-2*nbr_man->k; i++)
    for (j = 0; j < n; j++)
      nmod_mat_entry(B,2*nbr_man->k+i,j) = nmod_mat_entry(nbr_man->U,i,j);

  // The Gram matrix on this basis.
  nmod_mat_gram(temp,B, nbr_man->quot_gram,p*p);

  // Verify all is well.
  for ( i = 0; i < 2*nbr_man->k; i++)
    for ( j = 2*nbr_man->k; j < n; j++)
      assert(nmod_mat_entry(temp,i,j) % (p*p) == 0);

  nmod_mat_clear(temp);
#endif // DEBUG

  nmod_mat_clear(gram);
  nmod_mat_clear(X_new);
  nmod_mat_clear(Z_new);  
  nmod_mat_clear(B);
#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_clear(u, nbr_man->GF);
  free(excluded);
#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_clear(z, nbr_man->GF);
  free(paired);
#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_clear(x, nbr_man->GF);
#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_clear(basis_t, nbr_man->GF);
#ifdef __linux__
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif // __linux__
  fq_nmod_mat_clear(basis, nbr_man->GF);

  return;
}

void nbr_data_update_skew_matrix(nbr_data_t nbr_man, slong* row, slong* col)
{
  bool done;
  fq_nmod_t one;

  fq_nmod_init(one, nbr_man->GF);
  fq_nmod_one(one, nbr_man->GF);

#ifdef DEBUG_LEVEL_FULL
  printf("skew_matrix is: \n");
  fq_nmod_mat_print_pretty(nbr_man->p_skew, nbr_man->GF);
  printf("\n");
  printf("Updating row = %ld, col = %ld..\n", *row, *col);
#endif // DEBUG_LEVEL_FULL
  
  do {
    // Flag for determining whether we are done updating
    //  the skew matrix.
    done = true;
     
    // Increment value of the (row,col) position.
    fq_nmod_add(fq_nmod_mat_entry(nbr_man->p_skew,*row,*col), fq_nmod_mat_entry(nbr_man->p_skew,*row,*col), one, nbr_man->GF);
      
    // Update the coefficient of the skew matrix reflected
    //  across the anti-diagonal.
    fq_nmod_set(fq_nmod_mat_entry(nbr_man->p_skew, nbr_man->k-1-(*col), nbr_man->k-1-(*row)),
		fq_nmod_mat_entry(nbr_man->p_skew, *row, *col), nbr_man->GF);
    fq_nmod_neg(fq_nmod_mat_entry(nbr_man->p_skew, nbr_man->k-1-(*col), nbr_man->k-1-(*row)),
		fq_nmod_mat_entry(nbr_man->p_skew, nbr_man->k-1-(*col), nbr_man->k-1-(*row)), nbr_man->GF);
      
    // If we've rolled over, move on to the next position.
    if (fq_nmod_is_zero(fq_nmod_mat_entry(nbr_man->p_skew,*row,*col),nbr_man->GF)) {
      // The next column of our skew matrix.
      (*col)++;
      // Are we at the end of the column?
      if ((*row)+(*col) == nbr_man->k-1) {
	// Yes. Move to the next row.
	(*row)++;
	// And reset the column.
	*col = 0;
      }
      // Indicate we should repeat another iteration.
      done = false;
    }
  } while ((!done) && ((*row)+(*col) != nbr_man->k-1));

#ifdef DEBUG_LEVEL_FULL
  printf("skew_matrix is now: \n");
  fq_nmod_mat_print_pretty(nbr_man->p_skew, nbr_man->GF);
  printf("\n");
  printf("Updated row = %ld, col = %ld..\n", *row, *col);
#endif // DEBUG_LEVEL_FULL
  
  fq_nmod_clear(one, nbr_man->GF);
  return;
}

void nbr_data_update_skew_space(nbr_data_t nbr_man)
{
  slong i,j,l,n;
  slong p, val;
#ifdef DEBUG
  fmpz_mat_t B, temp, B_t, Bq;
#endif // DEBUG
  
  assert(nmod_mat_nrows(nbr_man->X) == nbr_man->k);
  assert(nmod_mat_nrows(nbr_man->Z) == nbr_man->k);
  assert(nbr_man->is_skew_init);

#ifdef DEBUG_LEVEL_FULL
  printf("X = ");
  nmod_mat_print(nbr_man->X);
  printf("Z = ");
  nmod_mat_print(nbr_man->Z);
  printf("skew = ");
  fq_nmod_mat_print(nbr_man->p_skew, nbr_man->GF);
#endif // DEBUG_LEVEL_FULL

  n = fmpz_mat_nrows(nbr_man->q);
  assert(n == fmpz_mat_ncols(nbr_man->q));
  
  p = fmpz_get_si(fq_nmod_ctx_prime(nbr_man->GF));
  // Update the skew space.
  for (i = 0; i < nbr_man->k; i++) {
    for (l = 0; l < n; l++)
      nmod_mat_entry(nbr_man->X_skew, i, l) = nmod_mat_entry(nbr_man->X, i, l);
    for (j = 0; j < nbr_man->k; j++) {
      val = nmod_poly_get_coeff_ui(fq_nmod_mat_entry(nbr_man->p_skew,i,j),0);
      for (l = 0; l < n; l++)
	nmod_mat_entry(nbr_man->X_skew,i,l) += p * val * nmod_mat_entry(nbr_man->Z,j,l);
    }
    for (l = 0; l < n; l++)
      nmod_mat_entry(nbr_man->X_skew, i, j) %= (p*p);
  }

#ifdef DEBUG_LEVEL_FULL
  printf("X_skew = ");
  nmod_mat_print(nbr_man->X_skew);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  fmpz_mat_init(B, n, n);
  fmpz_mat_init(B_t, n, n);
  fmpz_mat_init(Bq, n, n);
  fmpz_mat_init(temp, n, n);

  // Verify that X_skew is isotropic modulo p^2.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      fmpz_set_si(fmpz_mat_entry(B,i,j),nmod_mat_entry(nbr_man->X_skew, i, j));

  fmpz_mat_transpose(B_t, B);
  fmpz_mat_mul(Bq, B, nbr_man->q);
  // The Gram matrix on this basis.
  fmpz_mat_mul(temp, Bq, B_t);
  
  // Verify all is well.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < nbr_man->k; j++)
      assert(fmpz_get_si(fmpz_mat_entry(temp,i,j)) % (p*p) == 0);

  fmpz_mat_clear(B);
  fmpz_mat_clear(B_t);
  fmpz_mat_clear(Bq);
  fmpz_mat_clear(temp);
#endif // DEBUG
  
  return;
}

void nbr_data_get_next_neighbor(nbr_data_t nbr_man)
{
  slong row, col;

  // The starting position of the skew vector to update.
  row = 0;
  col = 0;
  // Update the skew matrix (only for k >= 2).
  if (nbr_man->skew_dim != 0) {
    nbr_data_update_skew_matrix(nbr_man, &row, &col);
  }

  // If we haven't rolled over, update the skew space and return...
  if (row + col < nbr_man->k-1) {
    nbr_data_update_skew_space(nbr_man);
    return;
  }

  // ...otherwise, get the next isotropic subspace modulo p.
  nbr_data_next_isotropic_subspace(nbr_man);

  // Lift the subspace if we haven't reached the end of the list.
  if (!(nbr_man->is_done)) {
    nbr_data_lift_subspace(nbr_man);
    nmod_mat_set(nbr_man->X_skew, nbr_man->X);
  }
  
  return;
}

void fmpz_quad_transform(fmpz_mat_t new_gram, const fmpz_mat_t gram, const fmpz_mat_t isom, slong scale)
{
  fmpz_mat_t isom_gram, isom_t;
  slong n;

  n = fmpz_mat_nrows(gram);
  assert(n == fmpz_mat_ncols(gram));
  
  fmpz_mat_init(isom_gram, n, n);
  fmpz_mat_init(isom_t, n, n);
  
  fmpz_mat_transpose(isom_t, isom);
  //  fmpz_mat_mul(isom_gram, isom, gram);
  // fmpz_mat_mul(new_gram, isom_gram, isom_t);
  fmpz_mat_mul(isom_gram, isom_t, gram);
  fmpz_mat_mul(new_gram, isom_gram, isom);

  fmpz_mat_scalar_divexact_si(new_gram, new_gram, scale*scale);
  
  fmpz_mat_clear(isom_gram);
  fmpz_mat_clear(isom_t);
  
  return;
}
// !! TODO - replace by the already implemented hnf in flint
void fmpz_mat_hermite_form(fmpz_mat_t H, const fmpz_mat_t s, slong d)
{
  slong row, pivot, col, j, n;
  fmpz_t a,b,x,y,g,q,q_a,q_b, prod1, prod2, scalar;
  fmpz_mat_t g_h_prime, b_primes;

  fmpz_init(a);
  fmpz_init(b);
  fmpz_init(x);
  fmpz_init(y);
  fmpz_init(g);
  fmpz_init(q);
  fmpz_init(q_a);
  fmpz_init(q_b);
  fmpz_init(prod1);
  fmpz_init(prod2);
  fmpz_init(scalar);

  n = fmpz_mat_nrows(s);
  assert(n == fmpz_mat_ncols(s));
  
  fmpz_mat_init(g_h_prime, 1, n);
  fmpz_mat_init(b_primes, 1, n);
  
  fmpz_mat_one(H);
  fmpz_mat_scalar_mul_si(H, H, d);
  for (row = 0; row < n; row++) {
    for (col = 0; col < n; col++) {
      fmpz_set(fmpz_mat_entry(b_primes,0,col), fmpz_mat_entry(s, row, col));
    }
    // Here we compute the HNF of H and this row and store it in H
    for (pivot = 0; pivot < n; pivot++) {
      fmpz_set(a,fmpz_mat_entry(H,pivot,pivot));
      fmpz_set(b,fmpz_mat_entry(b_primes,0,pivot));
      fmpz_xgcd(g, x, y, a, b);
      fmpz_divexact(q_a, a, g);
      fmpz_divexact(q_b, b, g);
      for (col = 0; col < n; col++) {
	fmpz_mul(prod1, x, fmpz_mat_entry(H, pivot, col));
	fmpz_mul(prod2, y, fmpz_mat_entry(b_primes, 0, col));
	fmpz_add(fmpz_mat_entry(g_h_prime,0,col),prod1,prod2);
      }
      for (col = 0; col < n; col++) {
	fmpz_mul(prod1, q_a, fmpz_mat_entry(b_primes, 0, col));
	fmpz_mul(prod2, q_b, fmpz_mat_entry(H, pivot, col));
	fmpz_sub(fmpz_mat_entry(b_primes,0,col),prod1,prod2);
      }
      for (j = pivot; j < n; j++) {
	// !! TODO - check if this should be divexact
	fmpz_tdiv_q(scalar, fmpz_mat_entry(b_primes, 0, j), fmpz_mat_entry(H, j, j));
	for (col = 0; col < n; col++) {
	  fmpz_mul(prod1, scalar, fmpz_mat_entry(H,j,col));
	  fmpz_sub(fmpz_mat_entry(b_primes,0,col), fmpz_mat_entry(b_primes,0,col), prod1);
	}
      }
      for (col = 0; col < n; col++)
	fmpz_set(fmpz_mat_entry(H, pivot, col), fmpz_mat_entry(g_h_prime,0,col));
    }
    for (pivot = n-1; pivot > 0; pivot--) {
      for (col = pivot; col < n; col++) {
	fmpz_tdiv_q(q, fmpz_mat_entry(H, pivot-1, col), fmpz_mat_entry(H, col, col));
	for (j = col; j < n; j++) {
	  fmpz_mul(prod1, q, fmpz_mat_entry(H, col, j));
	  fmpz_sub(fmpz_mat_entry(H, pivot-1, j), fmpz_mat_entry(H, pivot-1, j), prod1);
	}
      }
    }
  }

  fmpz_mat_clear(b_primes);
  fmpz_mat_clear(g_h_prime);

  fmpz_clear(scalar);
  fmpz_clear(prod2);
  fmpz_clear(prod1);
  fmpz_clear(a);
  fmpz_clear(b);
  fmpz_clear(x);
  fmpz_clear(y);
  fmpz_clear(g);
  fmpz_clear(q);
  fmpz_clear(q_a);
  fmpz_clear(q_b);
  fmpz_clear(prod1);
  fmpz_clear(prod2);
  return;
}

void nbr_data_build_neighbor(fmpz_mat_t nbr, fmpz_mat_t s, const nbr_data_t nbr_man)
{
  slong p, p2, p3;
  slong i, j, n;
  fmpz_mat_t hermite;

  p = fmpz_get_si(fq_nmod_ctx_prime(nbr_man->GF));
  p2 = p*p;
  p3 = p2*p;

  n = fmpz_mat_nrows(nbr_man->q);
  assert(n == fmpz_mat_ncols(nbr_man->q));
  
  fmpz_mat_init(hermite, n, n);
  // fill the isomtery by the X,Z,U
  // if lift_subspace was successful,
  // <X,X>, <X,U>,<Z,Z>,<Z,U> in p^2 and <X,Z> = 1 mod p^2

  // For now, we follow the magma implementation
  // start by scaling the basis

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      fmpz_set_si(fmpz_mat_entry(s,i,j),nmod_mat_entry(nbr_man->X_skew,i,j));

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < n; j++)
      fmpz_set_si(fmpz_mat_entry(s,nbr_man->k+i,j), p2 * nmod_mat_entry(nbr_man->Z,i,j));

  for (i = 0; i < n - 2  * nbr_man->k; i++)
    for (j = 0; j < n; j++)
      fmpz_set_si(fmpz_mat_entry(s,2*nbr_man->k+i,j), p * nmod_mat_entry(nbr_man->U,i,j));

#ifdef DEBUG_LEVEL_FULL
  printf("Before hnf, s = \n");
  fmpz_mat_print_pretty(s);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  // !! TODO - replace by the already implemented hnf in flint
  fmpz_mat_hermite_form(hermite, s, p3);

#ifdef DEBUG_LEVEL_FULL
  printf("After hnf, hermite = \n");
  fmpz_mat_print_pretty(hermite);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  fmpz_mat_transpose(s, hermite);

  // need to adjust determinant for s to be in SO
  // This transforms using the isometry (and rescales already)
  fmpz_quad_transform(nbr, nbr_man->q, s, p);

  fmpz_mat_clear(hermite);
  
  return;
}

bool nbr_data_has_ended(const nbr_data_t nbr_man)
{
  return nbr_man->is_done;
}


slong number_of_isotropic_subspaces(slong q, slong r, slong a, slong f, slong k, bool rad_cnt)
{
  fmpz_t q_z, prod1, prod;
  slong i, num;

  fmpz_init_set_si(q_z,q);
  fmpz_init(prod1);
  fmpz_init(prod);

  fmpz_pow_ui(prod, q_z, k*f);
  
  for (i = 1; i <= k; i++) {
    fmpz_pow_ui(prod1, q_z, r-i+1);
    fmpz_sub_si(prod1, prod1, 1);
    fmpz_mul(prod, prod, prod1);
    fmpz_pow_ui(prod1, q_z, r+a-i);
    fmpz_add_si(prod1, prod1, 1);
    fmpz_mul(prod, prod, prod1);
  }

  for (i = 1; i <= k; i++) {
    fmpz_pow_ui(prod1, q_z, i);
    fmpz_sub_si(prod1, prod1, 1);
    fmpz_divexact(prod, prod, prod1);
  }

  num = fmpz_get_si(prod);

  if (rad_cnt) {
    // in this case we also consider the isotropic subspaces in the radical
    // these are all the k-dimensionl spaces inside the radical
    fmpz_one(prod);
    for (i = 1; i <= k; i++) {
      fmpz_pow_ui(prod1, q_z, f-i+1);
      fmpz_sub_si(prod1, prod1, 1);
      fmpz_mul(prod, prod, prod1);
    }
    for (i = 1; i <= k; i++) {
      fmpz_pow_ui(prod1, q_z, i);
      fmpz_sub_si(prod1, prod1, 1);
      fmpz_divexact(prod, prod, prod1);
    }
    num += fmpz_get_si(prod);
  }
  
  fmpz_clear(q_z);
  fmpz_clear(prod1);
  fmpz_clear(prod);
  
  return num;
}

slong number_of_neighbors(const nbr_data_t nbr_man, bool rad_cnt)
{
  slong p, num;

  p = fmpz_get_si(fq_nmod_ctx_prime(nbr_man->GF));
  num = number_of_isotropic_subspaces(p, nbr_man->witt_index, nbr_man->aniso_dim,
				      nbr_man->rad_dim, nbr_man->k, rad_cnt);
  if (nbr_man->k == 2)
    num *= p;

  return num;
}

