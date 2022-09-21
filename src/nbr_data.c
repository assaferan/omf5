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

void nbr_data_init(nbr_data_t nbr_man, matrix_TYP* q, slong p_int, slong k)
{
  slong idx;
  fmpz_t p;
#ifdef DEBUG
  fq_nmod_t value;
#endif // DEBUG

  fmpz_init_set_si(p, p_int);
  fmpz_mat_init_set_matrix_TYP(nbr_man->q, q);
  fmpz_init(nbr_man->disc);
  fmpz_mat_det(nbr_man->disc, nbr_man->q);
  // (half) discriminant = det/2 when N = 5 
  fmpz_divexact_si(nbr_man->disc, nbr_man->disc, 2);

  fq_nmod_ctx_init(nbr_man->GF, p, 1, "1");
  fq_nmod_mat_init_set_fmpz_mat(nbr_man->b, nbr_man->q, nbr_man->GF);
  nmod_mat_init_set_fmpz_mat(nbr_man->quot_gram, nbr_man->q, p_int*p_int);

  fq_nmod_mat_init(nbr_man->vec, 1, N, nbr_man->GF);
#ifdef DEBUG
  fq_nmod_quad_isotropic_vector(nbr_man->vec, nbr_man->b, nbr_man->GF, 0, true);
#else
  fq_nmod_quad_isotropic_vector(nbr_man->vec, nbr_man->b, nbr_man->GF, 0, false);
#endif // DEBUG

#ifdef DEBUG
  fq_nmod_init(value, nbr_man->GF);
  fq_nmod_quad_evaluate(value, nbr_man->b, nbr_man->vec, nbr_man->GF);
  assert(fq_nmod_is_zero(value, nbr_man->GF));
#endif // DEBUG  

  fq_nmod_mat_init(nbr_man->p_std_gram, N, N, nbr_man->GF);
  fq_nmod_mat_init(nbr_man->p_basis, N, N, nbr_man->GF);
  fq_nmod_mat_init(nbr_man->p_skew, k, k, nbr_man->GF);

  nbr_man->is_skew_init = true;

#ifdef DEBUG
  fq_nmod_quad_decompose(nbr_man->p_std_gram, nbr_man->p_basis, nbr_man->b, nbr_man->GF, true);
#else
  fq_nmod_quad_decompose(nbr_man->p_std_gram, nbr_man->p_basis, nbr_man->b, nbr_man->GF, false);  
#endif // DEBUG

  fq_nmod_mpoly_ctx_init(nbr_man->p_q_std_ctx, N, ORD_DEGREVLEX, nbr_man->GF);
  fq_nmod_mpoly_init(nbr_man->p_q_std,nbr_man->p_q_std_ctx);
  fq_nmod_poly_set_fq_nmod_quad(nbr_man->p_q_std, nbr_man->p_std_gram, nbr_man->GF, nbr_man->p_q_std_ctx);

#ifdef DEBUG_LEVEL_FULL
  printf("Performed Witt Decomposition on\n");
  fq_nmod_mat_print_pretty(nbr_man->b,nbr_man->GF);
  printf("Resulting gram matrix is \n");
  fq_nmod_mat_print_pretty(nbr_man->p_std_gram,nbr_man->GF);
  printf("Resulting basis is \n");
  fq_nmod_mat_print_pretty(nbr_man->p_basis,nbr_man->GF);
#endif // DEBUG_LEVEL_FULL

  // Count the rows at the end of the matrix which are exactly zero.
  idx = N;
  while ((idx >= 1) && fq_nmod_mat_is_zero_row(nbr_man->p_std_gram,idx-1,nbr_man->GF)) idx--;

  // The dimension of the radical.
  nbr_man->rad_dim = N - idx;

  // Determine the dimension of the totally hyperbolic subspace.
  idx = 1;
  while ((idx <= N - nbr_man->rad_dim) && fq_nmod_is_zero(fq_nmod_mat_entry(nbr_man->p_std_gram,idx-1,idx-1),nbr_man->GF)) idx++;

  // Dimension of the anistotropic subspace.
  nbr_man->aniso_dim = N - nbr_man->rad_dim - idx + 1;

  // The number of hyperbolic planes in the Witt decomposition.
  nbr_man->witt_index = (idx - 1) / 2;

  pivot_data_init(nbr_man->pivots, N - nbr_man->rad_dim, nbr_man->aniso_dim, k);
  nbr_man->k = k;
  nbr_man->skew_dim = k*(k-1)/2;
  nbr_man->is_done = false;

  nbr_data_next_isotropic_subspace(nbr_man);
  nbr_data_lift_subspace(nbr_man);

  nmod_mat_init_set(nbr_man->X_skew, nbr_man->X);
  
#ifdef DEBUG
  fq_nmod_clear(value, nbr_man->GF);
#endif // DEBUG
  fmpz_clear(p);
  
  return;
}

void nbr_data_clear(nbr_data_t nbr_man)
{
  nmod_mat_clear(nbr_man->X_skew);
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
  slong i, j, pos;
  bool all_zero;
  fq_nmod_t one;
  
  if (nbr_man->pivots->num_params == 0) {
    // Move to the next pivot.
    nbr_man->pivots->pivot_ptr++;

    // If we've exceeded the list of pivots, we're done.
    if (nbr_man->pivots->pivot_ptr > nbr_man->pivots->total_len) {
      // Reset the pivot pointer so that we can loop through
      //  the isotropic subspaces again if needed.
      nbr_man->pivots->pivot_ptr = 0;
      fq_nmod_mat_clear(nbr_man->iso_subspace, nbr_man->GF);
      nbr_man->is_done = true;
      return;
    }
    // Initialize the new pivot.
    pivot_data_params_clear(nbr_man->pivots);
    nbr_data_params_init(nbr_man->pivots, nbr_man);
  }

  // The list of evaluation values.
  // fq_nmod_mat_init(eval_list, 1, N * nbr_man->k, nbr_man->GF);
  // fq_nmod_mat_zero(eval_list, nbr_man->GF);
  eval_list = (fq_nmod_struct**)malloc(N*nbr_man->k*sizeof(fq_nmod_struct*));
  for (i = 0; i < N * nbr_man->k; i++) {
    fq_nmod_init(eval_list[i], nbr_man->GF);
    fq_nmod_zero(eval_list[i], nbr_man->GF);
  }

  // Produce the isotropic subspace corresponding to the current
  //  parameters.

  for (i = 0; i < nbr_man->pivots->num_params; i++)
    // fq_nmod_set(fq_nmod_mat_entry(eval_list,0,i), nbr_man->pivots->params[i], nbr_man->GF);
    fq_nmod_set(eval_list[i], nbr_man->pivots->params[i], nbr_man->GF);

  // The basis for the current isotropic subspace.
  fq_nmod_mat_init(nbr_man->iso_subspace, nbr_man->k, N, nbr_man->GF);
  for (i = 0; i < nbr_man->k; i++) {
    for (j = 0; j < N; j++) {
      fq_nmod_mpoly_evaluate_all_fq_nmod(fq_nmod_mat_entry(nbr_man->iso_subspace,i,j),
					 fq_nmod_mpoly_mat_entry(nbr_man->pivots->p_isotropic_param,i,j),
					 eval_list, nbr_man->pivots->R);
    }
  }

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
    for (i = 0; i < nbr_man->pivots->num_params; i++) {
      fq_nmod_clear(nbr_man->pivots->params[i], nbr_man->GF);
    }
    free(nbr_man->pivots->params);
    nbr_man->pivots->num_params = 0;
  }
  
  // fq_nmod_mat_clear(eval_list, nbr_man->GF);
  for (i = 0; i < N * nbr_man->k; i++)
    fq_nmod_clear(eval_list[i], nbr_man->GF);
  free(eval_list);
  
  return;
}

void nbr_data_params_init(pivot_data_t pivots, const nbr_data_t nbr_man)
{
  slong* pivot;
  slong rank, row, col, remove_idx, remove_size, rows, i, j, r, c, data_idx, vec_idx;
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
  
  pivot = pivots->pivots[pivots->pivot_ptr-1];
  rank = nbr_man->k * N;

  // Initialize matrix that will determine parameterization.
  fq_nmod_mpoly_ctx_init(pivots->R, rank, ORD_DEGREVLEX, nbr_man->GF);
  /* fq_nmod_mpoly_mat_init(data, 1, rank, pivots->R); */
  /* for (i = 0; i < rank; i++) */
  /*   fq_nmod_mpoly_gen(fq_nmod_mpoly_mat_entry(data,0,i),i,pivots->R); */

  /* fq_nmod_mpoly_mat_init_set_fq_nmod_mpoly_vec(pivots->p_isotropic_param, data, nbr_man->k, N, pivots->R); */

  fq_nmod_mpoly_mat_init(pivots->p_isotropic_param, nbr_man->k, N, pivots->R);
  i = 0;
  for (row = 0; row < nbr_man->k; row++)
    for (col = 0; col < N; col++)
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
      remove[remove_idx++] = row*N+pivot[col];
    }

  // Clear the rows prior to the pivot positions (but not the radical).
  for (row = 0; row < nbr_man->k; row++)
    for (col = 0; col < pivot[row]; col++) {
      fq_nmod_mpoly_zero(fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param, row, col), pivots->R);
      remove[remove_idx++] = row*N+col;
    }

  // Check if one or more of the anisotropic coordinates need to be zero.
  for (row = 0; row < nbr_man->k; row++) {
    if (pivot[row] >= nbr_man->witt_index) {
      for (col = 0; col < nbr_man->aniso_dim; col++) {
	fq_nmod_mpoly_zero(fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param, row, N-1-nbr_man->rad_dim-col), pivots->R);
	remove[remove_idx++] = (row+1)*N - 1 - nbr_man->rad_dim-col;
      }
    }
  }

  // Here we will save the quadratic polynomials
  fq_nmod_mpoly_mat_init(data, 1, rows, pivots->R);

  // The matrix that we're going to echelonize.
  fq_nmod_mat_init(mat, rows, rank, nbr_man->GF);

  // fq_nmod_mpoly_mat_init(vec, 1, rows*N, pivots->R);
  vec = (fq_nmod_mpoly_struct**)malloc(rows*N*sizeof(fq_nmod_mpoly_struct*));
  for (i = 0; i < rows*N; i++)
    fq_nmod_mpoly_init(vec[i], pivots->R);
  
  fq_nmod_mat_init(l, 1, rank, nbr_man->GF);

  data_idx = 0;
  vec_idx = 0;
  row = 0;
  for (i = 0; i < nbr_man->k; i++)
    for (j = i; j < nbr_man->k; j++) {
      // The appropriate vector that we want to be isotropic.
      for (r = 0; r < N; r++) {
	fq_nmod_mpoly_set(vec[vec_idx], fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,i,r), pivots->R);
	if (i != j)
	  fq_nmod_mpoly_add(vec[vec_idx],
			    vec[vec_idx], fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,j,r), pivots->R);
	vec_idx++;
      }
      fq_nmod_mpoly_compose_fq_nmod_mpoly(f, nbr_man->p_q_std, vec, nbr_man->p_q_std_ctx, pivots->R);
      // Degree 2 terms are inhomogenous.
      fq_nmod_mpoly_quadratic_part(f, f, pivots->R);
      fq_nmod_mpoly_neg(fq_nmod_mpoly_mat_entry(data,0,data_idx++), f, pivots->R);
      // The other terms are linear
      // so fill in the matrix accordingly.
      fq_nmod_mpoly_linear_part(l, f, pivots->R);
      for (r = 0; r < rank; r++)
	fq_nmod_set(fq_nmod_mat_entry(mat, row, r), fq_nmod_mat_entry(l, 0, r), nbr_man->GF);
      // Move on to the next row.
      row++;
      // Silly question. Isn't row == i?
    }

#ifdef DEBUG_LEVEL_FULL
  printf("The matrix before echelon is mat = \n");
  fq_nmod_mat_print_pretty(mat, nbr_man->GF);
  printf("The last entry is the quadratic data = \n");
  for (i = 0; i < rows; i++) {
    fq_nmod_mpoly_print_pretty(fq_nmod_mpoly_mat_entry(data,0,i), var_names, pivots->R);
    printf(" ");
  }
  printf("\n");
#endif

  // Compute the Echelon form of the matrix.
  fq_nmod_mat_init(trans, rows, rows, nbr_man->GF);
  // fq_nmod_mat_rref(mat, trans, nbr_man->GF);
  fq_nmod_mat_rref_trans(mat, trans, nbr_man->GF);

#ifdef DEBUG_LEVEL_FULL
  printf("The matrix after echelon is mat = \n");
  fq_nmod_mat_print_pretty(mat, nbr_man->GF);
#endif

  // The evaluation list for replacing variables with their dependence
  //  relations.

  // fq_nmod_mpoly_mat_init(eval_list, 1, rank, pivots->R);
  eval_list = (fq_nmod_mpoly_struct**)malloc(rank*sizeof(fq_nmod_mpoly_struct*));
  for (i = 0; i < rank; i++) {
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
#endif // DEBUG_LEVEL_FULL

  vec_idx = 0;
  for (i = 0; i < nbr_man->k; i++)
    for (j = i; j < nbr_man->k; j++) {
      for (r = 0; r < N; r++) {
	fq_nmod_mpoly_set(vec[vec_idx], fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,i,r), pivots->R);
	if (i != j)
	  fq_nmod_mpoly_add(vec[vec_idx],
			    vec[vec_idx], fq_nmod_mpoly_mat_entry(pivots->p_isotropic_param,j,r), pivots->R);
	vec_idx++;
      }
  
#ifdef DEBUG_LEVEL_FULL
      printf("Substituting vec = ");
      for (i = 0; i < rows*N; i++)
	fq_nmod_mpoly_print_pretty(vec[i], var_names, pivots->R);
      printf(" in q_std = ");
      fq_nmod_mpoly_print_pretty(nbr_man->p_q_std, var_names, pivots->R);
#endif // DEBUG_LEVEL_FULL
      
      fq_nmod_mpoly_compose_fq_nmod_mpoly(f, nbr_man->p_q_std, vec, nbr_man->p_q_std_ctx, pivots->R);
#ifdef DEBUG_LEVEL_FULL
      printf(" yields f = ");
      fq_nmod_mpoly_print_pretty(f, var_names, pivots->R);
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

  for (i = 0; i < rank; i++)
    fq_nmod_mpoly_clear(eval_list[i], pivots->R);
  
  free(eval_list);
  
  fq_nmod_mpoly_clear(g, pivots->R);

  fq_nmod_mat_clear(trans, nbr_man->GF);
  
  fq_nmod_mat_clear(l, nbr_man->GF);
  fq_nmod_mpoly_clear(f, pivots->R);

  for (i = 0; i < rows*N; i++)
    fq_nmod_mpoly_clear(vec[i], pivots->R);
  free(vec);
  
  fq_nmod_mat_clear(mat,nbr_man->GF);
    
  free(remove);

  fq_nmod_mpoly_mat_clear(data, pivots->R);
  
  return;
}

void nmod_mat_gram(nmod_mat_t gram, const nmod_mat_t B, const nmod_mat_t q, slong n)
{
  nmod_mat_t Bq, B_t;
  slong i,j;
  
  nmod_mat_init(Bq, N, N, n);
  nmod_mat_init(B_t, N, N, n);

  nmod_mat_mul(Bq, B, q);
  nmod_mat_transpose(B_t, B);
  nmod_mat_mul(gram, Bq, B_t);

  // !! TODO - this isn't necessary only for debugging versus magma
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      nmod_mat_entry(gram,i,j) %= n;
  
  nmod_mat_clear(B_t);
  nmod_mat_clear(Bq);
  return;
}

void nmod_mat_gram_fmpz_mat(nmod_mat_t gram, const nmod_mat_t B, const fmpz_mat_t q)
{
  fmpz_mat_t B_fmpz, Bq, B_t, gram_fmpz;
  slong i,j;

  fmpz_mat_init(B_fmpz, N, N);
  fmpz_mat_init(gram_fmpz, N, N);
  fmpz_mat_init(Bq, N, N);
  fmpz_mat_init(B_t, N, N);

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      fmpz_set_si(fmpz_mat_entry(B_fmpz,i,j), nmod_mat_entry(B,i,j));

  fmpz_mat_mul(Bq, B_fmpz, q);
  fmpz_mat_transpose(B_t, B_fmpz);
  fmpz_mat_mul(gram_fmpz, Bq, B_t);

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
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
  slong num_pivots;
  slong i, j, l, u_idx;
  slong h_dim, delta, a, half, gram2, scalar;
  fq_nmod_mat_t basis, basis_t, x, z, u;
  bool* excluded;
  nmod_mat_t B, X_new, Z_new;
  // fmpz_mat_t B;
  nmod_mat_t gram;

#ifdef DEBUG
  nmod_mat_t temp;
#endif // DEBUG
  
  if (nbr_man->is_done) return;

  p = fmpz_get_si(fq_nmod_ctx_prime(nbr_man->GF));

  assert(nbr_man->pivots->pivot_ptr >= 1);

  // Get the pivots for the bases of the isotropic subspaces.
  pivots = nbr_man->pivots->pivots[nbr_man->pivots->pivot_ptr-1];
  num_pivots = nbr_man->pivots->pivot_lens[nbr_man->pivots->pivot_ptr-1];

#ifdef DEBUG_LEVEL_FULL
  printf("before lifting, p_basis is \n");
  fq_nmod_mat_print(nbr_man->p_basis, nbr_man->GF);
  printf("iso_subspace is \n");
  fq_nmod_mat_print(nbr_man->iso_subspace, nbr_man->GF);
  printf("pivots = ");
  for (i = 0; i < num_pivots; i++)
    printf("%ld", pivots[i]);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  fq_nmod_mat_init_set(basis, nbr_man->p_basis, nbr_man->GF);
  // Set up the correct basis vectors.
  for (i = 0; i < nbr_man->k; i++)
    for (j = pivots[i] + 1; j < N; j++)
      fq_nmod_mat_add_col(basis, pivots[i], j, fq_nmod_mat_entry(nbr_man->iso_subspace,i,j), nbr_man->GF);
  
#ifdef DEBUG_LEVEL_FULL
  printf("the correct basis vectors are \n");
  fq_nmod_mat_print(basis, nbr_man->GF);
#endif // DEBUG_LEVEL_FULL

  fq_nmod_mat_init(basis_t, fq_nmod_mat_ncols(basis, nbr_man->GF), fq_nmod_mat_nrows(basis, nbr_man->GF), nbr_man->GF);
  fq_nmod_mat_transpose(basis_t, basis, nbr_man->GF);
  
  // Extract our target isotropic subspace modulo p
  fq_nmod_mat_init(x, nbr_man->k, N, nbr_man->GF);
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
      fq_nmod_set(fq_nmod_mat_entry(x,i,j), fq_nmod_mat_entry(basis_t, pivots[i],j), nbr_man->GF);

#ifdef DEBUG_LEVEL_FULL
  printf("x = ");
  fq_nmod_mat_print(x, nbr_man->GF);
#endif // DEBUG_LEVEL_FULL

  // Extract the hyperbolic complement modulo p.
  fq_nmod_mat_init(z, nbr_man->k, N, nbr_man->GF);
  paired = (slong*)malloc((nbr_man->k)*sizeof(slong));
  h_dim = 2 * nbr_man->witt_index;
  for (i = 0 ; i < nbr_man->k; i++)
    paired[i] = h_dim - 1 - pivots[nbr_man->k-1-i];
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
      fq_nmod_set(fq_nmod_mat_entry(z,i,j), fq_nmod_mat_entry(basis_t, paired[i],j), nbr_man->GF);

#ifdef DEBUG_LEVEL_FULL
  printf("z = ");
  fq_nmod_mat_print(z, nbr_man->GF);
#endif // DEBUG_LEVEL_FULL

  // Extract the remaining basis vectors.
  excluded = (bool*)malloc(N * sizeof(bool));
  for (i = 0; i < N; i++)
    excluded[i] = true;
  for (i = 0; i < nbr_man->k; i++) {
    excluded[pivots[i]] = false;
    excluded[paired[i]] = false;
  }
  fq_nmod_mat_init(u, N-2*nbr_man->k, N, nbr_man->GF);
  u_idx = 0;
  for (i = 0; i < N; i++)
    if (excluded[i])
       for (j = 0; j < N; j++)
	 fq_nmod_set(fq_nmod_mat_entry(u,u_idx++,j), fq_nmod_mat_entry(basis_t, i,j), nbr_man->GF);

#ifdef DEBUG_LEVEL_FULL
  printf("u = ");
  fq_nmod_mat_print(u, nbr_man->GF);
#endif // DEBUG_LEVEL_FULL

  // Convert to coordinates modulo p^2.
  nmod_mat_init(nbr_man->X, nbr_man->k, N, p*p);
  nmod_mat_init(nbr_man->Z, nbr_man->k, N, p*p);
  nmod_mat_init(nbr_man->U, N-2*nbr_man->k, N, p*p);
  nmod_mat_init(B, N, N, p*p);
  // fmpz_mat_init(nbr_man->X, nbr_man->k, N);
  // fmpz_mat_init(nbr_man->Z, nbr_man->k, N);
  // fmpz_mat_init(nbr_man->U, N-2*nbr_man->k, N);
  // fmpz_mat_init(B, N, N);
  
  // Build the coordinate matrix.
  // !! TODO - the mod p is not necessary, good for debugging
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++) {
      nmod_mat_entry(B,i,j) = nmod_mat_entry(nbr_man->X,i,j) = nmod_poly_get_coeff_ui(fq_nmod_mat_entry(x,i,j),0) % p;
    }

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
       nmod_mat_entry(B,nbr_man->k+i,j) =
	 nmod_mat_entry(nbr_man->Z,i,j) = nmod_poly_get_coeff_ui(fq_nmod_mat_entry(z,i,j),0) % p;

  for (i = 0; i < N-2*nbr_man->k; i++)
    for (j = 0; j < N; j++)
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

  nmod_mat_init(gram, N, N, p*p);
  nmod_mat_gram(gram, B, nbr_man->quot_gram, p*p);

  // Lift Z so that it is in a hyperbolic pair with X modulo p^2.
  nmod_mat_init(Z_new, nbr_man->k, N, p*p);
  for (i = 0; i < nbr_man->k; i++) {
    for (j = 0; j < N; j++)
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
      for (l = 0; l < N; l++)
        nmod_mat_entry(Z_new,i,l) += a * nmod_mat_entry(nbr_man->Z,j,l);
    }
  }

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
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
    for (j = 0; j < N; j++)
      nmod_mat_entry(B,nbr_man->k+i,j) = nmod_mat_entry(nbr_man->Z,i,j);

  nmod_mat_init(temp, nbr_man->k, nbr_man->k, p*p);
  nmod_mat_gram(temp, B, nbr_man->quot_gram, p*p);
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < nbr_man->k; j++)
      // This is beacuse negative % is negative
      assert((nmod_mat_entry(temp, i, nbr_man->k+j) - ((i+j == nbr_man->k-1) ? 1 : 0)) % (p*p) == 0);	
#endif // DEBUG

  if (p == 2) {
    for (i = 0; i < nbr_man->k; i++)
      for (j = 0; j < N; j++)
	nmod_mat_entry(B, nbr_man->k+i,j) = nmod_mat_entry(nbr_man->Z,i,j);
    nmod_mat_clear(gram);
    nmod_mat_init(gram, N, N, p*p*p);
    nmod_mat_gram_fmpz_mat(gram, B, nbr_man->q);
  }
  // Lift X so that it is isotropic modulo p^2.
  nmod_mat_init(X_new, nbr_man->k, N, p*p);
  half = (p*p + 1)/2;
  for (i = 0; i < nbr_man->k; i++) {
    for (l = 0; l < N; l++)
      nmod_mat_entry(X_new,i,l) = nmod_mat_entry(nbr_man->X,i,l);
    gram2 = nmod_mat_entry(gram,i,i);
    gram2 = gram2/2 + (((nmod_mat_entry(gram,i,i) % 2) == 0) ? 0 : half);
    for (j = nbr_man->k-1-i; j < nbr_man->k; j++) {
      scalar = (i+j == nbr_man->k-1) ? gram2 : nmod_mat_entry(gram,i, nbr_man->k-1-j);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      for (l = 0; l < N; l++)
	nmod_mat_entry(X_new,i,l) += scalar * nmod_mat_entry(nbr_man->Z, j, l);
    }
  }

  // we are reducing to keep the sizes small
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
      nmod_mat_entry(nbr_man->X,i,j) = nmod_mat_entry(X_new,i,j) % (p*p);
  
#ifdef DEBUG_LEVEL_FULL
  printf("after setting <X,X> = 0\n");
  printf("X = ");
  nmod_mat_print(nbr_man->X);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  // Verify that X is isotropic modulo p^2.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
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
    for (l = 0; l < N; l++)
      nmod_mat_entry(Z_new, i, l) = nmod_mat_entry(nbr_man->Z, i, l);
    for (j = nbr_man->k-1-i; j < nbr_man->k; j++) {
      scalar = nmod_mat_entry(gram, nbr_man->k+i, 2*nbr_man->k-1-j);
      if (i + j == nbr_man->k-1)
	scalar = (scalar/2) + (((scalar % 2) == 0) ? 0 : half);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      for (l = 0; l < N; l++)
	nmod_mat_entry(Z_new,i,l) += scalar * nmod_mat_entry(nbr_man->X,j,l);
    }
  }

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
      nmod_mat_entry(nbr_man->Z,i,j) = nmod_mat_entry(Z_new,i,j) % (p*p);

#ifdef DEBUG_LEVEL_FULL
  printf("after setting <Z,Z> = 0\n");
  printf("Z = ");
  nmod_mat_print(nbr_man->Z);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  // Verify that Z is isotropic modulo p^2.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
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
    for (j = 0; j < N; j++)
      nmod_mat_entry(B, i, j) = nmod_mat_entry(nbr_man->X,i,j);
  
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
      nmod_mat_entry(B, nbr_man->k + i, j) = nmod_mat_entry(nbr_man->Z,i,j);

  nmod_mat_gram(gram,B, nbr_man->quot_gram,p*p);

  // Make U orthogonal to X+Z.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N - 2*nbr_man->k; j++) {
      // Clear components corresponding to X.
      scalar = nmod_mat_entry(gram, 2*nbr_man->k-1-i, 2*nbr_man->k+j);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      for (l = 0; l < N; l++)
	nmod_mat_entry(nbr_man->U,j,l) += scalar * nmod_mat_entry(nbr_man->X,i,l);

      // Clear components corresponding to Z.
      scalar = nmod_mat_entry(gram, nbr_man->k-1-i, 2*nbr_man->k+j);
      scalar = (scalar / (p*p) + 1)*p*p-scalar;
      if (scalar >= p*p)
	scalar -= p*p;
      for (l = 0; l < N; l++)
	nmod_mat_entry(nbr_man->U,j,l) += scalar * nmod_mat_entry(nbr_man->Z,i,l);
    }

  // make sure that all the entries of U are between 0 and p^2
  for (i = 0; i < N - 2*nbr_man->k; i++)
    for (j = 0; j < N; j++)
      nmod_mat_entry(nbr_man->U,i,j) = nmod_mat_entry(nbr_man->U,i,j) % (p*p);

#ifdef DEBUG_LEVEL_FULL
  printf("after setting <U,X+Z> = 0\n");
  printf("U = ");
  nmod_mat_print(nbr_man->U);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  // Verify that U is now orthogonal to X+Z.
  for (i = 0; i < N-2*nbr_man->k; i++)
    for (j = 0; j < N; j++)
      nmod_mat_entry(B,2*nbr_man->k+i,j) = nmod_mat_entry(nbr_man->U,i,j);

  // The Gram matrix on this basis.
  nmod_mat_gram(temp,B, nbr_man->quot_gram,p*p);

  // Verify all is well.
  for ( i = 0; i < 2*nbr_man->k; i++)
    for ( j = 2*nbr_man->k; j < N; j++)
      assert(nmod_mat_entry(temp,i,j) % (p*p) == 0);

  nmod_mat_clear(temp);
#endif // DEBUG

  nmod_mat_clear(gram);
  nmod_mat_clear(X_new);
  nmod_mat_clear(Z_new);  
  nmod_mat_clear(B);
  // fmpz_mat_clear(B)
  fq_nmod_mat_clear(u, nbr_man->GF);
  free(excluded);
  fq_nmod_mat_clear(z, nbr_man->GF);
  free(paired);
  fq_nmod_mat_clear(x, nbr_man->GF);
  fq_nmod_mat_clear(basis_t, nbr_man->GF);
  fq_nmod_mat_clear(basis, nbr_man->GF);

  return;
}

void nbr_data_update_skew_matrix(nbr_data_t nbr_man, slong row, slong col)
{
  bool done;
  fq_nmod_t one;

  fq_nmod_init(one, nbr_man->GF);
  fq_nmod_one(one, nbr_man->GF);
  
  do {
    // Flag for determining whether we are done updating
    //  the skew matrix.
    done = true;
     
    // Increment value of the (row,col) position.
    fq_nmod_add(fq_nmod_mat_entry(nbr_man->p_skew,row,col), fq_nmod_mat_entry(nbr_man->p_skew,row,col), one, nbr_man->GF);
      
    // Update the coefficient of the skew matrix reflected
    //  across the anti-diagonal.
    fq_nmod_set(fq_nmod_mat_entry(nbr_man->p_skew, nbr_man->k-1-col, nbr_man->k-1-row),
		fq_nmod_mat_entry(nbr_man->p_skew, row, col), nbr_man->GF);
    fq_nmod_neg(fq_nmod_mat_entry(nbr_man->p_skew, nbr_man->k-1-col, nbr_man->k-1-row),
		fq_nmod_mat_entry(nbr_man->p_skew, nbr_man->k-1-col, nbr_man->k-1-row), nbr_man->GF);
      
    // If we've rolled over, move on to the next position.
    if (fq_nmod_is_zero(fq_nmod_mat_entry(nbr_man->p_skew,row,col),nbr_man->GF)) {
      // The next column of our skew matrix.
      col++;
      // Are we at the end of the column?
      if (row+col == nbr_man->k-1) {
	// Yes. Move to the next row.
	row++;
	// And reset the column.
	col = 0;
      }
      // Indicate we should repeat another iteration.
      done = false;
    }
  } while ((!done) && (row+col != nbr_man->k-1));
  
  fq_nmod_clear(one, nbr_man->GF);
  return;
}

void nbr_data_update_skew_space(nbr_data_t nbr_man)
{
  slong i,j,l;
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

  p = fmpz_get_si(fq_nmod_ctx_prime(nbr_man->GF));
  // Update the skew space.
  for (i = 0; i < nbr_man->k; i++) {
    for (l = 0; l < N; l++)
      nmod_mat_entry(nbr_man->X_skew, i, l) = nmod_mat_entry(nbr_man->X, i, l);
    for (j = 0; j < nbr_man->k; j++) {
      val = nmod_poly_get_coeff_ui(fq_nmod_mat_entry(nbr_man->p_skew,i,j),0);
      for (l = 0; l < N; l++)
	nmod_mat_entry(nbr_man->X_skew,i,l) += p * val * nmod_mat_entry(nbr_man->Z,j,l);
    }
    for (l = 0; l < N; l++)
      nmod_mat_entry(nbr_man->X_skew, i, j) %= (p*p);
  }

#ifdef DEBUG_LEVEL_FULL
  printf("X_skew = ");
  nmod_mat_print(nbr_man->X_skew);
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG
  fmpz_mat_init(B, N, N);
  fmpz_mat_init(B_t, N, N);
  fmpz_mat_init(Bq, N, N);
  fmpz_mat_init(temp, N, N);

  // Verify that X_skew is isotropic modulo p^2.
  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
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
    nbr_data_update_skew_matrix(nbr_man, row, col);
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
    nmod_mat_init_set(nbr_man->X_skew, nbr_man->X);
  }
  else {
    nmod_mat_clear(nbr_man->X);
    nmod_mat_clear(nbr_man->Z);
    nmod_mat_clear(nbr_man->U);
  }
  
  return;
}

void fmpz_quad_transform(fmpz_mat_t new_gram, const fmpz_mat_t gram, const fmpz_mat_t isom, slong scale)
{
  fmpz_mat_t isom_gram, isom_t;

  fmpz_mat_init(isom_gram, N, N);
  fmpz_mat_init(isom_t, N, N);
  
  fmpz_mat_transpose(isom_t, isom);
  fmpz_mat_mul(isom_gram, isom, gram);
  fmpz_mat_mul(new_gram, isom_gram, isom_t);

  fmpz_mat_scalar_divexact_si(new_gram, new_gram, scale*scale);
  
  fmpz_mat_clear(isom_gram);
  fmpz_mat_clear(isom_t);
  
  return;
}
// !! TODO - replace by the already implemented hnf in flint
void fmpz_mat_hermite_form(fmpz_mat_t H, const fmpz_mat_t s, slong d)
{
  slong row, pivot, col, j;
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

  fmpz_mat_init(g_h_prime, 1, N);
  fmpz_mat_init(b_primes, 1, N);
  
  fmpz_mat_one(H);
  fmpz_mat_scalar_mul_si(H, H, d);
  for (row = 0; row < N; row++) {
    for (col = 0; col < N; col++) {
      fmpz_set(fmpz_mat_entry(b_primes,0,col), fmpz_mat_entry(s, row, col));
    }
    // Here we compute the HNF of H and this row and store it in H
    for (pivot = 0; pivot < N; pivot++) {
      fmpz_set(a,fmpz_mat_entry(H,pivot,pivot));
      fmpz_set(b,fmpz_mat_entry(b_primes,0,pivot));
      fmpz_xgcd(g, x, y, a, b);
      fmpz_divexact(q_a, a, g);
      fmpz_divexact(q_b, b, g);
      for (col = 0; col < N; col++) {
	fmpz_mul(prod1, x, fmpz_mat_entry(H, pivot, col));
	fmpz_mul(prod2, y, fmpz_mat_entry(b_primes, 0, col));
	fmpz_add(fmpz_mat_entry(g_h_prime,0,col),prod1,prod2);
      }
      for (col = 0; col < N; col++) {
	fmpz_mul(prod1, q_a, fmpz_mat_entry(b_primes, 0, col));
	fmpz_mul(prod2, q_b, fmpz_mat_entry(H, pivot, col));
	fmpz_sub(fmpz_mat_entry(b_primes,0,col),prod1,prod2);
      }
      for (j = pivot; j < N; j++) {
	// !! TODO - check if this should be divexact
	fmpz_tdiv_q(scalar, fmpz_mat_entry(b_primes, 0, j), fmpz_mat_entry(H, j, j));
	for (col = 0; col < N; col++) {
	  fmpz_mul(prod1, scalar, fmpz_mat_entry(H,j,col));
	  fmpz_sub(fmpz_mat_entry(b_primes,0,col), fmpz_mat_entry(b_primes,0,col), prod1);
	}
      }
      for (col = 0; col < N; col++)
	fmpz_set(fmpz_mat_entry(H, pivot, col), fmpz_mat_entry(g_h_prime,0,col));
    }
    for (pivot = N-1; pivot > 0; pivot--) {
      for (col = pivot; col < N; col++) {
	fmpz_tdiv_q(q, fmpz_mat_entry(H, pivot-1, col), fmpz_mat_entry(H, col, col));
	for (j = col; j < N; j++) {
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
  slong i, j;
  fmpz_mat_t hermite;

  p = fmpz_get_si(fq_nmod_ctx_prime(nbr_man->GF));
  p2 = p*p;
  p3 = p2*p;

  // fill the isomtery by the X,Z,U
  // if lift_subspace was successful,
  // <X,X>, <X,U>,<Z,Z>,<Z,U> in p^2 and <X,Z> = 1 mod p^2

  // For now, we follow thw magma implementation
  // start by scaling the basis

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
      fmpz_set_si(fmpz_mat_entry(s,i,j),nmod_mat_entry(nbr_man->X_skew,i,j));

  for (i = 0; i < nbr_man->k; i++)
    for (j = 0; j < N; j++)
      fmpz_set_si(fmpz_mat_entry(s,nbr_man->k+i,j), p2 * nmod_mat_entry(nbr_man->Z,i,j));

  for (i = 0; i < N - 2  * nbr_man->k; i++)
    for (j = 0; j < N; j++)
      fmpz_set_si(fmpz_mat_entry(s,2*nbr_man->k+i,j), p * nmod_mat_entry(nbr_man->U,i,j));

  fmpz_mat_init(hermite, N, N);

  // !! TODO - replace by the already implemented hnf in flint
  fmpz_mat_hermite_form(hermite, s, p3);

  fmpz_mat_transpose(s, hermite);

  // need to adjust determinant for s to be in SO
  // This transforms using the isometry (and rescales already)
  fmpz_quad_transform(nbr, nbr_man->q, s, p);

  fmpz_mat_clear(hermite);
  
  return;
}
