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
  
#ifdef DEBUG
  fq_nmod_clear(value, nbr_man->GF);
#endif // DEBUG
  fmpz_clear(p);
  
  return;
}

void nbr_data_clear(nbr_data_t nbr_man)
{
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
