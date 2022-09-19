#include <assert.h>

#include "carat/matrix.h"

#include "flint/fq_nmod.h"
#include "flint/fq_nmod_mat.h"
#include "flint/fq_nmod_mpoly.h"

#include "fq_nmod_mat.h"
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

  fq_nmod_mpoly_ctx_init(nbr_man->p_q_std_ctx, N*(N+1)/2, ORD_DEGREVLEX, nbr_man->GF);
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
