#ifndef __NBR_DATA_H__
#define __NBR_DATA_H__

#include "carat/matrix.h"

#include "flint/fq_nmod.h"
#include "flint/fq_nmod_vec.h"
#include "flint/fq_nmod_mat.h"
#include "flint/fq_nmod_mpoly.h"

#include "pivot_data.h"

typedef struct
{
  fq_nmod_ctx_t GF; // GF(p) the finite field for the reduction
  fq_nmod_mat_t vec; // an isotropic vector mod p
  fq_nmod_mat_t b; // reduction of the bilinear form mod p
  nmod_mat_t quot_gram; // reduction modulo p^2

  fmpz_mat_t q; // the quadratic form
  fmpz_t disc;

  fq_nmod_mat_t p_std_gram;
  fq_nmod_mat_t p_basis;
  fq_nmod_mpoly_t p_q_std;
  fq_nmod_mpoly_ctx_t p_q_std_ctx; // polynomial ring for p_q_std

  // dimension of the radical
  slong rad_dim;
  // dimension of the anisotropic subspace
  slong aniso_dim;
  // the Witt index (number of hyperbolic planes)
  slong witt_index;
  
  pivot_data_t pivots;
  slong k; // dimension of the isotropic subspace - in our case it will always be 2
  slong skew_dim;
  fq_nmod_mat_t p_skew;

  fq_nmod_mat_t iso_subspace;
  nmod_mat_t X, Z, U, X_skew;
  // fmpz_mat_t X, Z, U, X_skew;
  bool is_done;
  bool is_skew_init;
  
} nbr_data;

typedef nbr_data nbr_data_t[1];

void nbr_data_init(nbr_data_t nbr_man, matrix_TYP* q, slong p_int, slong k);

void nbr_data_clear(nbr_data_t nbr_man);
void nbr_data_params_init(pivot_data_t pivots, const nbr_data_t nbr_man);

void nbr_data_next_isotropic_subspace(nbr_data_t nbr_man);

void nbr_data_lift_subspace(nbr_data_t nbr_man);

/* void advance_nbr2_process(neighbor2_manager* nbr_man); */

/* /\* Compute one p-neighbour for Q_orig corresponding to vector x  */
/*  * On error, return NULL. */
/* *\/ */
/* matrix_TYP* build_nb2(neighbor2_manager* nbr_man); */

/* /\* get isotropic subspace, correposnding to the pivot matrix w. *\/ */

/* matrix_TYP* get_next_isotropic_subspace(neighbor2_manager* nbr_man); */

/* /\* update the pivot vector v *\/ */
/* void update_pivot(int* v, int p, int i); */

/* /\* get the next isotropic vector *\/ */
/* void update_isotropic_subspace(matrix_TYP*Q, int p, matrix_TYP* v); */

/* int has_ended(neighbor2_manager* nbr_man); */

#endif // __NBR_DATA_H__
