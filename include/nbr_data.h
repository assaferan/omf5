#ifndef __NBR_DATA_H__
#define __NBR_DATA_H__

#include "carat/matrix.h"

#include "flint/fq_nmod.h"
#include "flint/fq_nmod_vec.h"
#include "flint/fq_nmod_mat.h"
#include "flint/fq_nmod_mpoly.h"

#include "pivot_data.h"
#include "typedefs.h"

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
  bool is_iso_subspace_init;
  
} nbr_data;

typedef nbr_data nbr_data_t[1];

void nbr_data_init(nbr_data_t nbr_man, matrix_TYP* q, slong p_int, slong k);

void nbr_data_clear(nbr_data_t nbr_man);
void nbr_data_params_init(pivot_data_t pivots, const nbr_data_t nbr_man);

void nbr_data_next_isotropic_subspace(nbr_data_t nbr_man);

void nbr_data_lift_subspace(nbr_data_t nbr_man);

void nbr_data_build_neighbor(fmpz_mat_t nbr, fmpz_mat_t s, const nbr_data_t nbr_man);

void nbr_data_get_next_neighbor(nbr_data_t nbr_man);

bool nbr_data_has_ended(const nbr_data_t nbr_man);

#endif // __NBR_DATA_H__
