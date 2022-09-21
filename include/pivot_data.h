#ifndef __PIVOT_DATA_H__
#define __PIVOT_DATA_H__

#include "flint/flint.h"

#include "fq_nmod_mpoly_mat.h"

typedef struct
{
  slong** pivots;
  slong* pivot_lens;
  slong total_len;
  slong pivot_ptr;

  slong* free_vars;
  fq_nmod_t* params;

  slong num_free_vars;
  slong num_params;
  
  fq_nmod_mpoly_mat_t p_isotropic_param;

  fq_nmod_mpoly_ctx_t R; // polynomial ring for p_isotropic_param
  
} pivot_data;

typedef pivot_data pivot_data_t[1];

void pivot_data_init(pivot_data_t pivots, slong dim, slong aniso_dim, slong k);
void pivot_data_clear(pivot_data_t pivots);
void pivot_data_params_clear(pivot_data_t pivots);

#endif // __PIVOT_DATA_H__
