#include <assert.h>
#include "pivot_data.h"

void pivot_data_init(pivot_data_t pivots, slong dim, slong aniso, slong k)
{
  pivot_data_t pivots2;
  slong i,j;

  assert(k > 0);
  
  pivots->is_params_init = false;
  
  // Base case.
  if (k == 1) {
    pivots->total_len = dim - aniso;
    pivots->pivots = (slong**)malloc((dim-aniso)*sizeof(slong*));
    pivots->pivot_lens = (slong*)malloc((dim-aniso)*sizeof(slong));
    for (i = 0; i < dim - aniso; i++) {
      pivots->pivot_lens[i] = 1;
      pivots->pivots[i] = (slong*)malloc(sizeof(slong*));
      pivots->pivots[i][0] = i;
    }
    return;
  }

  // Retrieve lower-dimensional maximal pivots.
  pivot_data_init(pivots, dim-2, aniso, k-1);
  for (i = 0; i < pivots->total_len; i++)
    for (j = 0; j < pivots->pivot_lens[i]; j++)
      pivots->pivots[i][j]++;

  pivots->pivots = (slong**)realloc(pivots->pivots, 2*(pivots->total_len)*sizeof(slong*));
  pivots->pivot_lens = (slong*)realloc(pivots->pivot_lens, 2*(pivots->total_len)*sizeof(slong));
  for (i = 0; i < pivots->total_len; i++) {
    // we're already setting it to be longer, to save a realloc in the next for loop.
    pivots->pivots[i + pivots->total_len] = (slong*)malloc((pivots->pivot_lens[i]+1)*sizeof(slong*));
    for (j = 0; j < pivots->pivot_lens[i]; j++)
      pivots->pivots[i + pivots->total_len][j] = pivots->pivots[i][j];
    pivots->pivot_lens[i + pivots->total_len] = pivots->pivot_lens[i] + 1;
  }

  // Determine the first set of pivots.
  for (i = 0; i < pivots->total_len; i++) {
    pivots->pivots[i] = (slong*)realloc(pivots->pivots[i],(pivots->pivot_lens[i]+1)*sizeof(slong));
    for (j = pivots->pivot_lens[i]; j > 0 ; j--)
      pivots->pivots[i][j] = pivots->pivots[i][j-1];
    pivots->pivots[i][0] = 0;
    pivots->pivots[i + pivots->total_len][pivots->pivot_lens[i]] = dim-aniso-1;
    pivots->pivot_lens[i]++;
  }
  pivots->total_len *= 2;

  // Add additional pivots when we're not in the maximal case.
  if (2*k <= dim - aniso) {
    pivot_data_init(pivots2, dim-2, aniso, k);
    for (i = 0; i < pivots2->total_len; i++)
      for (j = 0; j < pivots2->pivot_lens[i]; j++)
	pivots2->pivots[i][j]++;
    pivots->pivots = (slong**)realloc(pivots->pivots, (pivots->total_len + pivots2->total_len)*sizeof(slong*));
    pivots->pivot_lens = (slong*)realloc(pivots->pivot_lens, (pivots->total_len + pivots2->total_len)*sizeof(slong));
    for (i = pivots->total_len - 1; i >= 0; i--) {
      pivots->pivots[i + pivots2->total_len] = pivots->pivots[i];
      pivots->pivot_lens[i + pivots2->total_len] = pivots->pivot_lens[i];
    }
    for (i = 0; i < pivots2->total_len; i++) {
      pivots->pivots[i] = (slong*)malloc(pivots2->pivot_lens[i]*sizeof(slong));
      pivots->pivot_lens[i] = pivots2->pivot_lens[i];
      for (j = 0; j < pivots2->pivot_lens[i]; j++)
	pivots->pivots[i][j] = pivots2->pivots[i][j];
    }
    pivots->total_len += pivots2->total_len;
    pivot_data_clear(pivots2);
  }
  
  return;
}

void pivot_data_clear(pivot_data_t pivots)
{
  slong i;
  
  for (i = 0; i < pivots->total_len; i++) {
    free(pivots->pivots[i]);
  }
  free(pivots->pivot_lens);
  free(pivots->pivots);
  return;
}

void pivot_data_params_clear(pivot_data_t pivots)
{
  slong i;
  
  free(pivots->free_vars);
  for (i = 0; i < pivots->num_params; i++)
    fq_nmod_clear(pivots->params[i], pivots->R->fqctx);
  free(pivots->params);
  pivots->num_params = 0;
  pivots->num_free_vars = 0;
  fq_nmod_mpoly_mat_clear(pivots->p_isotropic_param, pivots->R);
  fq_nmod_mpoly_ctx_clear(pivots->R);
  pivots->is_params_init = false;

  return;
}
