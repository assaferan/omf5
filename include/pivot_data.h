#ifndef __PIVOT_DATA_H__
#define __PIVOT_DATA_H__

#include "flint/flint.h"

typedef struct
{
  slong** pivots;
  slong* pivot_lens;
  slong total_len;
  slong pivot_ptr;
  
} pivot_data;

typedef pivot_data pivot_data_t[1];

void pivot_data_init(pivot_data_t pivots, slong nsing_dim, slong aniso_dim, slong k);
void pivot_data_clear(pivot_data_t pivots);

#endif // __PIVOT_DATA_H__
