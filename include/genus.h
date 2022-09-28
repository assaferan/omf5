#ifndef __GENUS_H__
#define __GENUS_H__

#include <carat/matrix.h>

#include "hash.h"

typedef struct
{
  hash_table_t genus;

  slong* dims;
  slong total_dim;
  
} omf_space;

typedef omf_space omf_space_t[1];

/* compute the genus of a quadratic form */
void get_genus_reps(hash_table_t table, matrix_TYP* Q);

#endif // __GENUS_H__
