#ifndef __GENUS_H__
#define __GENUS_H__

#include <carat/matrix.h>

#include "hash.h"
#include "isometry.h"
#include "spinor.h"
#include "square_matrix.h"

typedef struct
{
  hash_table_t genus_reps;
  spinor_t spinor;
  fmpz_t disc;
  
  slong* dims;
  slong* conductors;
  slong** num_auts;
  slong** lut_positions;
  slong num_conductors;

  isometry_t* isoms; // isometries corresponding to the genus representatives
  
} genus_struct;

typedef genus_struct genus_t[1];

/* compute the genus of a quadratic form */
void genus_init_square_matrix(genus_t genus, const square_matrix_t Q);

/* set the genus from a list */
void genus_init_set_square_matrix_vec(genus_t genus, const square_matrix_t* reps, int genus_size);

void genus_clear(genus_t genus);

#endif // __GENUS_H__
