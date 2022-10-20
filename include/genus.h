#ifndef __GENUS_H__
#define __GENUS_H__

#include <carat/matrix.h>

#include "hash.h"
#include "isometry.h"
#include "spinor.h"

typedef struct
{
  hash_table_t genus_reps;
  spinor_t spinor;
  
  slong* dims;
  slong* conductors;
  slong** num_auts;
  slong** lut_positions;
  slong num_conductors;

  isometry_t* isoms; // isometries corresponding to the genus representatives
  // matrix_TYP** isoms; // isometries corresponding to the genus representatives
  // slong* isom_denoms; // the denominators for the isometries
  
} genus_struct;

typedef genus_struct genus_t[1];

/* compute the genus of a quadratic form */
void genus_init(genus_t genus, matrix_TYP* Q);

void genus_clear(genus_t genus);

#endif // __GENUS_H__
