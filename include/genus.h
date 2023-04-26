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
void genus_init_square_matrix(genus_t genus, const square_matrix_t Q, int h);

/* set the genus from a list */
void genus_init_set_square_matrix_vec(genus_t genus, const square_matrix_t* reps, size_t genus_size);
void genus_init_set_square_matrix_vec_and_isoms(genus_t genus, const square_matrix_t* reps,
						const isometry_t* isoms, size_t genus_size);

void genus_init_empty(genus_t genus, size_t disc);

void genus_init_file(genus_t genus, const char* genus_fname, size_t disc, bool has_isoms);

void genus_clear(genus_t genus);

bool genus_is_trivial_cond(const genus_t genus);

#endif // __GENUS_H__
