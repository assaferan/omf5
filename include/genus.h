#ifndef __GENUS_H__
#define __GENUS_H__

#include "carat/matrix.h"

#include "hash.h"

/* a struct that saves the hash table (theta values),
 * as well as the isotropic vector and the form we work with
 */

struct nbrs_data_t {
  matrix_TYP* Q;
  matrix_TYP* v;
  int* th61;
};

typedef struct nbrs_data_t nbrs_data;

/* initialize the neighbors data */
void init_nbrs_data(nbrs_data* dtm, int p, int k);

/* identify the genus representative of Q (returns the index) */
int q61_id(matrix_TYP* Q, int* th61);

/* compute the genus of a quadratic form */
hash_table* get_genus_reps(matrix_TYP* Q);

void free_nbrs_data(nbrs_data* dtm);

#endif // __GENUS_H__
