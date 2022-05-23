#ifndef __NEIGHBOR_H__
#define __NEIGHBOR_H__

#include "carat/matrix.h"

struct neighbor_manager_t {
  matrix_TYP* iso_vec;
  matrix_TYP* v;
  matrix_TYP* w;
  matrix_TYP* b;
  matrix_TYP* Q;
  int iso_j;
  int i; // the i from Gonzalo's code
  int p;
};

typedef struct neighbor_manager_t neighbor_manager;

void init_nbr_process(neighbor_manager* nbr_man, matrix_TYP* Q, int p, int i);

void advance_nbr_process(neighbor_manager* nbr_man);

void free_nbr_process(neighbor_manager* nbr_man);

/* Compute one p-neighbour for Q_orig corresponding to vector x 
 * On error, return NULL.
*/
matrix_TYP* q61_nb(matrix_TYP* Q_orig, int p, matrix_TYP* x_mat);

/* find an isotropic vector for Q mod p */
/* return a row vector 1x5 */
matrix_TYP* get_isotropic_vector(matrix_TYP* Q, int p);

/* get isotropic vector , correposnding to the pivot vector w. */

matrix_TYP* get_next_isotropic_vector(neighbor_manager* nbr_man);

/* update the pivot vector v */
void update_pivot(int* v, int p, int i);

/* get the next isotropic vector */
void update_isotropic_vector(matrix_TYP*Q, int p, matrix_TYP* v);

int has_ended(neighbor_manager* nbr_man);

#endif // __NEIGHBOR_H__
