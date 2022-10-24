#ifndef __NEIGHBOR_H__
#define __NEIGHBOR_H__

#include <carat/matrix.h>

#include "isometry.h"
#include "square_matrix.h"
#include "typedefs.h"

typedef struct {
  vector_t iso_vec;
  vector_t v;
  vector_t w;
  vector_t b;
  square_matrix_t Q;
  Z64 iso_j;
  int i; // the i from Gonzalo's code
  Z64 p;
  bool first_iter;
} neighbor_manager;

typedef neighbor_manager neighbor_manager_t[1];

void nbr_process_init(neighbor_manager_t nbr_man, const square_matrix_t Q, Z64 p, int i);

void nbr_process_advance(neighbor_manager_t nbr_man);

void nbr_process_clear(neighbor_manager_t nbr_man);

/* Compute one p-neighbour for Q_orig corresponding to vector x 
*/
void nbr_process_build_nb(square_matrix_t Q, const neighbor_manager_t nbr_man);

void nbr_process_build_nb_and_isom(square_matrix_t Q, isometry_t s, const neighbor_manager_t nbr_man);

/* find an isotropic vector for Q mod p */
/* return a row vector 1x5 */
bool get_isotropic_vector(vector_t x, const square_matrix_t Q, Z64 p);

/* get isotropic vector , correposnding to the pivot vector w. */
bool get_next_isotropic_vector(neighbor_manager_t nbr_man);

/* update the pivot vector v */
void update_pivot(vector_t v, Z64 p, int i);

/* get the next isotropic vector */
void update_isotropic_vector(square_matrix_t Q, Z64 p, vector_t v);

bool nbr_process_has_ended(const neighbor_manager_t nbr_man);

#endif // __NEIGHBOR_H__
