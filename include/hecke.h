#ifndef __HECKE_H__
#define __HECKE_H__

#include "carat/matrix.h"

#include "genus.h"
#include "neighbor.h"

/*
int process_isotropic_vector(matrix_TYP* v, matrix_TYP* w_mat, matrix_TYP* Q,
			     int p, matrix_TYP* b, int* T, int* th61);
*/

int process_isotropic_vector(neighbor_manager* nbr_man, int* T,
			     int* th61, hash_table* genus);

int q61_nbs1(int* T, int p, int i, nbrs_data* init_orig, hash_table* genus);

#endif // __HECKE_H__
