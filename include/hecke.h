#ifndef __HECKE_H__
#define __HECKE_H__

#include "carat/matrix.h"

#include "genus.h"

int process_isotropic_vector(matrix_TYP* v, matrix_TYP* w_mat, matrix_TYP* Q,
			     int p, matrix_TYP* b, int* T, int* th61);

int q61_nbs1(int* T, int p, int i, nbrs_data* init_orig);

#endif // __HECKE_H__
