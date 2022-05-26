#ifndef __HECKE_H__
#define __HECKE_H__

#include "carat/matrix.h"

#include "genus.h"
#include "neighbor.h"

int process_isotropic_vector(neighbor_manager* nbr_man, int* T, hash_table* genus);

int process_neighbour_chunk(int* T, int p, int i, int gen_idx, hash_table* genus);

void hecke_col(int* T, int p, int gen_idx, hash_table* genus);

matrix_TYP* hecke_matrix(hash_table* genus, int p);

#endif // __HECKE_H__
