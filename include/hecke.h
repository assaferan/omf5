#ifndef __HECKE_H__
#define __HECKE_H__

#include "carat/matrix.h"

#include "genus.h"
#include "matrix_tools.h"
#include "neighbor.h"

int process_isotropic_vector(neighbor_manager* nbr_man, int* T, hash_table* genus, double* theta_time, double* isom_time, double* total_time, int* num_isom);

int process_neighbour_chunk(int* T, int p, int i, int gen_idx, hash_table* genus);

void hecke_col(int* T, int p, int gen_idx, hash_table* genus);

void get_hecke_ev(nf_elem_t e, hash_table* genus, eigenvalues* evs, int p, int ev_idx);

matrix_TYP* hecke_matrix(hash_table* genus, int p);

#endif // __HECKE_H__
