#ifndef __TESTS_H__
#define __TESTS_H__

#include "genus.h"
#include "typedefs.h"

typedef struct {

  Z64 Q_coeffs[15];
  const char* format;
  
  int num_conductors;
  int* dims;
  
} example_genus;

typedef example_genus example_genus_t[1];

typedef struct {

  int* num_forms;

  int num_ps[2];
  int* ps[2];
  
  int**** test_evs;
  
} example_evs;

typedef example_evs example_evs_t[1];

STATUS test_61();
STATUS test_69();
STATUS test_greedy_overflow();

STATUS compute_eigenvectors(const genus_t genus, bool nonlifts, int num_idxs, const int* idxs);

STATUS compute_eigenvalues(const genus_t genus, int p, bool nonlifts, int num_idxs, const int* idxs);

STATUS compute_eigenvalues_up_to(const genus_t genus, int prec, bool nonlifts, int num_idxs, const int* idxs);

STATUS compute_hecke_col(const genus_t genus, slong p, slong c);

STATUS compute_hecke_col_all_conds(const genus_t genus, slong p, int gen_idx);

STATUS compute_hecke_matrix(const genus_t genus, slong p, slong c);

STATUS compute_hecke_matrix_all_conds(const genus_t genus, slong p);

STATUS compute_first_hecke_matrix_all_conds(const genus_t genus);

void compute_genus(genus_t genus, const Z64* Q_coeffs, const char* format);

STATUS compute_first_large_hecke(const genus_t genus);

#endif // __TESTS_H__
