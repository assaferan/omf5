#ifndef __TESTS_H__
#define __TESTS_H__

#include "genus.h"
#include "typedefs.h"

typedef struct {

  int Q_coeffs[15];
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

STATUS compute_eigenvectors(const genus_t genus);

STATUS compute_eigenvalues(const genus_t genus, int p);

STATUS compute_eigenvalues_up_to(const genus_t genus, int prec);

STATUS compute_hecke_col(const genus_t genus, int p, int c);

STATUS compute_hecke_col_all_conds(const genus_t genus, int p, int gen_idx);

STATUS compute_first_hecke_matrix_all_conds(const genus_t genus);

void compute_genus(genus_t genus, const int* Q_coeffs, const char* format);

#endif // __TESTS_H__
