#ifndef __TESTS_H__
#define __TESTS_H__

#include "typedefs.h"

typedef struct {

  int Q_coeffs[15];
  const char* format;
  
  int num_conductors;
  int* dims;
  int* num_forms;

  int num_ps[2];
  int* ps[2];
  
  int**** test_evs;
  
} example_struct;

typedef example_struct example_t[1];

STATUS test_61();
STATUS test_69();
STATUS test_greedy_overflow();

STATUS compute_eigenvectors(const int* Q_coeffs, const char* inp_type);

STATUS compute_eigenvalues(const int* Q_coeffs, int p, const char* inp_type);

STATUS compute_eigenvalues_up_to(const int* Q_coeffs, int form_idx, int prec, const char* inp_type);

#endif // __TESTS_H__
