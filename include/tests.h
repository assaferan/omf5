#ifndef __TESTS_H__
#define __TESTS_H__

#include "typedefs.h"

STATUS test_61();
STATUS test_69();
STATUS test_greedy_overflow();

STATUS compute_eigenvectors(const int* Q_coeffs, const char* inp_type);

STATUS compute_eigenvalues(const int* Q_coeffs, int form_idx, int p, const char* inp_type);

STATUS compute_eigenvalues_up_to(const int* Q_coeffs, int form_idx, int prec, const char* inp_type);

#endif // __TESTS_H__
