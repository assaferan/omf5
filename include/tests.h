#ifndef __TESTS_H__
#define __TESTS_H__

typedef int STATUS;

#define SUCCESS 0
#define FAIL    1

STATUS test_61();
STATUS test_69();
STATUS test_greedy();

STATUS compute_eigenvectors(const int* Q_coeffs);

STATUS compute_eigenvalues(const int* Q_coeffs, int form_idx, int p);

STATUS compute_eigenvalues_up_to(const int* Q_coeffs, int form_idx, int prec);

#endif // __TESTS_H__
