#ifndef __TESTS_H__
#define __TESTS_H__

int test_61();
int test_69();

int compute_eigenvectors(int* Q_coeffs);

int compute_eigenvalues(int* Q_coeffs, int form_idx, int p);

int compute_eigenvalues_up_to(int* Q_coeffs, int form_idx, int prec);

#endif // __TESTS_H__
