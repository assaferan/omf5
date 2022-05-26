#ifndef __MATRIX_TOOLS_H__
#define __MATRIX_TOOLS_H__

#include "carat/matrix.h"

/* function to initialize a symmetric matrix from a vector of coefficients */
matrix_TYP* init_sym_matrix(int* coeff_vec);

/* function to print matrices nicely*/
void print_mat(matrix_TYP* Q);

/* swapping, assumes they do not point to the same thing!!! */
int swap(int** Q, int row1, int col1, int row2, int col2);

/* resymmetrizing a matrix after working with upper triangular part */
int resymmetrize(int **Q);

bravais_TYP* automorphism_group(matrix_TYP* Q);

matrix_TYP* is_isometric(matrix_TYP* Q1, matrix_TYP* Q2);

matrix_TYP* minkowski_reduce(matrix_TYP* Q);

void greedy(matrix_TYP* gram, matrix_TYP* s, int n, int dim);

#endif // __MATRIX_TOOLS_H__
