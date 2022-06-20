#ifndef __MATRIX_TOOLS_H__
#define __MATRIX_TOOLS_H__

#include <antic/nf.h>
#include <antic/nf_elem.h>

#include "carat/matrix.h"

/* function to initialize a symmetric matrix from a vector of coefficients */
matrix_TYP* init_sym_matrix(int* coeff_vec);

/* function to print matrices nicely*/
void print_mat(const matrix_TYP* Q);

/* swapping, assumes they do not point to the same thing!!! */
int swap(int** Q, int row1, int col1, int row2, int col2);

/* resymmetrizing a matrix after working with upper triangular part */
int resymmetrize(int **Q);

bravais_TYP* automorphism_group(matrix_TYP* Q);

matrix_TYP* is_isometric(matrix_TYP* Q1, matrix_TYP* Q2);

matrix_TYP* minkowski_reduce(matrix_TYP* Q);

void greedy(matrix_TYP* gram, matrix_TYP* s, int n, int dim);

struct eigenvalues_t {
  nf_t* nfs;
  nf_elem_t* eigenvals;
  nf_elem_t** eigenvecs; 
  int num;
  int dim;
};

typedef struct eigenvalues_t eigenvalues;

eigenvalues* get_eigenvalues(matrix_TYP* mat);

void free_eigenvalues(eigenvalues* evs);

#endif // __MATRIX_TOOLS_H__
