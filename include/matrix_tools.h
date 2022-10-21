#ifndef __MATRIX_TOOLS_H__
#define __MATRIX_TOOLS_H__

#include <antic/nf.h>
#include <antic/nf_elem.h>

#include "carat/matrix.h"

#include "flint/fmpz_mat.h"
#include "flint/fmpq_mat.h"
#include "flint/fmpq_poly.h"
#include "flint/fq_nmod.h"
#include "flint/fq_nmod_mat.h"

/* function to initialize a symmetric matrix from a vector of coefficients */
matrix_TYP* init_sym_matrix(const int* coeff_vec, const char* inp_type);

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

void eigenvalues_init(eigenvalues** evs, slong num, slong dim);

bool get_eigenvector_on_subspace(nf_elem_t* evec, nf_t nf, const fmpq_mat_t T, const fmpq_mat_t basis_W);

eigenvalues* get_eigenvalues(matrix_TYP* mat);

void free_eigenvalues(eigenvalues* evs);

void fmpq_mat_init_set_matrix_TYP(fmpq_mat_t M, const matrix_TYP* mat);
void fmpz_mat_init_set_matrix_TYP(fmpz_mat_t M, const matrix_TYP* mat);
void matrix_TYP_init_set_fmpz_mat(matrix_TYP** mat, const fmpz_mat_t M);
void nmod_mat_init_set_fmpz_mat(nmod_mat_t dest, const fmpz_mat_t mat, mp_limb_t n);
void fq_nmod_mat_init_set_fmpz_mat(fq_nmod_mat_t dest, const fmpz_mat_t mat, const fq_nmod_ctx_t F);

void restrict_mat(fmpq_mat_t res_T, const fmpq_mat_t T, const fmpq_mat_t basis_W);

void kernel_on(fmpq_mat_t ker, const fmpq_mat_t A, const fmpq_mat_t B);

#endif // __MATRIX_TOOLS_H__
