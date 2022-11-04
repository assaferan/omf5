#ifndef __SQUARE_MATRIX_H__
#define __SQUARE_MATRIX_H__

#include <carat/typedef.h>
#include <flint/fmpz_mat.h>

#include "typedefs.h"

#define QF_RANK 5

/*
in case we want more copmlicates stuff here

typedef struct {
  Z64 array[QF_RANK][QF_RANK];
} square_matrix;

typedef square_matrix square_matrix_t[1];
*/

typedef Z64 vector_t[QF_RANK];

typedef Z64 square_matrix_t[QF_RANK][QF_RANK];

void square_matrix_init(square_matrix_t mat);
void square_matrix_clear(square_matrix_t mat);

void square_matrix_init_set_symm(square_matrix_t mat, const int* coeff_vec, const char* inp_type);

void square_matrix_set(square_matrix_t dest, const square_matrix_t src);
void vector_set(vector_t dest, const vector_t src);

void square_matrix_set_fmpz_mat(square_matrix_t dest, const fmpz_mat_t src);
int square_matrix_set_matrix_TYP(square_matrix_t dest, matrix_TYP* src);
void fmpz_mat_init_set_square_matrix(fmpz_mat_t dest, const square_matrix_t src);
matrix_TYP* matrix_TYP_init_set_square_matrix(const square_matrix_t mat);

int vector_cmp(const vector_t vL, const vector_t vR);
bool square_matrix_is_equal(const square_matrix_t mat1, const square_matrix_t mat2);

void square_matrix_one(square_matrix_t mat);
void square_matrix_zero(square_matrix_t mat);

void vector_zero(vector_t vec);

bool square_matrix_is_one(const square_matrix_t mat);
bool square_matrix_is_zero(const square_matrix_t mat);
bool square_matrix_is_positive_definite(const square_matrix_t mat);

bool vector_is_zero(const vector_t vec);

bool square_matrix_is_bad_prime(const square_matrix_t mat, Z64 p);

void square_matrix_transpose(square_matrix_t tr, const square_matrix_t mat);
void square_matrix_add(square_matrix_t sum, const square_matrix_t matL, const square_matrix_t matR);
void square_matrix_mul(square_matrix_t prod, const square_matrix_t matL, const square_matrix_t matR);
void square_matrix_muleq_right(square_matrix_t prod, const square_matrix_t matR);
void square_matrix_muleq_left(square_matrix_t prod, const square_matrix_t matL);

void square_matrix_mul_vec_left(vector_t prod, const vector_t vec, const square_matrix_t mat);

void square_matrix_mul_scalar(square_matrix_t prod, const square_matrix_t mat, int scalar);

Z64 scalar_product(const vector_t v1, const vector_t v2);

int square_matrix_inv(square_matrix_t inv, const square_matrix_t mat, int denom);
void square_matrix_div_scalar(square_matrix_t quo, const square_matrix_t mat, int denom);

void vector_lin_comb(vector_t res, const vector_t v, const vector_t w, Z64 a_v, Z64 a_w);
void vector_mod_p(vector_t v, Z64 p);
void normalize_mod_p(vector_t v, Z64 p);

/* resymmetrizing a matrix after working with upper triangular part */
void square_matrix_resymmetrize(square_matrix_t Q);

/* swapping, assumes they do not point to the same thing!!! */
void square_matrix_swap_elts(square_matrix_t Q, int row1, int col1, int row2, int col2);

void square_matrix_print(const square_matrix_t mat);
void vector_print(const vector_t vec);

#endif // __SQUARE_MATRIX_H__
