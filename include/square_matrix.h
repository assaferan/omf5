#ifndef __SQUARE_MATRIX_H__
#define __SQUARE_MATRIX_H__

#include "typedefs.h"

#define N 5

/*
in case we want more copmlicates stuff here

typedef struct {
  Z64 array[N][N];
} square_matrix;

typedef square_matrix square_matrix_t[1];
*/

typedef Z64 square_matrix_t[N][N];

void square_matrix_init(square_matrix_t mat);
void square_matrix_clear(square_matrix_t mat);

void square_matrix_set(square_matrix_t dest, const square_matrix_t src);

void square_matrix_set_fmpz_mat(square_matrix_t dest, const fmpz_mat_t src);
int square_matrix_set_matrix_TYP(square_matrix_t dest, matrix_TYP* src);
void fmpz_mat_init_set_square_matrix(fmpz_mat_t dest, const square_matrix_t src);

bool square_matrix_is_equal(const square_matrix_t mat1, const square_matrix_t mat2);

void square_matrix_one(square_matrix_t mat);
void square_matrix_zero(square_matrix_t mat);

bool square_matrix_is_one(const square_matrix_t mat);
bool square_matrix_is_zero(const square_matrix_t mat);

void square_matrix_transpose(square_matrix_t tr, const square_matrix_t mat);
void square_matrix_mul(square_matrix_t prod, const square_matrix_t matL, const square_matrix_t matR);
int square_matrix_inv(square_matrix_t inv, const square_matrix_t mat, int denom);
void square_matrix_div_scalar(square_matrix_t quo, const square_matrix_t mat, int denom);

void square_matrix_print(const square_matrix_t mat);

#endif // __SQUARE_MATRIX_H__
