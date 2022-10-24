#ifndef __ISOMETRY_H__
#define __ISOMETRY_H__

#include <carat/matrix.h>
#include <flint/fmpz_mat.h>

#include "square_matrix.h"
#include "typedefs.h"

// This structure keeps all the isometry information (denominator and inverse) in one place

typedef struct
{
  square_matrix_t s;
  square_matrix_t s_inv;
  int denom;
  int inv_denom;
  
} isometry_struct;

typedef isometry_struct isometry_t[1];

void isometry_init(isometry_t isom);
void isometry_init_set(isometry_t dest, const isometry_t src);
void isometry_init_set_square_matrix(isometry_t isom, const square_matrix_t s, int denom);
void isometry_init_set_fmpz_mat(isometry_t, const fmpz_mat_t s, int denom);
void isometry_mul(isometry_t prod, const isometry_t s1, const isometry_t s2);
void isometry_mul_mat_left(isometry_t prod, const square_matrix_t s1, int s1_denom, const isometry_t s2);
void isometry_mul_mat_right(isometry_t prod, const isometry_t s1, const square_matrix_t s2, int s2_denom);
void isometry_inv(isometry_t inv, const isometry_t isom);
void isometry_transform_gram(square_matrix_t gtQg, isometry_t g, const square_matrix_t Q);
bool isometry_is_isom(isometry_t s, const square_matrix_t q1, const square_matrix_t q2);
void isometry_clear(isometry_t isom);

// bool is_isometry(matrix_TYP* s, matrix_TYP* q1, matrix_TYP* q2, int denom);

#endif // __ISOMETRY_H__
