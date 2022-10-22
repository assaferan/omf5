#ifndef __ISOMETRY_H__
#define __ISOMETRY_H__

#include <carat/matrix.h>
#include <flint/fmpz_mat.h>

#include "typedefs.h"

// This structure keeps all the isometry information (denominator and inverse) in one place

typedef struct
{
  
  matrix_TYP* s;
  matrix_TYP* s_inv;
  // This is not needed matrix_TYP already has kgv that stores the denominator
  // slong denom; 
  
} isometry_struct;

typedef isometry_struct isometry_t[1];

void isometry_init(isometry_t isom, matrix_TYP* s, slong denom);
void isometry_init_fmpz_mat(isometry_t, fmpz_mat_t s, slong denom);
void isometry_mul(isometry_t prod, isometry_t s1, isometry_t s2);
void isometry_mul_mat_left(isometry_t prod, matrix_TYP* s1, isometry_t s2);
void isometry_mul_mat_right(isometry_t prod, isometry_t s1, matrix_TYP* s2);
void isometry_inv(isometry_t inv, isometry_t isom);
matrix_TYP* isometry_transform_gram(isometry_t g, matrix_TYP* Q);
bool isometry_is_isom(isometry_t s, matrix_TYP* q1, matrix_TYP* q2);
void isometry_clear(isometry_t isom);

bool is_isometry(matrix_TYP* s, matrix_TYP* q1, matrix_TYP* q2, int denom);

#endif // __ISOMETRY_H__
