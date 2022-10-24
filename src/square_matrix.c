#include <carat/matrix.h>

#include "square_matrix.h"

void square_matrix_init(square_matrix_t mat)
{
  // no need to allocate any thing here
  return;
}

void square_matrix_clear(square_matrix_t mat)
{
  // we didn't allocate anything
  return;
}

void square_matrix_set(square_matrix_t dest, const square_matrix_t src)
{
  int i,j;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      dest[i][j] = src[i][j];

  return;
}

void square_matrix_set_fmpz_mat(square_matrix_t dest, const fmpz_mat_t src)
{
  int i,j;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      dest[i][j] = fmpz_get_si(fmpz_mat_entry(src,i,j));

  return;
}

int square_matrix_set_matrix_TYP(square_matrix_t dest, matrix_TYP* src)
{
  int i,j;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      dest[i][j] = src->array.SZ[i][j];

  return src->kgv;
}

void fmpz_mat_init_set_square_matrix(fmpz_mat_t dest, const square_matrix_t src)
{
  int i,j;
  fmpz_mat_init(dest, N, N);
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      fmpz_init_set_si(fmpz_mat_entry(dest,i,j), src[i][j]);

  return;
}

bool square_matrix_is_equal(const square_matrix_t mat1, const square_matrix_t mat2)
{
  int i,j;
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      if (mat1[i][j] != mat2[i][j])
	return false;

  return true;
}

void square_matrix_zero(square_matrix_t mat)
{
  int i,j;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      mat[i][j] = 0;

  return;
}

void square_matrix_one(square_matrix_t mat)
{
  int i;

  square_matrix_zero(mat);
  for (i = 0; i < N; i++)
    mat[i][i] = 1;

  return;
}

bool square_matrix_is_one(const square_matrix_t mat)
{
  int i,j;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      if (mat[i][j] != ((i == j) ? 1 : 0))
	return false;

  return true;
}

bool square_matrix_is_zero(const square_matrix_t mat)
{
  int i,j;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      if (mat[i][j] != 0)
	return false;

  return true;
}

// !! TODO - we could do this without copying by having a view
void square_matrix_transpose(square_matrix_t tr, const square_matrix_t mat)
{
  int i,j;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      tr[j][i] = mat[i][j];

  return;
}

// !! TODO - this could be made faster
void square_matrix_mul(square_matrix_t prod, const square_matrix_t matL, const square_matrix_t matR)
{
  int i,j,k;

  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) {
      prod[i][j] = 0;
      for (k = 0; k < N; k++)
	prod[i][j] += matL[i][k]*matR[k][j];
    }

  return;
}

int square_matrix_inv(square_matrix_t inv, const square_matrix_t mat, int denom)
{
  // !! TODO - at the moment we just use matrix_TYP*, see if we can do better
  matrix_TYP* s, s_inv;
  int inv_denom, i, j;

  s = init_mat(N,N,"");
  for (i = 0 ; i < N; i++)
    for (j = 0; j < N; j++)
      s->array.SZ[i][j] = mat[i][j];
  s->kgv = denom;

  s_inv = mat_inv(s);
  for (i = 0 ; i < N; i++)
    for (j = 0; j < N; j++)
      inv[i][j] = s_inv->array.SZ[i][j];

  inv_denom = s_inv->kgv;
  
  free_mat(s);
  free_mat(s_inv);

  return inv_denom;
}

void square_matrix_div_scalar(square_matrix_t quo, const square_matrix_t mat, int denom)
{
  int i,j;
  
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++) {
      assert(quo[i][j] % denom == 0);
      quo[i][j] = mat[i][j] / denom;
    }

  return;
}

void square_matrix_print(const square_matrix_t mat)
{
  int i,j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
      printf("%4ld ", mat[i][j]);
    printf("\n");
  }

  return;
}
