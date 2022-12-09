#include <assert.h>

#include <carat/autgrp.h>
#include <carat/symm.h>

#include "isometry.h"
#include "matrix_tools.h"

void isometry_init(isometry_t isom)
{
  // by default we initialize isometries to one
  square_matrix_one(isom->s);
  square_matrix_one(isom->s_inv);
  isom->denom = 1;
  isom->inv_denom = 1;

  return;
}

void isometry_init_set(isometry_t dest, const isometry_t src)
{
  square_matrix_set(dest->s, src->s);
  square_matrix_set(dest->s_inv, src->s_inv);
  dest->denom = src->denom;
  dest->inv_denom = src->inv_denom;

  return;
}

void isometry_init_set_square_matrix(isometry_t isom, const square_matrix_t s, int denom)
{
  // maybe we want to save time copying - check later
  square_matrix_set(isom->s,s);
  isom->denom = denom;
  isom->inv_denom = square_matrix_inv(isom->s_inv, s, isom->denom);

  return;
}

void isometry_init_set_fmpz_mat(isometry_t isom, const fmpz_mat_t s, int denom)
{

  square_matrix_set_fmpz_mat(isom->s, s);
  isom->denom = denom;
  isom->inv_denom = square_matrix_inv(isom->s_inv, isom->s, isom->denom);

  return;
}

void isometry_mul(isometry_t prod, const isometry_t s1, const isometry_t s2)
{

  square_matrix_mul(prod->s, s1->s, s2->s);
  square_matrix_mul(prod->s_inv, s2->s_inv, s1->s_inv);
  prod->denom = (s1->denom) * (s2->denom);
  prod->inv_denom = (s1->inv_denom) * (s2->inv_denom);

  return;
}

void isometry_muleq_right(isometry_t sL, const isometry_t sR)
{
  square_matrix_muleq_right(sL->s, sR->s);
  square_matrix_muleq_left(sL->s_inv, sR->s_inv);
  sL->denom = (sL->denom) * (sR->denom);
  sL->inv_denom = (sL->inv_denom) * (sR->inv_denom);

  return;
}

void isometry_muleq_left(isometry_t sR, const isometry_t sL)
{
  square_matrix_muleq_left(sR->s, sL->s);
  square_matrix_muleq_right(sR->s_inv, sL->s_inv);
  sR->denom = (sL->denom) * (sR->denom);
  sR->inv_denom = (sL->inv_denom) * (sR->inv_denom);

  return;
}

void isometry_mul_mat_left(isometry_t prod, const square_matrix_t s1, int s1_denom, const isometry_t s2)
{
  
  square_matrix_mul(prod->s, s1, s2->s);
  prod->denom = s1_denom * (s2->denom);
  prod->inv_denom = square_matrix_inv(prod->s_inv, prod->s, prod->denom);

  return;
}

void isometry_mul_mat_right(isometry_t prod, const isometry_t s1, const square_matrix_t s2, int s2_denom)
{
  
  square_matrix_mul(prod->s, s1->s, s2);
  prod->denom = s1->denom * (s2_denom);
  prod->inv_denom = square_matrix_inv(prod->s_inv, prod->s, prod->denom);
  
  return;
}

void isometry_mul_mateq_left(isometry_t sR, const square_matrix_t sL, int sL_denom)
{
  square_matrix_muleq_left(sR->s, sL);
  sR->denom = sL_denom * (sR->denom);
  sR->inv_denom = square_matrix_inv(sR->s_inv, sR->s, sR->denom);

  return;
}
void isometry_mul_mateq_right(isometry_t sL, const square_matrix_t sR, int sR_denom)
{
  square_matrix_muleq_right(sL->s, sR);
  sL->denom = (sL->denom) * (sR_denom);
  sL->inv_denom = square_matrix_inv(sL->s_inv, sL->s, sL->denom);

  return;
}

void isometry_inv(isometry_t inv, const isometry_t isom)
{
  
  square_matrix_set(inv->s, isom->s_inv);
  square_matrix_set(inv->s_inv, isom->s);
  inv->denom = isom->inv_denom;
  inv->inv_denom = isom->denom;

  return;
}

void isometry_clear(isometry_t isom)
{
  // did not allocate any memory
  return;
}

// !! TODO - we can make this operation faster by doing it all at once
void isometry_transform_gram(square_matrix_t gtQg, const isometry_t g, const square_matrix_t Q)
{
  int i,j,k,l;
  square_matrix_t /* gt, Qg, */ gtQg2;

  /*
  square_matrix_transpose(gt, g->s);
  square_matrix_mul(Qg, Q, g->s);
  */
  
  for (i = 0; i < QF_RANK; i++)
    for (j = i; j < QF_RANK; j++) {
      gtQg2[i][j] = 0;
      for (k = 0; k < QF_RANK; k++)
	for (l = 0; l < QF_RANK; l++)
	  gtQg2[i][j] += Q[k][l] * (g->s)[k][i] * (g->s)[l][j];
    }

  // square_matrix_mul(gtQg, gt, Qg);
  
  for (i = 0; i < QF_RANK; i++) {
    for (j = 0; j < i; j++)
      gtQg[i][j] = gtQg2[j][i];
    for (j = i; j < QF_RANK; j++)
      gtQg[i][j] = gtQg2[i][j];
  }

  //  assert(square_matrix_is_equal(gtQg, gtQg2));
  
  // this might overflow, as denom is only an int
  // square_matrix_div_scalar(gtQg, gtQg, (g->denom)*(g->denom));
  square_matrix_div_scalar(gtQg, gtQg, g->denom);
  square_matrix_div_scalar(gtQg, gtQg, g->denom);

  return;
}

bool isometry_is_isom(isometry_t s, const square_matrix_t q1, const square_matrix_t q2)
{
  square_matrix_t gram;

#ifdef DEBUG_LEVEL_FULL
  printf("s = \n");
  square_matrix_print(s->s);
  printf("\n");
  printf("q1 = \n");
  square_matrix_print(q1);
  printf("\n");
  printf("q2 = \n");
  square_matrix_print(q2);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  isometry_transform_gram(gram, s, q1);
  
  return square_matrix_is_equal(gram,q2);
}

bool is_isometric(isometry_t s, const square_matrix_t q1, const square_matrix_t q2)
{
  int i, Qmax;
  matrix_TYP *SV1, *SV2, *isom, *Q1, *Q2;
  square_matrix_t s_sm, s_tr;
  int options[6] = {0};

  assert(square_matrix_is_positive_definite(q1));
  assert(square_matrix_is_positive_definite(q2));
  
  Q1 = matrix_TYP_init_set_square_matrix(q1);
  Q2 = matrix_TYP_init_set_square_matrix(q2);
  
  Qmax = Q1->array.SZ[0][0];
  for(i=1;i < Q1->cols;i++) {
    if(Q1->array.SZ[i][i] > Qmax)
      Qmax = Q1->array.SZ[i][i];
  }
  SV1 = short_vectors(Q1, Qmax, 0, 0, 0, &i);
  SV2 = short_vectors(Q2, Qmax, 0, 0, 0, &i);

  isom = isometry(&Q1, &Q2, 1, SV1, SV2, NULL, 0, options);
  if (isom != NULL) {
    square_matrix_set_matrix_TYP(s_sm, isom);
    square_matrix_transpose(s_tr, s_sm);
    isometry_init_set_square_matrix(s, s_tr, 1);
    free_mat(isom);
  }

  free_mat(SV1);
  free_mat(SV2);
  free_mat(Q1);
  free_mat(Q2);

  return (isom != NULL);
}

void isometry_swap_vecs(isometry_t isom, int col1, int col2)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  int denom;
#endif // DEBUG
  int row;
  
  for (row = 0; row < QF_RANK; row++)
    square_matrix_swap_elts(isom->s, row, col1, row, col2);
  for (row = 0; row < QF_RANK; row++)
    square_matrix_swap_elts(isom->s_inv, col1, row, col2, row);

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG
  
  return;
}

void isometry_add_vec(isometry_t isom, int target_col, Z64 c, int src_col)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  int denom;
#endif // DEBUG
  int row;

  for (row = 0; row < QF_RANK; row++)
    isom->s[row][target_col] += c * isom->s[row][src_col];

  for (row = 0; row < QF_RANK; row++)
    isom->s_inv[src_col][row] -= c * isom->s_inv[target_col][row];

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG
  
  return;
}

void isometry_replace_vec(isometry_t isom, int col, const vector_t x)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  int denom;
#endif // DEBUG
  int row;
  
  for (row = 0; row < QF_RANK; row++)
    isom->s[row][col] = isom->s[row][0];
  for (row = 0; row < QF_RANK; row++)
    isom->s[row][0] = x[row];
  isom->inv_denom = isom->s[col][0];

  // and updating s_inv
  isom->s_inv[0][0] = 0;
  isom->s_inv[col][0] = 1;
  isom->s_inv[0][col] = 1;
  for (row = 1; row < QF_RANK; row++)
    isom->s_inv[row][col] = -isom->s[row][0];
  if (col != 0)
    isom->s_inv[col][col] = -isom->s[0][0];

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG
  
  return;
}

void isometry_vec_scalar_div(isometry_t isom, int col, int scalar)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  int denom;
#endif // DEBUG
  
  int i,j;

  isom->denom *= scalar;
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      if (j != col)
	isom->s[i][j] *= scalar;

  for (i = 0; i < QF_RANK; i++)
    isom->s_inv[col][i] *= scalar;

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG

  return;
}

void isometry_vec_scalar_mul(isometry_t isom, int col, int scalar)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  int denom;
#endif // DEBUG
  
  int i,j;

  isom->inv_denom *= scalar;
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      if (j != col)
	isom->s_inv[j][i] *= scalar;

  for (i = 0; i < QF_RANK; i++)
    isom->s[i][col] *= scalar;

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG

  return;
}

bool isometry_is_equal(const isometry_t isom1, const isometry_t isom2)
{
  square_matrix_t scaled_isom1, scaled_isom2;

  square_matrix_mul_scalar(scaled_isom1, isom1->s, isom2->denom);
  square_matrix_mul_scalar(scaled_isom2, isom2->s, isom1->denom);

  return square_matrix_is_equal(scaled_isom1, scaled_isom2);
}
