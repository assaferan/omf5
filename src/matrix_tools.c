#include "carat/autgrp.h"
#include "carat/reduction.h"
#include "carat/symm.h"

#include "matrix_tools.h"

/* Right now working with matrices as in CARAT package
 * using matrix_TYP.
 * May consider using instead fixed size matrices
 * (without heap allocation, always using a constant size) 
 * may also consider other implementations such as flint or eigen
*/

/* For now, we use N = 5 hard coded */

#define N 5

/* function to initialize a symmetric matrix */
matrix_TYP* init_sym_matrix(int* coeff_vec)
{
  int row, col, idx;
  matrix_TYP* Q_mat;
  int** Q;
  
  Q_mat = init_mat(N,N,"");
  Q = Q_mat->array.SZ;
  row = 0;
  col = 0;
  for (idx = 0; idx < N*(N+1)/2; idx++)
    {
      Q[row][col] = coeff_vec[idx];
      row++;
      if (row > col) {
	row = 0;
	col++;
      }
    }
  for (row = 0; row < N-1; row++)
    for (col = row+1; col < N; col++)
      Q[col][row] = Q[row][col];

  return Q_mat;
}

/* function to print matrices nicely
 */

void print_mat(matrix_TYP* Q)
{
  int row, col;
  for (row = 0; row < Q->rows; row++) {
    for (col = 0; col < Q->cols; col++) {
      printf("%4d", Q->array.SZ[row][col]);
    }
    printf("\n");
  }
}

/* swapping, assumes they do not point to the same thing!!! */

int swap(int** Q, int row1, int col1, int row2, int col2)
{
  if ((row1 == row2) && (col1 == col2)) {
    printf("trying to swap an element with itself.\n");
    return -1;
  }
  Q[row1][col1] ^= Q[row2][col2];
  Q[row2][col2] ^= Q[row1][col1];
  Q[row1][col1] ^= Q[row2][col2];
  return 0;
}

/* resymmetrizing a matrix after working with upper triangular part */

int resymmetrize(int **Q)
{
  int row, col;
  
  for (row = 0; row < N-1; row++)
    for (col = row+1; col < N; col++)
      Q[col][row] = Q[row][col];

  return 0;
}

bravais_TYP* automorphism_group(matrix_TYP* Q)
{
  bravais_TYP* grp;
  int i, Qmax;
  matrix_TYP* SV;
  int options[6] = {0};
  
  Qmax = Q->array.SZ[0][0];
  for(i=1;i < Q->cols;i++) {
    if(Q->array.SZ[i][i] > Qmax)
      Qmax = Q->array.SZ[i][i];
  }
  SV = short_vectors(Q, Qmax, 0, 0, 0, &i);

  grp = autgrp(&Q, 1, SV, NULL, 0, options);

  free_mat(SV);
  
  return grp;
}

matrix_TYP* is_isometric(matrix_TYP* Q1, matrix_TYP* Q2)
{
  int i, Qmax;
  matrix_TYP *SV1, *SV2, *isom;
  int options[6] = {0};
  
  Qmax = Q1->array.SZ[0][0];
  for(i=1;i < Q1->cols;i++) {
    if(Q1->array.SZ[i][i] > Qmax)
      Qmax = Q1->array.SZ[i][i];
  }
  SV1 = short_vectors(Q1, Qmax, 0, 0, 0, &i);
  SV2 = short_vectors(Q2, Qmax, 0, 0, 0, &i);

  isom = isometry(&Q1, &Q2, 1, SV1, SV2, NULL, 0, options);
  
  free_mat(SV1);
  free_mat(SV2);

  return isom;
}

matrix_TYP* minkowski_reduce(matrix_TYP* Q)
{
  matrix_TYP *T1, *T2, *red1, *red2;

  /* printf("Trying to minkowski reduce: \n"); */
  /* print_mat(Q); */
  
  T1 = init_mat(N, N, "");
  red1 = pair_red(Q, T1);

  /* printf("After pair reduce, got: \n"); */
  /* print_mat(red1); */
  
  T2 = init_mat(N, N, "1");
  red2 = mink_red(red1, T2);

  /* printf("After Minkowski reduce, got: \n"); */
  /* print_mat(red1); */
  
  free_mat(T1);
  free_mat(T2);
  free_mat(red1);

  return red2;
}
