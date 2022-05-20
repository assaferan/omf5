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
