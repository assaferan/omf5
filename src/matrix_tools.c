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
  
#ifdef DEBUG_LEVEL_FULL
  printf("After pair reduce, got: \n");
  print_mat(red1);
#endif // DEBUG_LEVEL_FULL
  
  T2 = init_mat(N, N, "1");
  red2 = mink_red(red1, T2);

#ifdef DEBuG_LEVEL_FULL
  printf("After Minkowski reduce, got: \n");
  print_mat(red1);
#endif // DEBUG_LEVEL_FULL
  
  free_mat(T1);
  free_mat(T2);
  free_mat(red1);

  return red2;
}

matrix_TYP* adjugate(matrix_TYP* Q, int dim)
{
  // We will use this only for dim <= 4
  // and we write it down explicitly for each case
  // !! TODO - This is not the most effective way to do this
  matrix_TYP* adj;
  int** a, **b;

  adj = init_mat(dim, dim, "");
  a = Q->array.SZ;
  b = adj->array.SZ;
  
  switch(dim) {
  case 1:
    b[0][0] = 1;
    break;
  case 2:
    b[0][0] = a[1][1];
    b[1][1] = a[0][0];
    b[0][1] = -a[0][1];
    b[1][0] = -a[1][0];
    break;
  case 3:
    b[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
    b[1][0] = - a[1][0]*a[2][2] + a[1][2]*a[2][0];
    b[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
    b[0][1] = - a[0][1]*a[2][2] + a[2][1]*a[0][2];
    b[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
    b[2][1] = - a[0][0]*a[2][1] + a[0][1]*a[2][0];
    b[0][2] = a[0][1]*a[1][2] - a[1][1]*a[0][2];
    b[1][2] = - a[0][0]*a[1][2] + a[1][0]*a[0][2];
    b[2][2] = a[1][1]*a[0][0] - a[1][0]*a[0][1];
    break;
  case 4:
    b[0][0] = a[1][1]*(a[2][2]*a[3][3]-a[2][3]*a[3][2]);
    b[0][0] -= a[1][2]*(a[2][1]*a[3][3]-a[2][3]*a[3][1]);
    b[0][0] += a[1][3]*(a[2][1]*a[3][2]-a[2][2]*a[3][1]);
    b[0][1] = -a[1][0]*(a[2][2]*a[3][3]-a[2][3]*a[3][2]);
    b[0][1] += a[1][2]*(a[2][0]*a[3][3]-a[2][3]*a[3][0]);
    b[0][1] -= a[1][3]*(a[2][0]*a[3][2]-a[2][2]*a[3][0]);
    b[0][2] = a[1][0]*(a[2][1]*a[3][3]-a[2][3]*a[3][1]);
    b[0][2] -= a[1][1]*(a[2][0]*a[3][3]-a[2][3]*a[3][0]);
    b[0][2] += a[1][3]*(a[2][0]*a[3][1]-a[2][1]*a[3][0]);
    b[0][3] = -a[1][0]*(a[2][1]*a[3][2]-a[2][2]*a[3][1]);
    b[0][3] += a[1][1]*(a[2][0]*a[3][2]-a[2][2]*a[3][0]);
    b[0][3] -= a[1][2]*(a[2][0]*a[3][1]-a[2][1]*a[3][0]);
    b[1][0] = -a[0][1]*(a[2][2]*a[3][3]-a[2][3]*a[3][2]);
    b[1][0] += a[2][1]*(a[0][2]*a[3][3]-a[3][2]*a[0][3]);
    b[1][0] -= a[3][1]*(a[0][2]*a[2][3]-a[2][2]*a[0][3]);
    b[1][1] = a[0][0]*(a[2][2]*a[3][3]-a[2][3]*a[3][2]);
    b[1][1] -= a[0][2]*(a[2][0]*a[3][3]-a[2][3]*a[3][0]);
    b[1][1] += a[0][3]*(a[2][0]*a[3][2]-a[2][2]*a[3][0]);
    b[1][2] = -a[0][0]*(a[2][1]*a[3][3]-a[2][3]*a[3][1]);
    b[1][2] += a[0][1]*(a[2][0]*a[3][3]-a[3][0]*a[2][3]);
    b[1][2] -= a[0][3]*(a[2][0]*a[3][1]-a[3][0]*a[2][1]);
    b[1][3] = a[0][0]*(a[2][1]*a[3][2]-a[3][1]*a[2][2]);
    b[1][3] -= a[0][1]*(a[2][0]*a[3][2]-a[3][0]*a[2][2]);
    b[1][3] += a[0][2]*(a[2][0]*a[3][1]-a[2][1]*a[3][0]);
    b[2][0] = a[0][1]*(a[1][2]*a[3][3]-a[3][2]*a[1][3]);
    b[2][0] -= a[1][1]*(a[0][2]*a[3][3]-a[3][2]*a[0][3]);
    b[2][0] += a[3][1]*(a[0][2]*a[1][3]-a[1][2]*a[0][3]);
    b[2][1] = -a[0][0]*(a[1][2]*a[3][3]-a[3][2]*a[1][3]);
    b[2][1] += a[1][0]*(a[0][2]*a[3][3]-a[0][3]*a[3][2]);
    b[2][1] -= a[3][0]*(a[0][2]*a[1][3]-a[0][3]*a[1][2]);
    b[2][2] = a[0][0]*(a[1][1]*a[3][3]-a[1][3]*a[3][1]);
    b[2][2] -= a[0][1]*(a[1][0]*a[3][3]-a[1][3]*a[3][0]);
    b[2][2] += a[0][3]*(a[1][0]*a[3][1]-a[1][1]*a[3][0]);
    b[2][3] = -a[0][0]*(a[1][1]*a[3][2]-a[1][2]*a[3][1]);
    b[2][3] += a[0][1]*(a[1][0]*a[3][2]-a[3][0]*a[1][2]);
    b[2][3] -= a[0][2]*(a[1][0]*a[3][1]-a[3][0]*a[1][1]);
    b[3][0] = -a[0][1]*(a[1][2]*a[2][3]-a[2][2]*a[1][3]);
    b[3][0] += a[1][1]*(a[0][2]*a[2][3]-a[2][2]*a[0][3]);
    b[3][0] -= a[2][1]*(a[0][2]*a[1][3]-a[1][2]*a[0][3]);
    b[3][1] = a[0][0]*(a[1][2]*a[2][3]-a[1][3]*a[2][2]);
    b[3][1] -= a[1][0]*(a[0][2]*a[2][3]-a[0][3]*a[2][2]);
    b[3][1] += a[2][0]*(a[0][2]*a[1][3]-a[1][2]*a[0][3]);
    b[3][2] = -a[0][0]*(a[1][1]*a[2][3]-a[2][1]*a[1][3]);
    b[3][2] += a[1][0]*(a[0][1]*a[2][3]-a[0][3]*a[2][1]);
    b[3][2] -= a[2][0]*(a[0][1]*a[1][3]-a[0][3]*a[1][1]);
    b[3][3] = a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]);
    b[3][3] -= a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]);
    b[3][3] += a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
    break;
  default:
    printf("Error! Trying to find adjugate of matrix of size %d!\n", dim);
  }

  return adj;
}

int* voronoi_bounds(int dim)
{
  int* bounds, i;
  
  bounds = (int *)malloc(dim*sizeof(int));
  for (i = 0; i < dim; i++)
    bounds[i] = 1;
  return bounds; 
}

matrix_TYP* transform(matrix_TYP* g, matrix_TYP* Q)
{
  return mat_mul(tr_pose(g), mat_mul(Q, g));
}

void closest_lattice_vector(matrix_TYP* Q, matrix_TYP* iso, int dim)
{
  matrix_TYP *H_int, *v_int, *y_int, *g, *min_g, *x_gram;

  int **q;
  int *voronoi, *x, *x_min, *x_max, *x_num, *x_closest;
  int i,j, det, tmp, num_xs, x_idx, min_dist;

  v_int = init_mat(1, dim-1, "");
  g = init_mat(dim-1, dim-1, "");
  H_int = adjugate(Q, dim-1);
  
  q = Q->array.SZ;
  
  for (i = 0; i < dim-1; i++) {
    v_int->array.SZ[0][i] = q[i][dim-1];
  }

#ifdef DEBUG_LEVEL_FULL
  printf("H_int = \n");
  print_mat(H_int);

  printf("v_int = \n");
  print_mat(v_int);
#endif // DEBUG_LEVEL_FULL
  
  y_int = mat_mul(v_int, tr_pose(H_int));

#ifdef DEBUG_LEVEL_FULL
  printf("y_int = \n");
  print_mat(y_int);
#endif // DEBUG_LEVEL_FULL
  
  voronoi = voronoi_bounds(dim-1);
  x = (int *)malloc((dim-1)*sizeof(int));
  x_min = (int *)malloc((dim-1)*sizeof(int));
  x_max = (int *)malloc((dim-1)*sizeof(int));
  x_num = (int *)malloc((dim-1)*sizeof(int));
  x_closest = (int *)malloc((dim-1)*sizeof(int));
  
  det = 0;
  for (i = 0; i < dim - 1; i++)
    det += H_int->array.SZ[0][i]*q[i][0];
  det = abs(det);

  for ( i = 0; i < dim-1; i++) {
    tmp =  y_int->array.SZ[0][i] - det*voronoi[i];
    x_min[i] = ((tmp >= 0) ? tmp+det-1 : tmp)/det;
  }
  for ( i = 0; i < dim-1; i++) {
    tmp =  y_int->array.SZ[0][i] + det*voronoi[i];
    x_max[i] = ((tmp >= 0) ? tmp : tmp-det+1)/det;
  }
  
  for (i = 0; i < dim-1; i++)
    x_num[i] = x_max[i] - x_min[i] + 1;
  
  num_xs = 1;
  for (i = 0; i < dim-1; i++)
    num_xs *= x_num[i];
  // This should be infinity
  min_dist = INT_MAX;
  for (x_idx = 0; x_idx < num_xs; x_idx++) {
    tmp = x_idx;
    for ( i = 0; i < dim-1; i++) {
      j = dim-2-i;
      x[j] = x_min[j] + (tmp % x_num[j]);
      tmp /= x_num[j];
    }
    for ( i = 0; i < dim-1; i++)
      g->array.SZ[i][dim-1] = -x[i];
    x_gram = transform(g,Q);
    if (x_gram->array.SZ[dim-1][dim-1] < min_dist) {
      min_dist = x_gram->array.SZ[dim-1][dim-1];
      min_g = g;
      x_closest = x;
    }
  }

#ifdef DEBUG_LEVEL_FULL
  printf("x_closest = \n");
  for (i = 0; i < dim-1; i++)
    printf("%4d ", x_closest[i]);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  iso = mat_mul(iso, min_g);
  Q = transform(min_g, Q);
#ifdef DEBUG_LEVEL_FULL
  printf("returning isometry: \n");
  print_mat(iso);
  printf("transformed gram to: \n");
  print_mat(Q);
#endif // DEBUG_LEVEL_FULL
  return;
}
