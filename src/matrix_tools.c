#include <assert.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>

#include <antic/nf.h>
#include <antic/nf_elem.h>

#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly_factor.h>

#include "carat/autgrp.h"
#include "carat/reduction.h"
#include "carat/symm.h"

#include "matrix_tools.h"
#include "nf_mat.h"
#include "typedefs.h"

/* Right now working with matrices as in CARAT package
 * using matrix_TYP.
 * May consider using instead fixed size matrices
 * (without heap allocation, always using a constant size) 
 * may also consider other implementations such as flint or eigen
*/

/* For now, we use N = 5 hard coded */

#define N 5

/* function to initialize a symmetric matrix */
matrix_TYP* init_sym_matrix_A(const int* coeff_vec)
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

matrix_TYP* init_sym_matrix_GG(const int* coeff_vec)
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
      col++;
      if (col >= N) {
	row++;
	col = row;
      }
    }
  for (row = 0; row < N; row++)
    for (col = 0; col < row; col++)
      Q[row][col] = Q[col][row];

  for (row = 0; row < N; row++)
    Q[row][row] *= 2;

#ifdef DEBUG
  print_mat(Q_mat);
#endif // DEBUG
  
  return Q_mat;
  
}

matrix_TYP* init_sym_matrix(const int* coeff_vec, const char* alg)
{
  if (strcmp(alg,"A") == 0)
    return init_sym_matrix_A(coeff_vec);
  if (strcmp(alg,"GG") == 0)
    return init_sym_matrix_GG(coeff_vec);

  assert(FALSE);
  return NULL;
}

/* function to print matrices nicely
 */

int num_digits(int a, int base)
{
  int num, b;

  num = 0;
  b = a;
  while (b) {
    b /= base;
    num++;
  }
  if (num == 0)
    num++;

  // the minus sign
  if (a < 0)
    num++;
  
  return num;
}

void print_mat(const matrix_TYP* Q)
{
  int row, col, width;
  int* widths;

  widths = (int *)malloc(Q->cols * sizeof(int));

  // checking width
  for (col = 0; col < Q->cols; col++) {
    widths[col] = 0;
    for (row = 0; row < Q->rows; row++) {
      width = num_digits(Q->array.SZ[row][col], 10);
      if (widths[col] < width)
	widths[col] = width;
    }
  }
  
  for (row = 0; row < Q->rows; row++) {
    for (col = 0; col < Q->cols; col++) {
      printf("%*d", widths[col] + 1, Q->array.SZ[row][col]);
    }
    printf("\n");
  }
  
  free(widths);
  return;
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

Z64* voronoi_bounds(int dim)
{
  Z64* bounds;
  int i;
  
  bounds = (Z64 *)malloc(dim*sizeof(Z64));
  for (i = 0; i < dim; i++)
    bounds[i] = 1;
  return bounds; 
}

matrix_TYP* transform(matrix_TYP* g, matrix_TYP* Q)
{
  return mat_mul(tr_pose(g), mat_mul(Q, g));
}

matrix_TYP* transform_eq(matrix_TYP* g, matrix_TYP* Q)
{
  int i,j;
  
  Q = mat_muleq(Q,g);
  for (i = 0; i < N-1; i++)
    for (j = i+1; j < N; j++) {
      swap(Q->array.SZ, i, j, j, i);
    }
  return mat_muleq(Q, g);
  //  return mat_mul(tr_pose(g), mat_muleq(Q, g));
}

// right now using integer arithmetic. Should probably run faster with floating point arithmetic
void closest_lattice_vector(matrix_TYP* Q, matrix_TYP* iso, int dim)
{
  matrix_TYP *H_int, *g, *min_g, *x_gram;

  int **q;
  Z64 *voronoi, *x, *x_min, *x_max, *x_num, *x_closest;
  int i,j, num_xs, x_idx, min_dist, n;

  Z64 *y_int, *v_int;
  Z64 tmp, det;

  n = Q->rows;

#ifdef DEBUG_LEVEL_FULL
  printf("finding closest_lattice_vector with gram:\n");
  print_mat(Q);
#endif // DEBUG_LEVEL_FULL
  
  // v_int = init_mat(1, dim-1, "");
  v_int = (Z64 *)malloc((dim-1)*sizeof(Z64));
  g = init_mat(n, n, "1");
  min_g = init_mat(n, n, "");
  g->flags.Symmetric = g->flags.Scalar = g->flags.Diagonal = 0;
  H_int = adjugate(Q, dim-1);
  
  q = Q->array.SZ;
  
  for (i = 0; i < dim-1; i++) {
    v_int[i] = q[i][dim-1];
  }

#ifdef DEBUG_LEVEL_FULL
  printf("H_int = \n");
  print_mat(H_int);

  printf("v_int = \n");
  //  print_mat(v_int);
  for (i = 0; i < dim-1; i++)
    printf("%lld ", v_int[i]);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  y_int = (Z64 *)malloc((dim-1)*sizeof(Z64));
  // y_int = mat_mul(v_int, tr_pose(H_int));
  
  for (i = 0; i < dim-1; i++)
    y_int[i] = 0;

  for (i = 0; i < dim-1; i++)
    for (j = 0; j < dim-1; j++)
      y_int[i] += H_int->array.SZ[i][j] * v_int[j];

#ifdef DEBUG_LEVEL_FULL
  printf("y_int = \n");
  //  print_mat(y_int);
  for (i = 0; i < dim-1; i++)
    printf("%lld ", y_int[i]);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  voronoi = voronoi_bounds(dim-1);
  x = (Z64 *)malloc((dim-1)*sizeof(Z64));
  x_min = (Z64 *)malloc((dim-1)*sizeof(Z64));
  x_max = (Z64 *)malloc((dim-1)*sizeof(Z64));
  x_num = (Z64 *)malloc((dim-1)*sizeof(Z64));
  x_closest = (Z64 *)malloc((dim-1)*sizeof(Z64));
  
  det = 0;
  for (i = 0; i < dim - 1; i++) {
    tmp = H_int->array.SZ[0][i];
    det += tmp*q[i][0];
  }
  det = llabs(det);

  for ( i = 0; i < dim-1; i++) {
    tmp =  y_int[i] - det*voronoi[i];
    x_min[i] = ((tmp >= 0) ? tmp+det-1 : tmp)/det;
  }
  for ( i = 0; i < dim-1; i++) {
    tmp =  y_int[i] + det*voronoi[i];
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
      min_g = copy_mat(g);
      for (j = 0; j < dim-1; j++)
	x_closest[j] = x[j];
    }
  }

#ifdef DEBUG_LEVEL_FULL
  printf("x_closest = \n");
  for (i = 0; i < dim-1; i++)
    printf("%4lld ", x_closest[i]);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  iso = mat_muleq(iso, min_g);
  Q = transform_eq(min_g, Q);
#ifdef DEBUG_LEVEL_FULL
  printf("returning isometry: \n");
  print_mat(iso);
  printf("transformed gram to: \n");
  print_mat(Q);
#endif // DEBUG_LEVEL_FULL

  //  free_mat(v_int);
  free(v_int);
  free_mat(g);
  free_mat(min_g);
  free_mat(H_int);
  free(y_int);
  free(x);
  free(x_min);
  free(x_max);
  free(x_num);
  free(x_closest);
  return;
}

void _swap(int* pos, int i, int j)
{
  pos[i] ^= pos[j];
  pos[j] ^= pos[i];
  pos[i] ^= pos[j];

  return;
}

void sort(int* values, int* pos, int n)
{
  int cmp, i, j, k;
  int sorting_network_5_idx1[9] = {1,3,1,2,1,3,2,4,3};
  int sorting_network_5_idx2[9] = {2,4,3,5,2,4,3,5,4};
  
  switch(n) {
  case 2:
    if (values[0] > values[1]) {
      _swap(pos, 0, 1);
    }
    break;
    
  case 3:
    cmp = values[0] > values[1];
    cmp |= (values[0] > values[2]) << 1;
    cmp |= (values[1] > values[2]) << 2;
    switch(cmp) {
    case 0:
      // 0 <= 1 <= 2
      break;
    case 1:
      // 1 < 0 <=2
      _swap(pos, 0, 1);
      break;
    case 2:
      // 0 <= 1 <= 2, 0 > 2 can't really happen
      printf("Error! Shouldn't happen!\n");
      break;
    case 3:
      // 1 <= 2 < 0
      _swap(pos, 0, 1);
      _swap(pos, 1, 2);
      break;
    case 4:
      // 0 <= 2 < 1
      _swap(pos, 1, 2);
      break;
    case 5:
      // 0 > 1 > 2, 0 <=2
      printf("Error! Shouldn't happen!\n");
      break;
    case 6:
      // 2 < 0 <= 1
      _swap(pos, 0, 1);
      _swap(pos, 0, 2);
      break;
    case 7:
      // 0 > 1 > 2
      _swap(pos, 0, 2);
      break;
    }
    break;
    
  case 4:
    // This might be improved - currently 5 comparisons and 10 swaps.
    // It is the minimal number of comparisons, but we might be able to save on swaps.
    if (values[0] > values[1]) {
      _swap(pos, 0, 1);
      _swap(values, 0, 1);
    }
    if (values[2] > values[3]) {
      _swap(pos, 2, 3);
      _swap(values, 2, 3);
    }
    if (values[0] > values[2]) {
      _swap(pos, 0, 2);
      _swap(values, 0, 2);
    }
    if (values[1] > values[3]) {
      _swap(pos, 1, 3);
      _swap(values, 1, 3);
    }
    if (values[1] > values[2]) {
      _swap(pos, 1, 2);
      _swap(values, 1, 2);
    }
    break;

  case 5:
    // This might be improved - currently 9 comparisons and 18 swaps.
    // It is the minimal number of comparisons, but we might be able to save on swaps.
    for (i = 0; i < 9; i++) {
      j = sorting_network_5_idx1[i]-1;
      k = sorting_network_5_idx2[i]-1;
      if (values[j] > values[k]) {
	_swap(pos, j, k);
	_swap(values, j, k);
      }
    }
    break;
  }
  
  return;
}

void updatePerm(matrix_TYP* isom, int* perm)
{
  matrix_TYP* temp;
  int i,j,n;

  n = isom->cols;
  temp = init_mat(n,n,"");
  for ( i = 0; i < n; i++)
    for ( j = 0; j < n; j++)
      temp->array.SZ[i][j] = isom->array.SZ[i][j];
  for ( i = 0; i < n; i++)
    for ( j = 0; j < n; j++)
      isom->array.SZ[perm[i]][j] = temp->array.SZ[i][j];

  free_mat(temp);
  
  return;
}


// to avoid recursive template instantiation,
// we supply a parameter defining the level of recursion
// and use only this part of the matrices
// All containers will have size n, but we will only use dim entries
void greedy(matrix_TYP* gram, matrix_TYP* s, int n, int dim)
{
  matrix_TYP* tmp, *iso;
  int* perm_norm;
  int* perm;
  int i;
  
#ifdef DEBUG_LEVEL_FULL
  matrix_TYP *s0, *q0;
#endif // DEBUG_LEVEL_FULL

#ifdef DEBUG_LEVEL_FULL
  s0 = copy_mat(s);
  q0 = copy_mat(gram);
#endif // DEBUG_LEVEL_FULL
  
  if (dim == 1) return;

  perm_norm = (int*)malloc(dim*sizeof(int));
  perm = (int*)malloc(n*sizeof(int));
  
  do {
    for (i = 0; i < dim; i++) {
      perm_norm[i] = gram->array.SZ[i][i];
      perm[i] = i;
    }
    sort(perm_norm, perm, dim);

    // this is to make sure we do not touch these rows
    for ( i = dim; i < n; i++)
      perm[i] = i;

     // temp isometry
    tmp = init_mat(n,n,"1");
    tmp->flags.Scalar = tmp->flags.Diagonal = tmp->flags.Symmetric = 0;
    
    updatePerm(tmp, perm);

    // update isometry s = s*tmp;
    s = mat_muleq(s, tmp);
    
    // update gram
    gram = transform_eq(tmp, gram);

    free_mat(tmp);
    
/* #ifdef DEBUG_LEVEL_FULL */
/*     assert((s0.inverse()*s).transform(q0) == gram); */
/* #endif // DEBUG_LEVEL_FULL */

    // !! - TODO - do we really need iso here
    // or could we simply pass s?
    iso = init_mat(n,n,"1");
    iso->flags.Scalar = iso->flags.Diagonal = iso->flags.Symmetric = 0;
    greedy(gram, iso, n, dim-1);

    //    s = s*iso;
    s = mat_muleq(s, iso);

    free_mat(iso);
    // !! TODO - one can use subgram to save computations
    // This transformation already happens inside greedy(dim-1)
    //     gram = iso.transform(gram);

/* #ifdef DEBUG_LEVEL_FULL */
/*     assert((s0.inverse()*s).transform(q0) == gram); */
/* #endif // DEBUG_LEVEL_FULL */
    
    closest_lattice_vector(gram, s, dim);

/* #ifdef DEBUG_LEVEL_FULL */
/*     assert((s0.inverse()*s).transform(q0) == gram); */
/* #endif // DEBUG_LEVEL_FULL */
    
  } while (gram->array.SZ[dim-1][dim-1] < gram->array.SZ[dim-2][dim-2]);

  free(perm);
  free(perm_norm);
  
  return;

}

// assumes res is already initialized and is in the correct length
void fmpz_vec_mul_mat(fmpz_t* res, const fmpz_t* v, const fmpz_mat_t mat)
{
  slong i,j;
  fmpz_t *vT;

  vT = (fmpz_t*)malloc(fmpz_mat_ncols(mat) * sizeof(fmpz_t));
  for (j = 0; j < fmpz_mat_ncols(mat); j++) {
    fmpz_init(vT[j]);
    fmpz_zero(vT[j]);
  }
  
  for (i = 0; i < fmpz_mat_nrows(mat); i++)
    for (j = 0; j < fmpz_mat_ncols(mat); j++)
      fmpz_addmul(vT[j], v[i], fmpz_mat_entry(mat, i, j));

  for (j = 0; j < fmpz_mat_ncols(mat); j++) {
    fmpz_set(res[j], vT[j]);
    fmpz_clear(vT[j]);
  }
  
  free(vT);
  return;
}

void nf_elem_vec_mul_fmpq_mat(nf_elem_t* res, const nf_elem_t* v, const fmpq_mat_t mat, const nf_t nf)
{
  slong i,j;
  nf_elem_t tmp;

  assert(res != v);

  nf_elem_init(tmp, nf);
 
  for (j = 0; j < fmpq_mat_ncols(mat); j++) {
    nf_elem_zero(res[j], nf);
  }
  
  for (i = 0; i < fmpq_mat_nrows(mat); i++)
    for (j = 0; j < fmpq_mat_ncols(mat); j++) {
      nf_elem_scalar_mul_fmpq(tmp, v[i], fmpq_mat_entry(mat, i, j), nf);
      nf_elem_add(res[j], res[j], tmp, nf);
    }

  nf_elem_clear(tmp, nf);
  return;
}

slong get_pivot(const nf_elem_t* evec, slong dim, const nf_t nf)
{
  slong pivot;

  pivot = 0;
  while (nf_elem_is_zero(evec[pivot], nf)) {
    pivot++;
    assert(pivot < dim);
  }
  return pivot;
}

void check_eigenvector(const nf_elem_t* evec, const fmpq_mat_t M, const nf_t nf)
{
  nf_elem_t* evec_M;
  slong pivot, dim, i;
  nf_elem_t res1, res2;

  nf_elem_init(res1, nf);
  nf_elem_init(res2, nf);

  dim = fmpq_mat_nrows(M);
  pivot = get_pivot(evec, dim, nf);

  evec_M = (nf_elem_t*)malloc(dim * sizeof(nf_elem_t));
  for (i = 0; i < dim; i++)
    nf_elem_init(evec_M[i], nf);
  
  nf_elem_vec_mul_fmpq_mat(evec_M, evec, M, nf);
  
  for (i = 0; i < dim; i++) {
    nf_elem_mul(res1, evec[pivot], evec_M[i], nf);
    nf_elem_mul(res2, evec_M[pivot], evec[i], nf);
    assert(nf_elem_equal(res1, res2, nf));
  }
  
  for (i = 0; i < dim; i++)
    nf_elem_clear(evec_M[i], nf);
  free(evec_M);
  nf_elem_clear(res1, nf);
  nf_elem_clear(res2, nf);
  
  return;
}

// evec and nf should already be allocated and initialized
void mat_irred_charpoly_eigenvector(nf_elem_t* evec, const fmpq_mat_t mat, const fmpq_poly_t f, const nf_t nf)
{
  slong i, j, n;
  fmpz_t denom, g;
  fmpq_t coeff, denom_scale, scale;
  nf_elem_t a, b, tmp;
  nf_elem_t *c, *u;
  fmpz_mat_t scaled_mat;
  fmpz_t *v, *vv;
  BOOL is_zero;
  
  n = fmpq_poly_degree(f);

  assert( (fmpq_mat_nrows(mat) == n) && (fmpq_mat_ncols(mat) == n) );
  
  if (n == 1) {
    nf_elem_one(evec[0], nf);
    return;
  }

  nf_elem_init(a, nf);
  nf_elem_init(b, nf);
  nf_elem_init(tmp, nf);
  fmpq_init(coeff);
  fmpq_init(denom_scale);
  fmpq_init(scale);
  fmpz_init(denom);
  fmpz_init(g);
  fmpz_mat_init(scaled_mat, n, n);

  nf_elem_gen(a, nf);
  nf_elem_inv(b, a, nf);

  c = (nf_elem_t*) malloc(n*sizeof(nf_elem_t));
  u = (nf_elem_t*) malloc(n*sizeof(nf_elem_t));
  
  for (i = 0; i < n; i++)
    nf_elem_init(c[i], nf);

  for (i = 0; i < n; i++)
    nf_elem_init(u[i], nf);

  v = (fmpz_t*) malloc(n*sizeof(fmpz_t));
  vv = (fmpz_t*) malloc(n*sizeof(fmpz_t));

  for (i = 0; i < n; i++)
    fmpz_init(v[i]);
  
  for (i = 0; i < n; i++)
    fmpz_init(vv[i]);

  // c[0] = -b*f(0)
  fmpq_poly_get_coeff_fmpq(coeff, f, 0);
  fmpq_neg(coeff, coeff);
  nf_elem_scalar_mul_fmpq(c[0], b, coeff, nf);

  for (i = 1; i < n; i++) {
    // c[i] = (c[i-1] - coeff(f,i))*b
    fmpq_poly_get_coeff_fmpq(coeff, f, i);
    nf_elem_sub_fmpq(c[i], c[i-1], coeff, nf);
    nf_elem_mul(c[i], c[i], b, nf);
  }

  // scaling
  // denom := LCM([Denominator(x): x in Eltseq(T)]);
  // T := Matrix(S, denom*T);
  fmpq_mat_get_fmpz_mat_matwise(scaled_mat, denom, mat);
  // denom_scale = 1/denom
  fmpq_set_fmpz(denom_scale, denom);
  fmpq_inv(denom_scale, denom_scale);
  
  for (i = 0; i < n; i++)
    fmpz_zero(v[i]);

  srand(time(NULL));
  do {
    i = rand()/((RAND_MAX + 1u)/n);
    fmpz_add_ui(v[i],v[i],1);
    for (j = 0; j < n; j++)
      nf_elem_scalar_mul_fmpz(evec[j], c[0], v[j], nf);
    for (j = 0; j < n; j++)
      fmpz_set(vv[j], v[j]);
    fmpq_set(scale, denom_scale);
    for (i = 1; i < n; i++) {
      fmpz_vec_mul_mat(vv, vv, scaled_mat);
      for (j = 0; j < n; j++)
	nf_elem_set_fmpz(u[j], vv[j], nf);
      if (!(fmpq_is_one(denom_scale))) {
	fmpz_zero(g);
	for (j = 0; j < n; j++)
	  fmpz_gcd(g,g,vv[j]);
	if (!(fmpz_is_one(g))) {
	  for (j = 0; j < n; j++)
	    fmpz_divexact(vv[j], vv[j], g);
	  fmpq_mul_fmpz(scale, scale, g);
	}
	for (j = 0; j < n; j++)
	  nf_elem_set_fmpz(u[j], vv[j], nf);
	for (j = 0; j < n; j++)
	  nf_elem_scalar_mul_fmpq(u[j],u[j],scale,nf);
	fmpq_mul(scale, scale, denom_scale);
      }
      else {
	for (j = 0; j < n; j++)
	  nf_elem_set_fmpz(u[j], vv[j], nf);
      }
      for (j = 0; j < n; j++) {
	nf_elem_mul(tmp, c[i], u[j], nf);
	nf_elem_add(evec[j], evec[j], tmp, nf);
      }
    }
    is_zero = TRUE;
    for (j = 0; j < n; j++) {
      if (!(nf_elem_is_zero(evec[j], nf))) {
	is_zero = FALSE;
	break;
      }
    }
  } while (is_zero);

#ifdef DEBUG
  check_eigenvector(evec, mat, nf);
#endif // DEBUG
  
  for (i = 0; i < n; i++)
    fmpz_clear(vv[i]);
  free(vv);
  for (i = 0; i < n; i++)
    fmpz_clear(v[i]);
  free(v);
  for (i = 0; i < n; i++)
    nf_elem_clear(u[i], nf);
  free(u);
  for (i = 0; i < n; i++)
    nf_elem_clear(c[i], nf);
  free(c);
  fmpz_mat_clear(scaled_mat);
  fmpz_clear(g);
  fmpz_clear(denom);
  fmpq_clear(scale);
  fmpq_clear(denom_scale);
  fmpq_clear(coeff);
  nf_elem_clear(tmp, nf);
  nf_elem_clear(b, nf);
  nf_elem_clear(a, nf);
  
  return;
}

void check_restrict(const fmpq_mat_t T_W, const fmpq_mat_t T, const fmpq_mat_t W)
{
  fmpq_mat_t WT, T_W_W;

  fmpq_mat_init(WT, fmpq_mat_nrows(W), fmpq_mat_ncols(T));
  fmpq_mat_init(T_W_W, fmpq_mat_nrows(T_W), fmpq_mat_ncols(W));
  
  fmpq_mat_mul(WT, W, T);
  fmpq_mat_mul(T_W_W, T_W, W);

  assert(fmpq_mat_equal(WT, T_W_W));

  fmpq_mat_clear(WT);
  fmpq_mat_clear(T_W_W);
  return;
}

// This is what we do in magma. Could be faster to work with the integral matrices whenever we can
// assumes res_T is initialized
void restrict_mat(fmpq_mat_t res_T, const fmpq_mat_t T, const fmpq_mat_t basis_W)
{
  fmpq_mat_t B_t;
  fmpq_mat_t BA, BA_t;
  fmpq_mat_t B_BA_t, B_B_t;
  
  if (fmpq_mat_nrows(basis_W) == 0) {
    return;
  }

  fmpq_mat_init(BA, fmpq_mat_nrows(basis_W), fmpq_mat_ncols(T));
  fmpq_mat_init(BA_t, fmpq_mat_ncols(T), fmpq_mat_nrows(basis_W));
  fmpq_mat_mul(BA, basis_W, T);
  // we transpose, since we solve for the other direction
  fmpq_mat_transpose(BA_t, BA);
  fmpq_mat_init(B_t, fmpq_mat_ncols(basis_W), fmpq_mat_nrows(basis_W));
  fmpq_mat_transpose(B_t, basis_W);
  // transofrming to square matrices so flint will be able to solve it
  fmpq_mat_init(B_B_t, fmpq_mat_ncols(B_t), fmpq_mat_ncols(B_t));
  fmpq_mat_mul(B_B_t, basis_W, B_t);
  fmpq_mat_init(B_BA_t, fmpq_mat_ncols(B_t), fmpq_mat_ncols(BA_t));
  fmpq_mat_mul(B_BA_t, basis_W, BA_t);
  
  // fmpq_mat_solve_fraction_free(res_T, B_t, BA_t);
  fmpq_mat_solve_fraction_free(res_T, B_B_t, B_BA_t);

  // assert (!fmpq_mat_is_zero(res_T));
  
  fmpq_mat_transpose(res_T, res_T);

#ifdef DEBUG
  check_restrict(res_T, T, basis_W);
#endif // DEBUG
  

  fmpq_mat_clear(B_B_t);
  fmpq_mat_clear(B_BA_t);
  fmpq_mat_clear(BA_t);
  fmpq_mat_clear(BA);
  fmpq_mat_clear(B_t);

  return;
}

// initializes nf and evec (assumes evec is allocated)
void get_eigenvector(nf_elem_t* evec, nf_t nf, const fmpq_mat_t T, const fmpq_mat_t basis_W)
{
  slong i, dim_W, dim;
  fmpq_mat_t T_W;
  fmpq_poly_t f;
  nf_elem_t* evec_W;

  dim_W = fmpq_mat_nrows(basis_W);
  dim = fmpq_mat_ncols(basis_W);

  assert ((dim == fmpq_mat_nrows(T)) && (dim == fmpq_mat_ncols(T)) );
  
  evec_W = (nf_elem_t*) malloc(dim_W * sizeof(nf_elem_t));

  fmpq_mat_init(T_W, dim_W, dim_W);
  restrict_mat(T_W, T, basis_W);
  fmpq_poly_init(f);
  fmpq_mat_charpoly(f, T_W);
  nf_init(nf, f);
  for (i = 0; i < dim; i++)
    nf_elem_init(evec[i], nf);
  for (i = 0; i < dim_W; i++)
    nf_elem_init(evec_W[i], nf);
  
  mat_irred_charpoly_eigenvector(evec_W, T_W, f, nf);
  nf_elem_vec_mul_fmpq_mat(evec, evec_W, basis_W, nf);

#ifdef DEBUG
  check_eigenvector(evec, T, nf);
#endif // DEBUG
  
  for (i = 0; i < dim_W; i++)
    nf_elem_clear(evec_W[i], nf);
  free(evec_W);
  fmpq_mat_clear(T_W);
  
  return;
}

void fmpq_poly_factor(fmpq_poly_factor_t factors, const fmpq_poly_t f)
{
  slong i;
  fmpz_poly_t f_int;
  fmpz_poly_factor_t f_int_fac;

  fmpz_poly_init(f_int);
  fmpz_poly_factor_init(f_int_fac);
 
  fmpq_poly_get_numerator(f_int, f);
  fmpz_poly_factor_zassenhaus(f_int_fac, f_int);
  factors->num = f_int_fac->num;
  factors->p = (fmpq_poly_struct*)malloc(factors->num * sizeof(fmpq_poly_struct));
  factors->exp = (slong*)malloc(factors->num * sizeof(slong));
  for (i = 0; i < factors->num; i++) {
    factors->exp[i] = f_int_fac->exp[i];
    fmpq_poly_init(&(factors->p[i]));
    fmpq_poly_set_fmpz_poly(&(factors->p[i]), &(f_int_fac->p[i]));
  }
  
  fmpz_poly_factor_clear(f_int_fac);
  fmpz_poly_clear(f_int);
  return;
}

void fmpq_poly_factor_free(fmpq_poly_factor_t factors)
{
  free(factors->p);
  free(factors->exp);
}

void fmpq_poly_evaluate_fmpq_mat(fmpq_mat_t res, const fmpq_poly_t poly, const fmpq_mat_t a)
{
  slong i, n, r;
  fmpq_t coeff;
  fmpq_mat_t scalar;
  fmpz_mat_t num;
  fmpz_t denom, mat_gcd;

  fmpz_init(denom);
  fmpz_init(mat_gcd);
  fmpq_init(coeff);
  r = fmpq_mat_nrows(a);
  assert (r == fmpq_mat_ncols(a));
  fmpq_mat_init(scalar, r, r);
  fmpz_mat_init(num, r, r);
  
  fmpq_mat_zero(res);
  n = fmpq_poly_degree(poly);
  for (i = n; i >= 0; i--) {
    fmpq_mat_mul(res, res, a);
    fmpq_poly_get_coeff_fmpq(coeff, poly, i);
    fmpq_mat_one(scalar);
    fmpq_mat_scalar_mul_fmpq(scalar, scalar, coeff);
    fmpq_mat_add(res, res, scalar);
  }

  if (!fmpq_mat_is_zero(res)) {
    // reducing size of coefficients
    fmpq_mat_get_fmpz_mat_matwise(num, denom, res);
    fmpz_mat_content(mat_gcd, num);
    fmpq_mat_scalar_div_fmpz(res, res, mat_gcd);
  }

  fmpz_mat_clear(num);
  fmpq_mat_clear(scalar);
  fmpz_clear(denom);
  fmpz_clear(mat_gcd);
  fmpq_clear(coeff);
      
  return;
}

/* // initializes the kernel */
/* void fmpq_mat_left_kernel(fmpq_mat_t kernel, const fmpq_mat_t mat) */
/* { */
/*   fmpq_mat_t echelon, trans; */
/*   slong rank, row, col; */

/*   fmpq_mat_init_set(echelon, mat); */
/*   fmpq_mat_init(trans, fmpq_mat_nrows(mat), fmpq_mat_ncols(mat)); */

/*   rank = fmpq_mat_rref_classical(echelon, echelon); */
  
/*   // getting the zero rows */
/*   fmpq_mat_init(kernel, fmpq_mat_nrows(mat)-rank, fmpq_mat_nrows(mat)); */

/*   // !! TODO - we could optimize this by using windows */
/*   for (row = rank; row < fmpq_mat_nrows(mat); row++) */
/*     for (col = 0; col < fmpq_mat_nrows(mat); col++) */
/*       fmpq_set(fmpq_mat_entry(kernel,row-rank,col), fmpq_mat_entry(trans,row,col)); */
  
/*   fmpq_mat_clear(echelon); */
/*   fmpq_mat_clear(trans); */
   
/*   return; */
/* } */

// void fmpq_mat_left_kernel(fmpq_mat_t ker, const fmpq_mat_t mat)
void fmpq_mat_kernel(fmpq_mat_t ker, const fmpq_mat_t mat)
{
  fmpz_mat_t num, int_ker;
  fmpz_t den, mat_gcd;
  slong rank, row, col;

  fmpz_init(den);
  fmpz_init(mat_gcd);
  fmpz_mat_init(num, fmpq_mat_nrows(mat), fmpq_mat_ncols(mat));
  fmpz_mat_init(int_ker, fmpq_mat_ncols(mat), fmpq_mat_ncols(mat));

  fmpq_mat_get_fmpz_mat_matwise(num, den, mat);
  rank = fmpz_mat_nullspace(int_ker, num);

  fmpz_mat_content(mat_gcd, int_ker);
  fmpz_mat_scalar_divexact_fmpz(int_ker, int_ker, mat_gcd);
  
  //fmpq_mat_init(ker, fmpq_mat_ncols(mat), rank);
  fmpq_mat_init(ker, rank, fmpq_mat_ncols(mat)); 

  // fmpq_mat_set_fmpz_mat(ker, int_ker);
  for (row = 0; row < fmpq_mat_ncols(mat); row++)
    for (col = 0; col < rank; col++)
      fmpq_set_fmpz(fmpq_mat_entry(ker,col,row), fmpz_mat_entry(int_ker,row,col));
  
  fmpz_mat_clear(int_ker);
  fmpz_mat_clear(num);
  fmpz_clear(mat_gcd);
  fmpz_clear(den);
}

// not sure now if we are doing left kernel or right kernel. 
// void fmpq_mat_kernel(fmpq_mat_t ker, const fmpq_mat_t mat)
void fmpq_mat_left_kernel(fmpq_mat_t ker, const fmpq_mat_t mat)
{
  fmpq_mat_t tr_mat;

  fmpq_mat_init(tr_mat, fmpq_mat_ncols(mat), fmpq_mat_nrows(mat));
  
  fmpq_mat_transpose(tr_mat, mat);
  //  fmpq_mat_left_kernel(ker, tr_mat);
  fmpq_mat_kernel(ker, tr_mat);

  fmpq_mat_clear(tr_mat);
  
  return;
}

void print_content_and_coeff_size(const fmpq_mat_t A, const char* name)
{
  fmpz_mat_t num;
  fmpz_t mat_gcd, denom;
  mp_bitcnt_t coeff_size;
  slong i,j;

  fmpz_init(mat_gcd);
  fmpz_init(denom);
  
  fmpz_mat_init(num, fmpq_mat_nrows(A), fmpq_mat_ncols(A));
  fmpq_mat_get_fmpz_mat_matwise(num, denom, A);
  fmpz_mat_content(mat_gcd, num);
  coeff_size = 0;
  for (i = 0; i < fmpz_mat_nrows(num); i++)
    for (j = 0; j < fmpz_mat_ncols(num); j++)
      if (fmpz_bits(fmpz_mat_entry(num, i,j)) > coeff_size)
	coeff_size = fmpz_bits(fmpz_mat_entry(num, i,j));
  printf("content of %s is %ld, size of coefficients is %lu\n", name, fmpz_get_si(mat_gcd), coeff_size);
  
  fmpz_clear(mat_gcd);
  fmpz_clear(denom);
  fmpz_mat_clear(num);
  return;
}

// initializes ker
void kernel_on(fmpq_mat_t ker, const fmpq_mat_t A, const fmpq_mat_t B)
{
  fmpq_mat_t ker_A;
#ifdef DEBUG
  print_content_and_coeff_size(A, "A");
  print_content_and_coeff_size(B, "B");
#endif // DEBUG
  
  // fmpq_mat_kernel(ker, A);
  fmpq_mat_left_kernel(ker_A, A);

#ifdef DEBUG
  print_content_and_coeff_size(ker_A, "ker(A)");
#endif // DEBUG

  fmpq_mat_init(ker, fmpq_mat_nrows(ker_A), fmpq_mat_ncols(B));
  fmpq_mat_mul(ker, ker_A, B);

#ifdef DEBUG
  print_content_and_coeff_size(ker, "ker(A)*B");
#endif // DEBUG

  fmpq_mat_clear(ker_A);
  return;
}

void fmpq_mat_init_set_matrix_TYP(fmpq_mat_t M, const matrix_TYP* mat)
{
  slong row, col;
  
  fmpq_mat_init(M, mat->rows, mat->cols);

  for (row = 0; row < mat->rows; row++)
    for (col = 0; col  < mat->cols; col++)
      fmpq_set_si(fmpq_mat_entry(M, row, col), mat->array.SZ[row][col], 1);

  return;
}

eigenvalues* get_eigenvalues(matrix_TYP* mat)
{
  fmpz_mat_t M;
  fmpz_poly_t cp;
  fmpz_poly_factor_t cp_fac;
  fmpq_poly_t nf_poly;
  int row, col, i, j;
  eigenvalues* evs;
  nf_mat_t M_K, ker, lambda;
  
  fmpz_mat_init(M, mat->rows, mat->cols);

  for (row = 0; row < mat->rows; row++)
    for (col = 0; col  < mat->cols; col++)
      *fmpz_mat_entry(M, row, col) = mat->array.SZ[row][col];

  fmpz_poly_init(cp);
  fmpz_mat_charpoly(cp, M);
#ifdef DEBUG_LEVEL_FULL
  fmpz_poly_print_pretty(cp, "x");
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  fmpz_poly_factor_init(cp_fac);
  fmpz_poly_factor_zassenhaus(cp_fac, cp);

#ifdef DEBUG_LEVEL_FULL
  for (i = 0; i < cp_fac->num; i++) {
    fmpz_poly_print_pretty(&(cp_fac->p[i]), "x");
    printf("\n");
  }
#endif // DEBUG_LEVEL_FULL

  evs = (eigenvalues*)malloc(sizeof(eigenvalues));
  
  evs->num = cp_fac->num;
  evs->dim = mat->rows;
  evs->nfs = (nf_t*)malloc(cp_fac->num * sizeof(nf_t));
  evs->eigenvals = (nf_elem_t*)malloc(cp_fac->num * sizeof(nf_elem_t));
  evs->eigenvecs = (nf_elem_t**)malloc(cp_fac->num * sizeof(nf_elem_t*));
  
  for (i = 0; i < cp_fac->num; i++) {
    evs->eigenvecs[i] = (nf_elem_t*)malloc(mat->rows * sizeof(nf_elem_t));
    fmpq_poly_init(nf_poly);
    fmpq_poly_set_fmpz_poly(nf_poly, &(cp_fac->p[i]));
    nf_init(evs->nfs[i], nf_poly);
    nf_elem_init(evs->eigenvals[i], evs->nfs[i]);
    nf_elem_gen(evs->eigenvals[i], evs->nfs[i]);
    nf_mat_init(M_K, &(evs->nfs[i]), mat->rows, mat->cols);
    nf_mat_init(lambda, &(evs->nfs[i]), mat->rows, mat->cols);
    fmpz_mat_get_nf_mat(M_K, M);
    nf_mat_one(lambda);
    nf_mat_scalar_mul_nf(lambda, lambda, evs->eigenvals[i]);
    nf_mat_sub(M_K, M_K, lambda);
    nf_mat_kernel(ker, M_K);
#ifdef DEBUG
    printf("ker = \n");
    nf_mat_print_pretty(ker);
#endif // DEBUG
    for (j = 0; j < mat->rows; j++) {
      nf_elem_init(evs->eigenvecs[i][j], evs->nfs[i]);
      nf_elem_set(evs->eigenvecs[i][j], *nf_mat_entry(ker, 0, j), evs->nfs[i]);
    }
    nf_mat_clear(ker);
    nf_mat_clear(M_K);
    nf_mat_clear(lambda);
  }

  fmpz_mat_clear(M);
  fmpz_poly_factor_clear(cp_fac);
  fmpz_poly_clear(cp);
  
  return evs;
}

void eigenvalues_init(eigenvalues** evs, slong num, slong dim)
{
  slong i;
  
  (*evs)->num = num;
  (*evs)->dim = dim;
  (*evs)->nfs = (nf_t*)malloc(num * sizeof(nf_t));
  (*evs)->eigenvals = (nf_elem_t*)malloc(num * sizeof(nf_elem_t));
  (*evs)->eigenvecs = (nf_elem_t**)malloc(num * sizeof(nf_elem_t*));
  for (i = 0; i < num; i++)
    (*evs)->eigenvecs[i] = (nf_elem_t*)malloc(dim * sizeof(nf_elem_t));
  
  return;
}

void free_eigenvalues(eigenvalues* evs)
{
  int i, j;

  for (i = 0; i < evs->num; i++) {
    for (j = 0; j < evs->dim; j++) {
      nf_elem_clear(evs->eigenvecs[i][j], evs->nfs[i]);
    }
    free(evs->eigenvecs[i]);
    nf_elem_clear(evs->eigenvals[i], evs->nfs[i]);
    nf_clear(evs->nfs[i]);
  }

  free(evs->eigenvecs);
  free(evs->eigenvals);
  free(evs->nfs);
  
  free(evs);
  
  return;
}
