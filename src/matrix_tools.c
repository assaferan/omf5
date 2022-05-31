#include <limits.h>

#include <antic/nf.h>
#include <antic/nf_elem.h>

#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly_factor.h>

#include "carat/autgrp.h"
#include "carat/reduction.h"
#include "carat/symm.h"

#include "matrix_tools.h"
#include "nf_mat.h"

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

void print_mat(matrix_TYP* Q)
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

void closest_lattice_vector(matrix_TYP* Q, matrix_TYP* iso, int dim)
{
  matrix_TYP *H_int, *v_int, *y_int, *g, *min_g, *x_gram;

  int **q;
  int *voronoi, *x, *x_min, *x_max, *x_num, *x_closest;
  int i,j, det, tmp, num_xs, x_idx, min_dist, n;

  n = Q->rows;

#ifdef DEBUG_LEVEL_FULL
  printf("finding closest_lattice_vector with gram:\n");
  print_mat(Q);
#endif // DEBUG_LEVEL_FULL
  
  v_int = init_mat(1, dim-1, "");
  g = init_mat(n, n, "1");
  min_g = init_mat(n, n, "");
  g->flags.Symmetric = g->flags.Scalar = g->flags.Diagonal = 0;
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
      min_g = copy_mat(g);
      for (j = 0; j < dim-1; j++)
	x_closest[j] = x[j];
    }
  }

#ifdef DEBUG_LEVEL_FULL
  printf("x_closest = \n");
  for (i = 0; i < dim-1; i++)
    printf("%4d ", x_closest[i]);
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

  free_mat(v_int);
  free_mat(g);
  free_mat(min_g);
  free_mat(H_int);
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
