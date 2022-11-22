#include <assert.h>
#include <inttypes.h>

#include <carat/matrix.h>
#include <carat/symm.h>

#include "arith.h"
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

void square_matrix_init_set_symm_A(square_matrix_t Q, const Z64* coeff_vec)
{
  int row, col, idx;
  
  row = 0;
  col = 0;
  for (idx = 0; idx < QF_RANK*(QF_RANK+1)/2; idx++)
    {
      Q[row][col] = coeff_vec[idx];
      row++;
      if (row > col) {
	row = 0;
	col++;
      }
    }
  for (row = 0; row < QF_RANK-1; row++)
    for (col = row+1; col < QF_RANK; col++)
      Q[col][row] = Q[row][col];

  return;
}

void square_matrix_init_set_symm_GG(square_matrix_t Q, const Z64* coeff_vec)
{
  int row, col, idx;

  row = 0;
  col = 0;
  for (idx = 0; idx < QF_RANK*(QF_RANK+1)/2; idx++)
    {
      Q[row][col] = coeff_vec[idx];
      col++;
      if (col >= QF_RANK) {
	row++;
	col = row;
      }
    }
  for (row = 0; row < QF_RANK; row++)
    for (col = 0; col < row; col++)
      Q[row][col] = Q[col][row];

  for (row = 0; row < QF_RANK; row++)
    Q[row][row] *= 2;

#ifdef DEBUG
  square_matrix_print(Q);
#endif // DEBUG
  
  return;
}

void square_matrix_init_set_symm(square_matrix_t mat, const Z64* coeff_vec, const char* alg)
{
  if (strcmp(alg,"A") == 0)
    return square_matrix_init_set_symm_A(mat, coeff_vec);
  if (strcmp(alg,"GG") == 0)
    return square_matrix_init_set_symm_GG(mat, coeff_vec);

  assert(false);
  return;
}

void square_matrix_set(square_matrix_t dest, const square_matrix_t src)
{
  int i,j;

  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      dest[i][j] = src[i][j];

  return;
}

void vector_set(vector_t dest, const vector_t src)
{
  int i;
  for (i = 0; i < QF_RANK; i++)
    dest[i] = src[i];

  return;
}

void square_matrix_set_fmpz_mat(square_matrix_t dest, const fmpz_mat_t src)
{
  int i,j;

  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      dest[i][j] = fmpz_get_si(fmpz_mat_entry(src,i,j));

  return;
}

int square_matrix_set_matrix_TYP(square_matrix_t dest, matrix_TYP* src)
{
  int i,j;
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      dest[i][j] = src->array.SZ[i][j];

  return src->kgv;
}

void fmpz_mat_init_set_square_matrix(fmpz_mat_t dest, const square_matrix_t src)
{
  int i,j;
  fmpz_mat_init(dest, QF_RANK, QF_RANK);
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      fmpz_init_set_si(fmpz_mat_entry(dest,i,j), src[i][j]);

  return;
}

matrix_TYP* matrix_TYP_init_set_square_matrix(const square_matrix_t mat)
{
  int i,j;
  matrix_TYP* nmat;

  nmat = init_mat(QF_RANK, QF_RANK, "");
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      nmat->array.SZ[i][j] = mat[i][j];

  return nmat;
}

int vector_cmp(const vector_t vL, const vector_t vR)
{
  int i;
  
  for (i = 0; i < QF_RANK; i++) {
    if (vL[i] > vR[i])
      return 1;
    if (vL[i] < vR[i])
      return -1;
  }

  return 0;
}

bool square_matrix_is_equal(const square_matrix_t mat1, const square_matrix_t mat2)
{
  int i,j;
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      if (mat1[i][j] != mat2[i][j])
	return false;

  return true;
}

void square_matrix_zero(square_matrix_t mat)
{
  int i,j;

  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      mat[i][j] = 0;

  return;
}

void vector_zero(vector_t vec)
{
  int i;

  for (i = 0; i < QF_RANK; i++)
    vec[i] = 0;

  return;
}

void square_matrix_one(square_matrix_t mat)
{
  int i;

  square_matrix_zero(mat);
  for (i = 0; i < QF_RANK; i++)
    mat[i][i] = 1;

  return;
}

bool square_matrix_is_one(const square_matrix_t mat)
{
  int i,j;

  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      if (mat[i][j] != ((i == j) ? 1 : 0))
	return false;

  return true;
}

bool square_matrix_is_zero(const square_matrix_t mat)
{
  int i,j;

  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      if (mat[i][j] != 0)
	return false;

  return true;
}

bool vector_is_zero(const vector_t vec)
{
  int i;

  for (i = 0; i < QF_RANK; i++)
    if (vec[i] != 0)
      return false;
      
  return true;
}

bool square_matrix_is_positive_definite(const square_matrix_t mat)
{
  matrix_TYP* Q;
  bool is_pos_def;

  Q = matrix_TYP_init_set_square_matrix(mat);
  is_pos_def = definite_test(Q);

  free_mat(Q);
  
  return is_pos_def;
}

bool square_matrix_is_bad_prime(const square_matrix_t mat, Z64 p)
{
  fmpz_mat_t Q;
  fmpz_t det;
  bool is_bad;
  Z64 disc;

  fmpz_init(det);
  fmpz_mat_init_set_square_matrix(Q, mat);
  fmpz_mat_det(det, Q);
  disc = fmpz_get_si(det) / 2;
  
  is_bad = (disc % p == 0);

  fmpz_clear(det);
  fmpz_mat_clear(Q);
  return is_bad;
}

// !! TODO - we could do this without copying by having a view
void square_matrix_transpose(square_matrix_t tr, const square_matrix_t mat)
{
  int i,j;

  assert(tr != mat);
  
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      tr[j][i] = mat[i][j];

  return;
}

void square_matrix_add(square_matrix_t sum, const square_matrix_t matL, const square_matrix_t matR)
{
  int i,j;

  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      sum[i][j] = matL[i][j]+matR[i][j];

  return;
}

// !! TODO - this could be made faster
void square_matrix_mul(square_matrix_t prod, const square_matrix_t matL, const square_matrix_t matR)
{
  int i,j,k;

  assert((prod != matR) && (prod != matL));
  
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++) {
      prod[i][j] = 0;
      for (k = 0; k < QF_RANK; k++)
	prod[i][j] += matL[i][k]*matR[k][j];
    }

  return;
}

void square_matrix_muleq_right(square_matrix_t matL, const square_matrix_t matR)
{
  square_matrix_t prod;

  square_matrix_mul(prod, matL, matR);
  square_matrix_set(matL, prod);

  return;
}
void square_matrix_muleq_left(square_matrix_t matR, const square_matrix_t matL)
{
  square_matrix_t prod;

  square_matrix_mul(prod, matL, matR);
  square_matrix_set(matR, prod);

  return;
}

void square_matrix_mul_vec_left(vector_t prod, const vector_t vec, const square_matrix_t mat)
{
  int i, j;

  for (j = 0; j < QF_RANK; j++)
    prod[j] = 0;
  
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      prod[j] += vec[i]*mat[i][j];
  
  return;
}

int square_matrix_inv(square_matrix_t inv, const square_matrix_t mat, int denom)
{
  // !! TODO - at the moment we just use matrix_TYP*, see if we can do better
  matrix_TYP *s, *s_inv;
  int inv_denom, i, j;

  s = init_mat(QF_RANK,QF_RANK,"");
  for (i = 0 ; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      s->array.SZ[i][j] = mat[i][j];
  s->kgv = denom;

  s_inv = mat_inv(s);
  for (i = 0 ; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      inv[i][j] = s_inv->array.SZ[i][j];

  inv_denom = s_inv->kgv;
  
  free_mat(s);
  free_mat(s_inv);

  return inv_denom;
}

void square_matrix_div_scalar(square_matrix_t quo, const square_matrix_t mat, int denom)
{
  int i,j;
  
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++) {
      assert(quo[i][j] % denom == 0);
      quo[i][j] = mat[i][j] / denom;
    }

  return;
}

void vector_lin_comb(vector_t res, const vector_t v, const vector_t w, Z64 a_v, Z64 a_w)
{
  int i;

  for (i = 0; i < QF_RANK; i++)
    res[i] = a_v * v[i] + a_w * w[i];

  return;
}

void vector_mod_p(vector_t v, Z64 p)
{
  int i;
  for (i = 0; i < QF_RANK; i++) {
    v[i] %= p;
    if (v[i] < 0)
      v[i] += p;
  }

  return;
}

// assumes all coordinates of v are between 0 and p-1
void normalize_mod_p(vector_t v, Z64 p)
{
  int pivot, i;
  Z64 x,y;

  for (pivot = 0; (pivot < QF_RANK) && (v[pivot] == 0); ) pivot++;

  assert(pivot != QF_RANK);

  gcdext(v[pivot], p, &x, &y);

  assert((x*v[pivot]-1)%p == 0);

  for (i = 0; i < QF_RANK; i++) {
    v[i] = (x * v[i]) % p;
    if (v[i] < 0)
      v[i] += p;
  }
  
  return;
}
  

Z64 scalar_product(const vector_t v1, const vector_t v2)
{
  int i;
  Z64 prod = 0;

  for (i = 0; i < QF_RANK; i++)
    prod += v1[i]*v2[i];

  return prod;
}

void square_matrix_resymmetrize(square_matrix_t Q)
{
  int row, col;
  
  for (row = 0; row < QF_RANK-1; row++)
    for (col = row+1; col < QF_RANK; col++)
      Q[col][row] = Q[row][col];

  return;
}

void square_matrix_swap_elts(square_matrix_t Q, int row1, int col1, int row2, int col2)
{
  assert((row1 != row2) || (col1 != col2));

  Q[row1][col1] ^= Q[row2][col2];
  Q[row2][col2] ^= Q[row1][col1];
  Q[row1][col1] ^= Q[row2][col2];
  
  return;
}

void square_matrix_print(const square_matrix_t mat)
{
  int i,j;

  printf("[");
  for (i = 0; i < QF_RANK; i++) {
    printf("[");
    for (j = 0; j < QF_RANK; j++) {
      printf("%" PRId64 " ", mat[i][j]);
      if (j != QF_RANK - 1)
	printf(",");
    }
    printf("]");
    if (i != QF_RANK-1)
      printf(",");
  }
  printf("]");

  return;
}

void square_matrix_print_pretty(const square_matrix_t mat)
{
  int i,j;

  for (i = 0; i < QF_RANK; i++) {
    for (j = 0; j < QF_RANK; j++)
      printf("%4" PRId64 " ", mat[i][j]);
    printf("\n");
  }

  return;
}


void vector_print(const vector_t vec)
{
  int i;
  for (i = 0; i < QF_RANK; i++)
    printf("%4" PRId64 " ", vec[i]);
  printf("\n");

  return;
}
