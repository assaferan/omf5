#include "nf_mat.h"

void nf_mat_init(nf_mat_t mat, const nf_t* nf, slong nrows, slong ncols)
{
  slong row, col;

  mat->nf = nf;
  mat->rows = nrows;
  mat->cols = ncols;

  mat->array = (nf_elem_t**)malloc(nrows*sizeof(nf_elem_t*));
  for (row = 0; row < nrows; row++) {
    mat->array[row] = (nf_elem_t*)malloc(ncols*sizeof(nf_elem_t));
    for (col = 0; col < ncols; col++)
      nf_elem_init(mat->array[row][col], *(mat->nf));
  }

  return;
}

void nf_mat_set(nf_mat_t mat1, const nf_mat_t mat2)
{
  slong row, col;
  
  for (row = 0; row < mat1->rows; row++)
    for (col = 0; col < mat1->cols; col++)
      nf_elem_set(mat1->array[row][col], mat2->array[row][col], *(mat1->nf));
  
  return;
}

void nf_mat_init_set(nf_mat_t mat, const nf_mat_t src)
{

  nf_mat_init(mat, src->nf, src->rows, src->cols);
  nf_mat_set(mat, src);

  return;
}

nf_elem_t* nf_mat_entry(const nf_mat_t mat, slong i, slong j)
{
  return mat->array[i] + j;
}

void nf_mat_zero(nf_mat_t mat)
{
  slong row, col;
  
  for (row = 0; row < mat->rows; row++)
    for (col = 0; col < mat->cols; col++)
      nf_elem_zero(mat->array[row][col], *(mat->nf));

  return;
}

void nf_mat_one(nf_mat_t mat)
{
  slong row, max;

  nf_mat_zero(mat);

  max = (mat->rows < mat->cols) ? mat->rows : mat->cols;
  
  for (row = 0; row < max; row++)
    nf_elem_one(mat->array[row][row], *(mat->nf));

  return;
}

void fmpz_mat_get_nf_mat(nf_mat_t B, const fmpz_mat_t A)
{
  slong row, col;

  for (row = 0; row < fmpz_mat_nrows(A); row++)
    for (col = 0; col < fmpz_mat_ncols(A); col++)
      nf_elem_set_fmpz(B->array[row][col], fmpz_mat_entry(A,row,col), *(B->nf));

  return;
}

void fmpq_mat_get_nf_mat(nf_mat_t B, const fmpq_mat_t A)
{
  slong row, col;

  for (row = 0; row < fmpq_mat_nrows(A); row++)
    for (col = 0; col < fmpq_mat_ncols(A); col++)
      nf_elem_set_fmpq(B->array[row][col], fmpq_mat_entry(A,row,col), *(B->nf));

  return;
}

void nf_mat_clear(nf_mat_t mat)
{
  slong row, col;

  for (row = 0; row < mat->rows; row++) {
    for (col = 0; col < mat->cols; col++)
      nf_elem_clear(mat->array[row][col], *(mat->nf));
    free(mat->array[row]);
  }

  free(mat->array);
  
  return;
}

int nf_mat_print(const nf_mat_t mat)
{
  slong row, col;
  int num_printed = 0;

  num_printed = printf("%ld %ld  ", mat->rows, mat->cols);
  for (row = 0; row < mat->rows; row++) {
    for (col = 0; col < mat->cols; col++) {
      nf_elem_print_pretty(mat->array[row][col], *(mat->nf), "a");
      num_printed += printf(" ");
    }
  }
  
  return num_printed;
}

int nf_mat_print_pretty(const nf_mat_t mat)
{
  slong row, col;
  int num_printed = 0;

  num_printed = printf("[\n");
  for (row = 0; row < mat->rows; row++) {
    num_printed += printf("[");
    for (col = 0; col < mat->cols; col++) {
      nf_elem_print_pretty(mat->array[row][col], *(mat->nf), "a");
      num_printed += printf(" ");
    }
    num_printed += printf("]\n");
  }
  num_printed += printf("]\n");
  
  return num_printed;
}

void swap_rows(nf_mat_t mat, int row1, int row2)
{
  nf_elem_t tmp;
  slong col;

  nf_elem_init(tmp, *(mat->nf));
  nf_elem_zero(tmp, *(mat->nf));
  
  for (col = 0; col < mat->cols; col++) {
    nf_elem_set(tmp, mat->array[row1][col], *(mat->nf));
    nf_elem_set(mat->array[row1][col], mat->array[row2][col], *(mat->nf));
    nf_elem_set(mat->array[row2][col], tmp, *(mat->nf));
  }

  nf_elem_clear(tmp, *(mat->nf));
  return;
}

void multiply_row(nf_mat_t mat, int row, nf_elem_t scalar)
{
  slong col;
  
  for (col = 0; col < mat->cols; col++)
    nf_elem_mul(mat->array[row][col], mat->array[row][col], scalar, *(mat->nf));

  return;
}

slong nf_mat_rref_classical(nf_mat_t echelon, nf_mat_t trans)
{
  slong row, col;
  slong pivot_row, pivot_col, row_max;
  nf_elem_t max_val, scalar, factor, prod;

  trans->nf = echelon->nf;
  nf_elem_init(max_val, *(echelon->nf));
  nf_elem_init(factor, *(echelon->nf));
  nf_elem_init(scalar, *(echelon->nf));
  nf_elem_init(prod, *(echelon->nf));

  nf_mat_one(trans);

  pivot_row = pivot_col = 0;
  nf_elem_zero(max_val, *(echelon->nf));

#ifdef DEBUG_LEVEL_FULL
  printf("Computing row echelon form of\n");
  nf_mat_print_pretty(echelon);
#endif // DEBUG_LEVEL_FULL
  
  while ((pivot_row < echelon->rows) && (pivot_col < echelon->cols)) {
    for (row_max = pivot_row; row_max < echelon->rows;) {
      nf_elem_set(max_val,echelon->array[row_max][pivot_col], *(echelon->nf));
      if (nf_elem_is_zero(max_val, *(echelon->nf)))
	row_max++;
      else
	break;
    }
    if (nf_elem_is_zero(max_val, *(echelon->nf))) {
      pivot_col++;
    }
    else {
#ifdef DEBUG_LEVEL_FULL
      printf("swapping rows %ld and %ld.\n", pivot_row, row_max);
#endif // DEBUG_LEVEL_FULL
      swap_rows(echelon, pivot_row, row_max);
      swap_rows(trans, pivot_row, row_max);
      nf_elem_inv(scalar, max_val, *(echelon->nf));
#ifdef DEBUG_LEVEL_FULL
      printf("multiplying row %ld by ", pivot_row);
      nf_elem_print_pretty(scalar, *(echelon->nf), "x");
      printf(".\n");
#endif // DEBUG_LEVEL_FULL
      multiply_row(echelon, pivot_row, scalar);
      multiply_row(trans, pivot_row, scalar);
#ifdef DEBUG_LEVEL_FULL
      printf("zeroing the pivot column.\n");
#endif // DEBUG_LEVEL_FULL

      for ( row = 0; row < pivot_row; row++) {
	nf_elem_set(factor, echelon->array[row][pivot_col], *(echelon->nf));
	nf_elem_zero(echelon->array[row][pivot_col], *(echelon->nf));

	for ( col = pivot_col + 1; col < echelon->cols; col++) {
	  nf_elem_mul(prod, factor, echelon->array[pivot_row][col], *(echelon->nf));
	  nf_elem_sub(echelon->array[row][col], echelon->array[row][col], prod, *(echelon->nf));
	}
	for (col = 0; col < trans->cols; col++) {
	  nf_elem_mul(prod, factor, trans->array[pivot_row][col], *(trans->nf));
	  nf_elem_sub(trans->array[row][col], trans->array[row][col], prod, *(trans->nf));
	}
      }
      
      for ( row = pivot_row+1; row < echelon->rows; row++) {
	nf_elem_set(factor, echelon->array[row][pivot_col], *(echelon->nf));
	nf_elem_zero(echelon->array[row][pivot_col], *(echelon->nf));

	for (col = pivot_col + 1; col < echelon->cols; col++) {
	  nf_elem_mul(prod, factor, echelon->array[pivot_row][col], *(echelon->nf));
	  nf_elem_sub(echelon->array[row][col], echelon->array[row][col], prod, *(echelon->nf));
	}
	for (col = 0; col < trans->cols; col++) {
	  nf_elem_mul(prod, factor, trans->array[pivot_row][col], *(trans->nf));
	  nf_elem_sub(trans->array[row][col], trans->array[row][col], prod, *(trans->nf));
	}
      }
      
#ifdef DEBUG_LEVEL_FULL
      printf("echelon = \n");
      nf_mat_print_pretty(echelon);
      printf("trans = \n");
      nf_mat_print_pretty(trans);
#endif // DEBUG_LEVEL_FULL
      
      pivot_row++;
      pivot_col++;
    }
  }

  nf_elem_clear(max_val, *(echelon->nf));
  nf_elem_clear(scalar, *(echelon->nf));
  nf_elem_clear(factor, *(echelon->nf));
  nf_elem_clear(prod, *(echelon->nf));
  
  return pivot_row;
}

void nf_mat_left_kernel(nf_mat_t kernel, const nf_mat_t mat)
{
  nf_mat_t echelon, trans;
  slong rank, row, col;

  nf_mat_init_set(echelon, mat);
  nf_mat_init(trans, mat->nf, mat->rows, mat->cols);

  rank = nf_mat_rref_classical(echelon, trans);
  
  // getting the zero rows
  nf_mat_init(kernel, mat->nf, mat->rows-rank, mat->rows);

  // !! TODO - we could optimize this by using windows
  for (row = rank; row < mat->rows; row++)
    for (col = 0; col < mat->rows; col++)
      nf_elem_set(kernel->array[row-rank][col], trans->array[row][col], *(kernel->nf));
  
  nf_mat_clear(echelon);
  nf_mat_clear(trans);
   
  return;
}

void nf_mat_transpose(nf_mat_t tr_mat, const nf_mat_t mat)
{
  slong row, col;

  for (row = 0; row < mat->rows; row++)
    for (col = 0; col < mat->cols; col++)
      nf_elem_set(tr_mat->array[row][col], mat->array[col][row], *(tr_mat->nf));

  return;
}

void nf_mat_kernel(nf_mat_t ker, const nf_mat_t mat)
{
  nf_mat_t tr_mat;

  nf_mat_init(tr_mat, mat->nf, mat->cols, mat->rows);
  
  nf_mat_transpose(tr_mat, mat);
  nf_mat_left_kernel(ker, tr_mat);

  nf_mat_clear(tr_mat);

  return;
}

/* Addition, scalar multiplication */

// Sets mat to the difference of mat1 and mat2, assuming that all three matrices have the same dimensions.
void nf_mat_sub(nf_mat_t mat, const nf_mat_t mat1, const nf_mat_t mat2)
{
  slong row, col;

  for (row = 0; row < mat->rows; row++)
    for (col = 0; col < mat->cols; col++)
      nf_elem_sub(mat->array[row][col],mat1->array[row][col],mat2->array[row][col], *(mat->nf));

  return;
}

// Sets rop to op multiplied by the number field element x, assuming that the two matrices have the same dimensions.
void nf_mat_scalar_mul_nf(nf_mat_t rop, const nf_mat_t op, const nf_elem_t x)
{
  slong row, col;

  for (row = 0; row < op->rows; row++)
    for (col = 0; col < op->cols; col++)
      nf_elem_mul(rop->array[row][col], op->array[row][col], x, *(rop->nf));

  return;
}
