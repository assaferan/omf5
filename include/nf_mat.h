#ifndef __NF_MAT_H__
#define __NF_MAT_H__

/*
 *
 * matrices over number fields
 * 
 */

#include <antic/nf.h>
#include <antic/nf_elem.h>

typedef struct {
  const nf_t* nf; // a pointer to the number field (not created)
  nf_elem_t** array;
  slong rows;
  slong cols;
} nf_mat_struct;

typedef nf_mat_struct nf_mat_t[1];

/* Memory management */

// Initialises a matrix with the given number of rows and columns for use. //
void nf_mat_init(nf_mat_t mat, const nf_t* nf, slong rows, slong cols);

// Clears the given matrix //
void nf_mat_clear(nf_mat_t mat);

/* Basic assignment and manipulation */

// Sets mat1 to a copy of mat2. The dimensions and number fields of mat1 and mat2 must be the same.
void nf_mat_set(nf_mat_t mat1, const nf_mat_t mat2);

// Initialises the matrix mat to the same size as src and sets it to a copy of src.
void nf_mat_init_set(nf_mat_t mat, const nf_mat_t src);

// Returns a reference to the entry of mat at row i and column j.
// This reference can be passed as an input or output variable to any function in the nf_elem module for direct manipulation.
// Both i and j must not exceed the dimensions of the matrix.
nf_elem_t* nf_mat_entry(const nf_mat_t mat, slong i, slong j);

// Sets all entries of mat to 0.
void nf_mat_zero(nf_mat_t mat);

// Sets mat to the unit matrix, having ones on the main diagonal and zeroes elsewhere. If mat is nonsquare, it is set to the truncation of a unit matrix.
void nf_mat_one(nf_mat_t mat);

// Sets the matrix rop to the transpose of the matrix op, assuming that their dimensions are compatible.
void nf_mat_transpose(nf_mat_t rop, const nf_mat_t op);

/* Addition, scalar multiplication */

// Sets mat to the difference of mat1 and mat2, assuming that all three matrices have the same dimensions,
// and the same number field.
void nf_mat_sub(nf_mat_t mat, const nf_mat_t mat1, const nf_mat_t mat2);

// Sets rop to op multiplied by the number field element x, assuming that the two matrices have the same dimensions, and the same number field.
void nf_mat_scalar_mul_nf(nf_mat_t rop, const nf_mat_t op, const nf_elem_t x);

/* Input and Output */

// Prints the given matrix to the stream stdout. The format is the number of rows, a space, the number of columns, two spaces, then a space separated list of coefficients, one row after the other.
// In case of success, returns a positive value; otherwise, returns a non-positive value.
int nf_mat_print(const nf_mat_t mat);

// Prints the given matrix to stdout.The format is an opening square bracket then on each line a row of the matrix, followed by a closing square bracket. Each row is written as an opening square bracket followed by a space separated list of coefficients followed by a closing square bracket.
// In case of success, returns a positive value; otherwise, returns a non-positive value.
int nf_mat_print_pretty(const nf_mat_t mat);

/* Conversions */

// Sets the entries of B as nf_elems corresponding to the entries of A.
void fmpz_mat_get_nf_mat(nf_mat_t B, const fmpz_mat_t A);

/* Echelon form */

// Replaces E by its reduced row echelon form, and returns the rank. Also returns the transformation T such that T*E is the result. Performs Gauss-Jordan elimination directly over the rational numbers. This algorithm is usually inefficient.
slong nf_mat_rref_classical(nf_mat_t E, nf_mat_t T);

// Returns in ker the kernel of the matrix mat (also initializes kernel)
void nf_mat_kernel(nf_mat_t ker, const nf_mat_t mat);

#endif // __NF_MAT_H__
