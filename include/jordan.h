#ifndef __JORDAN_H__
#define __JORDAN_H__

#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz_mat.h>

#include "square_matrix.h"

typedef struct {
  fmpq_mat_t* matrices;
  fmpq_mat_t* grams;
  ulong* exponents;
  ulong num_blocks;
} jordan_data;

typedef jordan_data jordan_data_t[1];

void jordan_data_init(jordan_data_t jordan, ulong n);

void jordan_data_clear(jordan_data_t jordan);

void inner_product(fmpq_t res, const fmpz_mat_t G,
		   const fmpq_mat_t S, ulong idx1, ulong idx2);

void jordan_decomposition(jordan_data_t jordan, const square_matrix_t q, const fmpz_t p);

#endif // __JORDAN_H__
