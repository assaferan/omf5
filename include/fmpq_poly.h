#ifndef __FMPQ_POLY_H__
#define __FMPQ_POLY_H__

#include <flint/fmpq_poly.h>
#include <flint/fmpq_mat.h>

#include "typedefs.h"

typedef struct
{
  fmpq_poly_struct *p;
  slong *exp;
  slong num;
}
fmpq_poly_factor_struct;

typedef fmpq_poly_factor_struct fmpq_poly_factor_t[1];

void fmpq_poly_factor(fmpq_poly_factor_t factors, const fmpq_poly_t f);
void fmpq_poly_factor_clear(fmpq_poly_factor_t factors);
bool fmpq_poly_is_irreducible(const fmpq_poly_t f);
void fmpq_poly_evaluate_fmpq_mat(fmpq_mat_t res, const fmpq_poly_t poly, const fmpq_mat_t a);

#endif // __FMPQ_POLY_H__
