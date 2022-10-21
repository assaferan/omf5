#include <assert.h>

#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpz_poly_factor.h>

#include "fmpq_poly.h"
#include "typedefs.h"

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

void fmpq_poly_factor_clear(fmpq_poly_factor_t factors)
{
  free(factors->p);
  free(factors->exp);
}

bool fmpq_poly_is_irreducible(const fmpq_poly_t f)
{
  fmpq_poly_factor_t fac;
  bool ret;

  fmpq_poly_factor(fac, f);
  ret = (fac->num == 1) && (fac->exp[0] == 1);
  
  fmpq_poly_factor_clear(fac);
  return ret;
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
