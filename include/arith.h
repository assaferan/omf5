#ifndef __ARITH_H__
#define __ARITH_H__

#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fq.h>

#include "carat/typedef.h"

#include "typedefs.h"

int fmpz_valuation_unit(fmpz_t a, const fmpz_t x, const fmpz_t p);

int fmpz_valuation(const fmpz_t x, const fmpz_t p);

int fmpq_valuation(const fmpq_t x, const fmpz_t p);

int hilbert_symbol(const fmpz_t x, const fmpz_t y, const fmpz_t p);

/* Recursive function for a temporary extended Euclidean algorithm. */
/* It uses pointers to return multiple values. */

int gcdext(int a, int b, int *x, int *y);

int rational_lt(rational a, rational b);

rational rational_sum(rational a, rational b);

void bernoulli_number(fmpq_t b, ulong n);

void binomial_coefficient(fmpz_t b, const fmpz_t n, const fmpz_t k);

slong kronecker_symbol(const fmpz_t a, const fmpz_t n);

void bernoulli_number_chi(fmpq_t b_chi, ulong n, const fmpz_t d);

bool fmpz_is_local_square(const fmpz_t a, const fmpz_t p);

bool fmpq_is_local_square(const fmpq_t a, const fmpz_t p);

bool fq_is_square(const fq_t a, const fq_ctx_t F);

void fq_sqrt(fq_t sqrt_a, const fq_t a, const fq_ctx_t F);

bool fmpq_is_integral(const fmpq_t r);

void fmpq_floor(fmpz_t res, const fmpq_t r);

int primes_up_to(int** ps, int bound);

#endif // __ARITH_H
