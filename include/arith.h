/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : arith.h
 *
 * Description: simple useful arithmetic functions
 *
 **********************************************************
 */

#ifndef __ARITH_H__
#define __ARITH_H__

// Required packages dependencies

#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fq.h>

#include <carat/typedef.h>

// Self dependencies

#include "typedefs.h"

/*************************************************************************
 *
 * Function: fmpz_valuation_unit
 *
 * Description: given an integer x and a prime p, return v = v_p(x) -
 *              the valuation of x at p, and the corresponding unit
 *              a = p^(-v) x. x should be nonzero.
 *
 * Arguments:
 *     + x (const fmpz_t) - the integer whose valuation is to be computed
 *     + p (const fmpz_t) - the prime
 *
 * Returns:
 *     + (int) - the valuation v = v_p(x)
 *     + a (fmpz_t) - the unit p^(-v) x
 *
 *************************************************************************
 */

int fmpz_valuation_unit(fmpz_t a, const fmpz_t x, const fmpz_t p);

/*************************************************************************
 *
 * Function: fmpz_valuation
 *
 * Description: given an integer x and a prime p, return v = v_p(x) -
 *              the valuation of x at p.
 *
 * Arguments:
 *     + x (const fmpz_t) - the integer whose valuation is to be computed
 *     + p (const fmpz_t) - the prime
 *
 * Returns:
 *     + (int) - the valuation v = v_p(x)
 *
 *************************************************************************
 */

int fmpz_valuation(const fmpz_t x, const fmpz_t p);

/*************************************************************************
 *
 * Function: fmpq_valuation
 *
 * Description: given a rational x and a prime p, return v = v_p(x) -
 *              the valuation of x at p.
 *
 * Arguments:
 *     + x (const fmpq_t) - the number whose valuation is to be computed
 *     + p (const fmpz_t) - the prime
 *
 * Returns:
 *     + (int) - the valuation v = v_p(x)
 *
 *************************************************************************
 */

int fmpq_valuation(const fmpq_t x, const fmpz_t p);

/*************************************************************************
 *
 * Function: hilbert_symbol
 *
 * Description: given two integers x,y and a prime p, return (x,y)_p -
 *              the Hilbert symbol of (x,y) at p.
 *
 * Arguments:
 *     + x,y (const fmpz_t) - the integers whose symbol is to be computed
 *     + p (const fmpz_t) - the prime
 *
 * Returns:
 *     + (int) - the Hilbert symbol (x,y)_p
 *
 *************************************************************************
 */

int hilbert_symbol(const fmpz_t x, const fmpz_t y, const fmpz_t p);

/* Recursive function for a temporary extended Euclidean algorithm. */
/* It uses pointers to return multiple values. */

/*************************************************************************
 *
 * Function: gcdext
 *
 * Description: given two integers a,b, compute d = gcd(a,b) and
 *              two integers x,y such that ax+by = d.
 *
 * Arguments:
 *     + a,b (Z64) - the integers whose gcd is to be computed
 *
 * Returns:
 *     + (Z64) - the gcd d = gcd(a,b)
 *     + x,y (Z64) - integers such that ax+by = d
 *
 *************************************************************************
 */

Z64 gcdext(Z64 a, Z64 b, Z64 *x, Z64 *y);

/*************************************************************************
 *
 * Function: lcm
 *
 * Description: given two integers a,b, compute their lcm -
 *              least common multiple.
 *
 * Arguments:
 *     + a,b (int) - the integers whose lcm is to be computed
 *
 * Returns:
 *     + (Z64) - the lcm m = lcm(a,b)
 *
 *************************************************************************
 */

int lcm(int a, int b);

/*************************************************************************
 *
 * Function: rational_lt
 *
 * Description: Compares two rationals a,b (rational is a type from CARAT)
 *
 * Arguments:
 *     + a,b (rational) - the numbers to be compared.
 *
 * Returns:
 *     + (int) - 1 if a < b
 *               0 if a >= b
 *
 *************************************************************************
 */

int rational_lt(rational a, rational b);

/*************************************************************************
 *
 * Function: rational_sum
 *
 * Description: Sums two rationals a,b (rational is a type from CARAT).
 *
 * Arguments:
 *     + a,b (rational) - the numbers to be summed
 *
 * Returns:
 *     + (rational) - a + b
 *
 *************************************************************************
 */

rational rational_sum(rational a, rational b);

/*************************************************************************
 *
 * Function: bernoulli_number
 *
 * Description: Compute the n-th Bernoulli number B_n
 *
 * Arguments:
 *     + n (ulong) - the index of the Bernoulli number to be computed
 *
 * Returns:
 *     + b (fmpq_t) - b = B_n is the n-th Bernoulli number
 *
 *************************************************************************
 */

void bernoulli_number(fmpq_t b, ulong n);

/*************************************************************************
 *
 * Function: binomial_coefficient
 *
 * Description: Computes the binomial coefficient (n choose k)
 *
 * Arguments:
 *     + n (fmpz_t) - number of objects to choose from
 *     + k (fmpz_t) - number of objects to choose
 *
 * Returns:
 *     + b (fmpz_t) - b = (n choose k) is the binomial coefficient.
 *
 *************************************************************************
 */

void binomial_coefficient(fmpz_t b, const fmpz_t n, const fmpz_t k);

/*************************************************************************
 *
 * Function: kronecker_symbol
 *
 * Description: Computes the Kronecker symbol (a|n)
 *
 * Arguments:
 *     + a (fmpz_t) - top of symbol
 *     + n (fmpz_t) - bottom of symbol
 *
 * Returns:
 *     + (slong) - the Kronecker symbol (a|n)
 *
 *************************************************************************
 */

slong kronecker_symbol(const fmpz_t a, const fmpz_t n);

/*************************************************************************
 *
 * Function: bernoulli_number_chi
 *
 * Description: Compute the n-th Bernoulli number B_{n,chi} where
 *              chi is the quadratic character corresponding to the
 *              quadratic field with discrminant d.
 *
 * Arguments:
 *     + n (ulong) - the index of the Bernoulli number to be computed
 *     + d (const fmpz_t) - specifies that chi = chi_d
 *
 * Returns:
 *     + b_chi (fmpq_t) - b_chi = B_{n,chi} is the n-th Bernoulli number.
 *
 *************************************************************************
 */

void bernoulli_number_chi(fmpq_t b_chi, ulong n, const fmpz_t d);

/*************************************************************************
 *
 * Function: fmpz_is_local_square
 *
 * Description: Checks if an integer is a square in Z_p.
 *
 * Arguments:
 *     + a (const fmpz_t) - the number to be checked
 *     + p (const fmpz_t) - the prime at which to localize
 *
 * Returns:
 *     + (bool) - true if a is a square in Z_p (exists x^2 = a)
 *
 *************************************************************************
 */

bool fmpz_is_local_square(const fmpz_t a, const fmpz_t p);

/*************************************************************************
 *
 * Function: fmpq_is_local_square
 *
 * Description: Checks if a rational is a square in Q_p.
 *
 * Arguments:
 *     + a (const fmpq_t) - the number to be checked
 *     + p (const fmpz_t) - the prime at which to localize
 *
 * Returns:
 *     + (bool) - true if a is a square in Q_p (exists x^2 = a)
 *
 *************************************************************************
 */

bool fmpq_is_local_square(const fmpq_t a, const fmpz_t p);

/*************************************************************************
 *
 * Function: fq_is_square
 *
 * Description: Checks if a finite field element is a square.
 *              Assuming the finite field F is a prime field F_p.
 *
 * Arguments:
 *     + a (const fq_t) - the element to be checked
 *     + F (const fq_ctx_t) - the finite field
 *
 * Returns:
 *     + (bool) - true if a is a square in F (exists x^2 = a)
 *
 *************************************************************************
 */

bool fq_is_square(const fq_t a, const fq_ctx_t F);

// After writing this function we found its equivalent in flint,
// hence we don't need to use it - finding a square root in a finite field.

// void fq_sqrt(fq_t sqrt_a, const fq_t a, const fq_ctx_t F);

/*************************************************************************
 *
 * Function: fmpq_is_integral
 *
 * Description: Checks if a ratinoal number is integral.
 *
 * Arguments:
 *     + r (const fmpq_t) - the number to be checked
 *
 * Returns:
 *     + (bool) - true if r is integral
 *
 *************************************************************************
 */

bool fmpq_is_integral(const fmpq_t r);

/*************************************************************************
 *
 * Function: fmpq_floor
 *
 * Description: Returns the floor of a rational number - [x].
 *
 * Arguments:
 *     + x (const fmpq_t) - the input number
 *
 * Returns:
 *     + res (fmpz_t) - floor of x = [x].
 *
 *************************************************************************
 */

void fmpq_floor(fmpz_t res, const fmpq_t x);

/*************************************************************************
 *
 * Function: primes_up_to
 *
 * Description: Generate the list of prime numbers up to a bound.
 *
 * Arguments:
 *     + bound (int) - the upper bound up to which we wish to compute
 *
 * Returns:
 *     + ps (int**) - the list of primes up to bound.
 *     + (int) - the length of the list (pi(bound))
 *
 *************************************************************************
 */

int primes_up_to(int** ps, int bound);

/*************************************************************************
 *
 * Function: primes_up_to_prime_to
 *
 * Description: Generate the list of prime numbers up to a bound,
 *              which do not divide some bad locus.
 *
 * Arguments:
 *     + bound (int) - the upper bound up to which we wish to compute
 *     + bad (int) - product of bad primes, to be avoided.
 *
 * Returns:
 *     + ps (int**) - the list of primes up to bound, prime to bad.
 *     + (int) - the length of the list (pi(bound))
 *
 *************************************************************************
 */

int primes_up_to_prime_to(int** ps, int bound, int bad);

#endif // __ARITH_H
