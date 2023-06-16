/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : arith.c
 *
 * Description: simple useful arithmetic functions
 *
 **********************************************************
 */

// System dependencies

#include <assert.h>

// Self dependencies

#include "arith.h"

// Lookup tables for quickly computing the Hilbert symbol

// When p is odd, only need v_p(x) mod 2 and whether u_p(x) is a square

int hilbert_lut_odd[16] = { 1, 1, 1, 1,
			    1, 1,-1,-1,
			    1,-1, 1,-1,
			    1,-1,-1, 1 };

// When p==2 we need to consider each of the inputs mod 8

int hilbert_lut_p2[64] = { 1, 1, 1, 1, 1, 1, 1, 1,
			   1,-1, 1,-1,-1, 1,-1, 1,
			   1, 1, 1, 1,-1,-1,-1,-1,
			   1,-1, 1,-1, 1,-1, 1,-1,
			   1,-1,-1, 1, 1,-1,-1, 1,
			   1, 1,-1,-1,-1,-1, 1, 1,
			   1,-1,-1, 1,-1, 1, 1,-1,
			   1, 1,-1,-1, 1, 1,-1,-1 };

// return the valuation of x at p, and put the non-divisible part in a.
int fmpz_valuation_unit(fmpz_t a, const fmpz_t x, const fmpz_t p)
{
  int a_val = 0;
  fmpz_t a_div_p, a_mod_p;
  
  fmpz_set(a, x);
  fmpz_init(a_div_p);
  fmpz_init(a_mod_p);

  // We repeatedly divide by p, until we obtain the unit
  
  fmpz_fdiv_qr(a_div_p, a_mod_p, a, p);
  while (fmpz_is_zero(a_mod_p)) {
    ++a_val;
    fmpz_set(a, a_div_p);
    fmpz_fdiv_qr(a_div_p, a_mod_p, a, p);
  }
  
  fmpz_clear(a_div_p);
  fmpz_clear(a_mod_p);

  return a_val;
}

// return the valuation of x at p
int fmpz_valuation(const fmpz_t x, const fmpz_t p)
{
  fmpz_t a;
  int a_val;

  // If x = 0, we would need to return infinity...
  assert(!fmpz_is_zero(x));
  fmpz_init(a);

  // simply apply previous function and forget the unit 
  a_val = fmpz_valuation_unit(a,x,p);
  
  fmpz_clear(a);

  return a_val;
}

// return the valuation of x at p
int fmpq_valuation(const fmpq_t x, const fmpz_t p)
{
  // If x = 0, we would need to return infinity...
  assert(!fmpq_is_zero(x));
  return fmpz_valuation(fmpq_numref(x), p) - fmpz_valuation(fmpq_denref(x), p);
}

// return the Hilbert symbol (x,y)_p
int hilbert_symbol(const fmpz_t x, const fmpz_t y, const fmpz_t p)
{
  fmpz_t a, b;
  int a_val, b_val;
  int aa, bb, index, a_notsqr, b_notsqr, hilbert;
  fmpz_t p_mod_4, a_mod_8, b_mod_8;

  fmpz_init(a);
  fmpz_init(b);
  
  a_val = fmpz_valuation_unit(a, x, p);
  b_val = fmpz_valuation_unit(b, y, p);

  if (fmpz_equal_si(p,2)) {
    fmpz_init(a_mod_8);
    fmpz_init(b_mod_8);
    fmpz_mod_ui(a_mod_8, a, 8);
    fmpz_mod_ui(b_mod_8, b, 8);
    aa = (fmpz_get_si(a_mod_8) >> 1) & 0x3;
    bb = (fmpz_get_si(b_mod_8) >> 1) & 0x3;
    index = ((a_val&0x1)<<5) | (aa << 3) | ((b_val&0x1)<<2) | bb;
    fmpz_clear(a_mod_8);
    fmpz_clear(b_mod_8);
    fmpz_clear(a);
    fmpz_clear(b);
    return hilbert_lut_p2[index];
  }

  // Using fmpz_jacobi is better, as it handles more endcases.
  /* a_notsqr = mpz_legendre(a.get_mpz_t(), p.get_mpz_t()) == -1; */
  /* b_notsqr = mpz_legendre(b.get_mpz_t(), p.get_mpz_t()) == -1; */
  a_notsqr = (fmpz_jacobi(a, p) == -1);
  b_notsqr = (fmpz_jacobi(b, p) == -1);

  fmpz_init(p_mod_4);
  index = ((a_val&0x1)<<3) | (a_notsqr<<2) | ((b_val&0x1)<<1) | b_notsqr;
  fmpz_mod_ui(p_mod_4, p, 4);
  if (((index & 0xa) == 0xa) && (fmpz_equal_si(p_mod_4,0x3))) {
    hilbert = -hilbert_lut_odd[index];
  }
  else {
    hilbert = hilbert_lut_odd[index];
  }

  fmpz_clear(p_mod_4);
  fmpz_clear(a);
  fmpz_clear(b);
  return hilbert;
}

// lcm for ints
// We simply use fmpz lcm, wrapping the conversions

int lcm(int a, int b)
{
  fmpz_t a_Z, b_Z, m_Z;
  int m;

  fmpz_init_set_si(a_Z, a);
  fmpz_init_set_si(b_Z, b);
  fmpz_init(m_Z);
  fmpz_lcm(m_Z, a_Z, b_Z);
  m = fmpz_get_si(m_Z);

  fmpz_clear(a_Z);
  fmpz_clear(b_Z);
  fmpz_clear(m_Z);
  return m;
}

/* Recursive function for a temporary extended Euclidean algorithm. */
/* It uses pointers to return multiple values. */
Z64 gcdext(Z64 a, Z64 b, Z64 *x, Z64 *y)
{
  Z64 _x, _y, gcd;
  
  if (a == 0) {
    *x = 0;
    *y = 1;
    return b;
  }

  gcd = gcdext(b % a, a, &_x, &_y);
 
  *x = _y - (b/a) * _x;
  *y = _x;

  /* we return the positive gcd */
  if (gcd < 0) {
    *x = -*x;
    *y = -*y;
    gcd = -gcd;
  }
  
  return gcd;
}

// return a < b for the CARAT type rational

int rational_lt(rational a, rational b)
{
  // we assume here that a and b are positive as they should be
  return a.z * b.n < b.z * a.n;
}

// return a + b for the CARAT type rational

rational rational_sum(rational a, rational b)
{
  Z64 d, dummy1, dummy2;
  rational c;

  c.z = a.z * b.n + b.z * a.n;
  c.n = a.n * b.n;

  d = gcdext(c.z, c.n, &dummy1, &dummy2);
  c.z /= d;
  c.n /= d;

  return c;
}

// We need twisted bernoulli numbers, so here is our implementation of it

// collect all of the numbers along the way.
// We are using Akiyama and Tanigawa's algorithm
// It's not the fastest, but it is one of the simplest.

// This should be the same as FLINT's arith_bernoulli_number_vec, maybe replace
void _bernoulli_up_to(fmpq_t* b, ulong n)
{
  fmpq_t* a;
  fmpz_t mult;
  fmpq_t diff;
  ulong i, j;

  fmpz_init(mult);
  fmpq_init(diff);
  
  a = (fmpq_t*)malloc((n+1)*sizeof(fmpq_t));
#ifdef DEBUG
  if (a == NULL) {
    printf("could not allocate vector of size %lu\n", n);
    exit(-1);
  }
  if (b == NULL) {
    printf("received a null pointer as argument!\n");
    exit(-1);
  }
#endif // DEBUG

  for (i = 0; i <= n; i++) {
    fmpq_init(a[i]);
    fmpq_set_si(a[i], 1, i+1);
  }
  fmpq_set(b[0], a[0]);

  for (i = 0; i < n; i++) {
    for (j = 0; j < n - i; j++) {
      fmpz_set_si(mult, j+1);
      fmpq_sub(diff, a[j], a[j+1]);
      fmpq_mul_fmpz(a[j], diff, mult);
    }
    fmpq_set(b[i+1],a[0]);
  }

  fmpz_clear(mult);
  fmpq_clear(diff);
  for (i = 0; i < n; i++)
    fmpq_clear(a[i]);
  free(a);
  return;
  
}

// Compute B_n - simply calls _bernoulli_up_to and pulls the last entry.

void bernoulli_number(fmpq_t b, ulong n)
{
  fmpq_t* b_vec;
  int i;

  b_vec = (fmpq_t*)malloc((n+1)*sizeof(fmpq_t));
  for (i = 0; i <= n; i++)
    fmpq_init(b_vec[i]);
  
  _bernoulli_up_to(b_vec, n);

  fmpq_set(b, b_vec[n]);

  for (i = 0; i <= n; i++)
    fmpq_clear(b_vec[i]);
  free(b_vec);
  
  return;
}

// Compute n choose k = n!/k!(n-k)!

void binomial_coefficient(fmpz_t b, const fmpz_t n, const fmpz_t k)
{
  fmpz_t i, n_k;
  
  fmpz_init(i); // this sets to 0
  fmpz_init(n_k);
  fmpz_set_si(b,1);

  fmpz_sub(n_k, n, k);
  // use symmetry if n < 2k
  if (fmpz_cmp(n_k, k) < 0)
    return binomial_coefficient(b,n,n_k);

  for (; fmpz_cmp(i,k) < 0; ) {
    fmpz_sub(n_k, n, i);
    // multiply by (n-k-i)
    fmpz_mul(b, b, n_k);
    fmpz_add_si(i,i,1);
    // divide by (i+1)
    fmpz_fdiv_q(b, b, i);
  }

  fmpz_clear(n_k);
  fmpz_clear(i);
  return;
}

// Compute the bernoulli polynomial, for twisted Bernoulli numbers

void _bernoulli_poly(fmpq_t* b, ulong n)
{
  fmpq_t* b_vec;
  ulong k;
  fmpz_t nn, kk, n_choose_k;

  fmpz_init_set_ui(nn, n);
  fmpz_init(kk);
  fmpz_init(n_choose_k);

  b_vec = (fmpq_t*)malloc((n+1)*sizeof(fmpq_t));
  for (k = 0; k <= n; k++)
    fmpq_init(b_vec[k]);
  _bernoulli_up_to(b_vec, n);
  
  for (k = 0; k <= n; k++)
    fmpq_set(b[k], b_vec[n-k]);
  
  for (k = 0; k <= n; k++) {
    fmpz_set_ui(kk,k);
    binomial_coefficient(n_choose_k, nn, kk);
    fmpq_mul_fmpz(b[k], b[k], n_choose_k);
  }

  for (k = 0; k <= n; k++)
    fmpq_clear(b_vec[k]);
  free(b_vec);
  fmpz_clear(n_choose_k);
  fmpz_clear(kk);
  fmpz_clear(nn);
  return;
}

// Compute the Kronecker symbol (a|n)
// TODO - could be spurious to FLINT's implementation of Jacobi symbol
// Current implementaion is recursive

slong kronecker_symbol(const fmpz_t a, const fmpz_t n)
{
  fmpz_t val_mpz, n_sqr, minus1, minus_n, minus_a, n_star, two, a_mod_n;
  slong val, n_prime, res;

  // extremal cases
  if (fmpz_is_zero(n)) return (fmpz_is_pm1(a) ? 1 : 0);
  if (fmpz_equal_si(n,-1)) return ((fmpz_cmp_si(a,0)) < 0 ? -1 : 1);
  if (fmpz_is_one(n)) return 1;
  if (fmpz_equal_si(n,2)) {
    if (fmpz_is_even(a)) return 0;
    fmpz_init(val_mpz);
    fmpz_mod_ui(val_mpz, a, 8);
    val = fmpz_get_si(val_mpz);
    fmpz_clear(val_mpz);
    val /= 2;
    if ((val == 0) || (val == 3))
      return 1;
    return -1;
  }
  if (fmpz_equal_si(a, -1)) {
    n_prime = fmpz_get_si(n);
    while (n_prime % 2 == 0) n_prime /= 2;
    return ((n_prime / 2) % 2 == 0) ? 1 : -1;
  }
  if (fmpz_equal_si(a,2) && fmpz_is_odd(n)) {
    fmpz_init(n_sqr);
    fmpz_pow_ui(n_sqr, n, 2);
    fmpz_fdiv_q_2exp(n_sqr, n_sqr, 3);
    res = (fmpz_is_even(n_sqr) ? 1 : -1);
    fmpz_clear(n_sqr);
    return res;
  }
  // multiplicativity
  if (fmpz_cmp_si(n,0) < 0) {
    fmpz_init_set_si(minus1, -1);
    fmpz_init(minus_n);
    fmpz_neg(minus_n, n);
    res = kronecker_symbol(a, minus1) * kronecker_symbol(a, minus_n);
    fmpz_clear(minus_n);
    fmpz_clear(minus1);
    return res;
  }
  if (fmpz_cmp_si(a,0) < 0) {
    fmpz_init_set_si(minus1, -1);
    fmpz_init(minus_a);
    fmpz_neg(minus_a, a);
    res = kronecker_symbol(minus1, n) * kronecker_symbol(minus_a, n);
    fmpz_clear(minus_a);
    fmpz_clear(minus1);
    return res;
  }

  // now may assume n >= 3, a >= 0
 
  // quadratic reciprocity
  if (fmpz_cmp(a, n) < 0) {
    n_prime = fmpz_get_si(n);
    while (n_prime % 2 == 0) n_prime /= 2;
    fmpz_init_set(n_star, n);
    if ((n_prime / 2) % 2 != 0)
      fmpz_neg(n_star, n);
    res = kronecker_symbol(n_star, a);
    fmpz_clear(n_star);
    return res;
  }

  // now we may also assume a ge n

  // if n = 2 mod 4, we can't reduce, use multiplicativity again
  if (fmpz_get_si(n) % 4 == 2) {
    fmpz_init(n_star);
    fmpz_init_set_si(two,2);
    fmpz_fdiv_q_2exp(n_star, n, 1);
    res = kronecker_symbol(a,n_star)*kronecker_symbol(a, two);
    fmpz_clear(two);
    fmpz_clear(n_star);
    return res;
  }
  // now we can reduce
  fmpz_init(a_mod_n);
  fmpz_mod(a_mod_n, a, n);
  res = kronecker_symbol(a_mod_n, n);
  fmpz_clear(a_mod_n);
  return res;
}


// B_{n, chi} where chi is the quadratic character corresponding to
// quadratic field with discrminant d (Is it? verify we are working
// with the correct discriminant (d or 4d maybe?))
void bernoulli_number_chi(fmpq_t b_chi, ulong n, const fmpz_t d)
{
  fmpq_t* b;
  fmpz_t a_pow, d_pow, a, nn;
  fmpq_t s, ba_pow;
  slong chi_a;
  ulong k;

  fmpz_init(a);
  fmpz_init(a_pow);
  fmpz_init(d_pow);
  fmpq_init(ba_pow);
  fmpz_init_set_ui(nn, n);
  fmpq_init(s);
  b = (fmpq_t*)malloc((n+1)*sizeof(fmpq_t));
  for (k = 0; k <= n; k++)
    fmpq_init(b[k]);
  _bernoulli_poly(b,n);

  fmpz_pow_ui(d_pow, d, n);

  fmpq_zero(b_chi);

  for (fmpz_zero(a); fmpz_cmp(a,d) < 0; fmpz_add_si(a,a,1)) {
    chi_a = kronecker_symbol(a, nn);
    fmpz_one(a_pow);
    fmpq_zero(s);
    for (k = 0; k <= n; k++) {
      fmpz_fdiv_q(d_pow, d_pow, d);
      fmpq_mul_fmpz(ba_pow, b[k], a_pow);
      fmpq_mul_fmpz(ba_pow, ba_pow, d_pow);
      fmpq_add(s, s, ba_pow);
      fmpz_mul(a_pow, a_pow, a);
    }
    fmpq_mul_si(s, s, chi_a);
    fmpq_add(b_chi, b_chi, s);
  }
  for (k = 0; k <= n; k++)
    fmpq_clear(b[k]);
  free(b);
  fmpq_clear(s);
  fmpz_clear(nn);
  fmpq_clear(ba_pow);
  fmpz_clear(d_pow);
  fmpz_clear(a_pow);
  fmpz_clear(a);
  return;
}

// check if a is a local square
bool fmpz_is_local_square(const fmpz_t a, const fmpz_t p)
{
  fmpz_t a0, a0_1;
  ulong v, i, w, ee;
  slong ok, ww;
  int res;

  if (fmpz_is_zero(a)) return true;
  v = fmpz_valuation(a,p);
  if (v % 2 == 1) return false;
  
  fmpz_init_set(a0,a);
  for (i = 0; i < v; i++)
    fmpz_fdiv_q(a0, a0, p);

  ok = (kronecker_symbol(a0,p) != -1);

  if (!fmpz_equal_si(p,2)) {
    fmpz_clear(a0);
    return ok;
  }
  fmpz_init(a0_1);
  fmpz_sub_si(a0_1, a0, 1);
  w = fmpz_valuation(a0_1,p);

  assert(w >= 1);
  ee = 2;

  while ((w < ee) && (w % 2 == 0)) {
    ww = (1+ (1<< (w/2)))*(1+ (1<< (w/2)));
    fmpz_fdiv_q_si(a0, a0, ww);
    fmpz_sub_si(a0_1, a0, 1);
    w = fmpz_valuation(a0_1, p);
  }
  fmpz_mod_ui(a0, a0, 8);
  res = ((w > ee) || ((w == ee) && fmpz_is_one(a0)));
  fmpz_clear(a0);
  fmpz_clear(a0_1);
  return res;
}

// check if a is a local square
bool fmpq_is_local_square(const fmpq_t a, const fmpz_t p)
{
  // if a = x/y, then a is a square if and only if either both x and y are squares or both are nonsquares
  return (fmpz_is_local_square(fmpq_numref(a), p)
	  == fmpz_is_local_square(fmpq_denref(a), p));
}

// Checks if a is a square in F, assumes F = F_p is a prime field
bool fq_is_square(const fq_t a, const fq_ctx_t F)
{
  fmpz_t a_lift;
  slong legendre;
  
  fmpz_init(a_lift);
  // Here we use the assumption, the trace is simply a lift of the element to the integers
  fq_trace(a_lift, a, F);
  legendre = kronecker_symbol(a_lift, fq_ctx_prime(F));

  fmpz_clear(a_lift);
  return (legendre == -1) ? false : true;
}

// Checks if r is integral
bool fmpq_is_integral(const fmpq_t r)
{
  return fmpz_divisible(fmpq_numref(r), fmpq_denref(r));
}

// returns the floor of r
void fmpq_floor(fmpz_t res, const fmpq_t r)
{
  fmpz_fdiv_q(res, fmpq_numref(r), fmpq_denref(r));
  return;
}

// Return prime up to bound
int primes_up_to(int** ps, int bound)
{
  int num_ps, p, i;
  int* sieve;
  
  num_ps = 0;
  sieve = malloc((bound+2) * sizeof(int));
  *ps = malloc(bound * sizeof(int));
  // since bound is small we construct the first primes using aristothenes sieve
  // replace by fmpz_nextprime if improving.
  for (i = 0; i <= bound+1; i++)
    sieve[i] = i;
  
  p = 2;
  while(p <= bound) {
    (*ps)[num_ps] = p;
    num_ps++;
    for (i = p; i <= bound; i+= p)
      sieve[i] = -1;
    while(sieve[p] == -1) p++;
  }

  free(sieve);
  return num_ps;
}

// Return prime up to bound, prime to bad
int primes_up_to_prime_to(int** ps, int bound, int bad)
{
  int num_ps, p, i;
  int* sieve;
  
  num_ps = 0;
  sieve = malloc((bound+2) * sizeof(int));
  *ps = malloc(bound * sizeof(int));
  // since bound is small we construct the first primes using aristothenes sieve
  // replace by fmpz_nextprime if improving.
  for (i = 0; i <= bound+1; i++)
    sieve[i] = i;
  
  p = 2;
  while(p <= bound) {
    if (bad % p != 0) {
      (*ps)[num_ps] = p;
      num_ps++;
    }
    for (i = p; i <= bound; i+= p)
      sieve[i] = -1;
    while(sieve[p] == -1) p++;
  }

  free(sieve);
  return num_ps;
}

// fq_sqrt seems to already exist

/* void fq_sqrt(fq_t sqrt_a, const fq_t a, const fq_ctx_t F) */
/* { */
/*   fmpz_t p, q; */
/*   fq_t z, b, c, r, t, t1; */
/*   int i, j, m, s, e; */
  
/*   if (fq_is_one(a,F)) { */
/*     fq_one(sqrt_a, F); */
/*     return; */
/*   } */
/*   if (fq_is_zero(a,F)) { */
/*     fq_zero(sqrt_a, F); */
/*     return; */
/*   } */

/*   assert(is_square(a,F)); */

/*   fmpz_init(p); */
/*   fmpz_init(q); */

/*   fmpz_set(p,fq_ctx_prime(F)); */
/*   fmpz_sub_si(q,p,1); */
  
/*   s = 0; */
/*   while (fmpz_is_even(q)) { */
/*     fmpz_div_q_2exp(q, q, 1); */
/*     s++; */
/*   } */

/*   if (s == 1) { */
/*     fmpz_add_si(q, p, 1); */
/*     fmpz_div_q_2exp(q, q, 4); */
/*     fq_pow(sqrt_a, a, q, F); */
/*     fmpz_clear(q); */
/*     fmpz_clear(p); */
/*     return; */
/*   } */
  
/*   fq_init(z,F); */
/*   fq_init(b,F); */
/*   fq_init(c,F); */
/*   fq_init(r,F); */
/*   fq_init(t,F); */
/*   fq_init(t1,F); */

/*   // looking for a non-square */
/*   fq_zero(z,F); */
/*   while (fq_is_square(z,F)) */
/*     fq_sub_one(z,z,F); */

/*   m = s; */
/*   fq_pow(c,z,q,F); */
/*   fq_pow(t,a,q,F); */
/*   fmpz_add_si(q,q,1); */
/*   fmpz_div_q_2exp(q,q,1); */
/*   fq_pow(r,a,q,F); */

/*   // !! TODO - this is poorly written should have the is_one condition in the outer while loop */
/*   while (true) { */
/*     if (fq_is_one(t,F)) { */
/*       fq_set(sqrt_a, r, F); */
/*       fq_clear(t1,F); */
/*       fq_clear(t,F); */
/*       fq_clear(r,F); */
/*       fq_clear(c,F); */
/*       fq_clear(b,F); */
/*       fq_clear(z,F); */
/*       fmpz_clear(s); */
/*       fmpz_clear(q); */
/*       fmpz_clear(p); */
/*       return; */
/*     } */
/*     i = 0; */
/*     fq_set(t1, t, F); */
/*     while (!fq_is_one(t1,F)) { */
/*       fq_sqr(t1,t1,F); */
/*       i++; */
/*     } */

/*     e = 1; */
/*     for (j = 0; j < m-i-1; j++) */
/*       e <<= 1; */

/*     fq_pow(b,c,e,F); */
/*     fq_mul(r,r,b,F); */
/*     fq_sqr(c,b,F); */
/*     fq_mul(t,t,c,F); */
/*     m = i; */
/*   } */

/*   f1_clear(t1,F); */
/*   fq_clear(t,F); */
/*   fq_clear(r,F); */
/*   fq_clear(c,F); */
/*   fq_clear(b,F); */
/*   fq_clear(z,F); */
/*   fmpz_clear(s); */
/*   fmpz_clear(q); */
/*   fmpz_clear(p); */

/*   assert(false); */
  
/*   return; */
/* } */
