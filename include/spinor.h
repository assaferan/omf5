#ifndef __SPINOR_H__
#define __SPINOR_H__

#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>

#include "isometry.h"
#include "typedefs.h"

typedef struct {

  fmpz_mat_t Q;
  
  nmod_mat_t* rads;
  nmod_t* primes;
  slong* pivots;
  slong* exps; // the powers to which the primes divide the discriminant
  fmpz_mat_t* lifts; // vectors spanning the radicals modulo these powers of the primes
  // The primes are those dividing the discriminant, except for the prime 2
  // when the radical of the lattice modulo 2 is the whole lattice modulo 2.

  slong num_primes;
  W64 twist;

  // Whether the radical modulo 2 is of dimension 3 and the quadratic form is
  // odd at a vector of it. The character at 2 is then that of the spinor norm.
  bool odd_on_rad_2;

  
} spinor;

typedef spinor spinor_t[1];

void spinor_init_fmpz(spinor_t spinor, const fmpz_t disc);
void spinor_init_fmpz_mat(spinor_t spinor, const fmpz_mat_t q);
void spinor_init_square_matrix(spinor_t spinor, const square_matrix_t q);
void spinor_clear(spinor_t spinor);

bool spinor_has_prime(const spinor_t spinor, ulong p);

W64 spinor_norm_fmpz_mat(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom);

W64 spinor_norm(const spinor_t spinor, matrix_TYP* mat, int denom);
W64 spinor_norm_isom(const spinor_t spinor, const isometry_t isom);
W64 spinor_norm_cd_fmpz_mat(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom);
W64 spinor_norm_zas_fmpz_mat(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom);

#endif // __SPINOR_H__
