#ifndef __SPINOR_H__
#define __SPINOR_H__

#include <flint/fq_nmod.h>
#include <flint/fq_nmod_mat.h>

#include "typedefs.h"

typedef struct {

  fq_nmod_mat_t* rads;
  fq_nmod_ctx_t* fields;

  slong num_primes;
  W64 twist;
  
} spinor;

typedef spinor spinor_t[1];

void spinor_init(spinor_t spinor, const fmpz_mat_t q);
void spinor_clear(spinor_t spinor);

W64 spinor_norm(const spinor_t spinor, const fmpz_mat_t mat, const fmpz_t denom);

#endif // __SPINOR_H__
