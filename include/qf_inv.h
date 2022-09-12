#ifndef __QF_INV_H__
#define __QF_INV_H__

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

// a data structure for invariants of a quadratic form

typedef struct {
  slong num_stored;
  slong num_alloc;

  ulong n;
  slong I;
  fmpz_t* primes;
  slong* symbols;
} qf_inv;

typedef qf_inv qf_inv_t[1];

void qf_inv_init(qf_inv_t hasse, ulong n);

void qf_inv_clear(qf_inv_t hasse);

void qf_inv_add(qf_inv_t hasse, const fmpz_t p, slong s);

slong hasse(fmpz_mat_t D, fmpz_t p);

void invariants(fmpz_t prod, qf_inv_t F, const fmpz_mat_t Q);

int witt_to_hasse(fmpz_mat_t hasse,
		  const fmpz_t det,
		  const qf_inv_t finite);

#endif // __QF_INV_H__
