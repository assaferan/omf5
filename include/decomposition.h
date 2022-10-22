#ifndef __DECOMPOSITION_H__
#define __DECOMPOSITION_H__

#include <flint/fmpq_mat.h>

#include "genus.h"

// this struct holds data for the decomposition of the spaces as Hecke modules

typedef struct {
  fmpq_mat_t** bases; // bases for the vector spaces, by conductor
  fmpq_mat_t** hecke;
  slong* num; // number of elements in the decomposition, by conductor
  slong num_primes;
  slong num_conductors;
} decomposition_struct;

typedef decomposition_struct decomposition_t[1];

void decomposition_init(decomposition_t decomp, slong num);
void decomposition_clear(decomposition_t decomp);
void decompose(decomposition_t decomp, const genus_t genus, slong c);

#endif // __DECOMPOSITION_H__
