#ifndef __MASS_H__
#define __MASS_H__

#include <carat/matrix.h>

#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>

#include "jordan.h"
#include "square_matrix.h"

void local_factor(fmpq_t f, const fmpq_mat_t g, const fmpz_t p);

void diagonal_join(fmpq_mat_t joined, const jordan_data_t jordan);

/* compute the mass of the genus represented by q */
void get_mass(fmpq_t mass, const square_matrix_t q);

#endif // __MASS_H__
