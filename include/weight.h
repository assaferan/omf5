#ifndef __WEIGHT_H__
#define __WEIGHT_H__

#include <flint/fmpq_mat.h>

#include "aut_grp.h"
#include "genus.h"
#include "square_matrix.h"

isometry_t** hecke_matrices_isometries(const genus_t genus, int p);

void invariant_subspace(fmpq_mat_t sub, const aut_grp_t grp);

#endif // __WEIGHT_H__
