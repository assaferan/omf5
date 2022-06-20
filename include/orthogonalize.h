#ifndef __ORTHOGONALIZE_H__
#define  __ORTHOGONALIZE_H__

#include <flint/fmpz_mat.h>

void orthogonalize_gram(fmpz_mat_t D, const fmpz_mat_t Q);

#endif //  __ORTHOGONALIZE_H__
