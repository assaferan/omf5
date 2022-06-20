#ifndef __GENUS_H__
#define __GENUS_H__

#include <flint/fmpq.h>

#include "carat/matrix.h"

#include "hash.h"

/* compute the genus of a quadratic form */
hash_table* get_genus_reps(matrix_TYP* Q);

#endif // __GENUS_H__
