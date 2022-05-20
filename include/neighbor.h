#ifndef __NEIGHBOR_H__
#define __NEIGHBOR_H__

#include "carat/matrix.h"

/* Compute one p-neighbour for Q_orig corresponding to vector x 
 * On error, return NULL.
*/
matrix_TYP* q61_nb(matrix_TYP* Q_orig, int p, matrix_TYP* x_mat);

#endif // __NEIGHBOR_H__
