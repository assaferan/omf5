#ifndef __ARITH_H__
#define __ARITH_H__

#include "carat/typedef.h"

/* Recursive function for a temporary extended Euclidean algorithm. */
/* It uses pointers to return multiple values. */

int gcdext(int a, int b, int *x, int *y);

int rational_lt(rational a, rational b);

rational rational_sum(rational a, rational b);

#endif // __ARITH_H
