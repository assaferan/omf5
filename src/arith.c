#include "arith.h"

/* Recursive function for a temporary extended Euclidean algorithm. */
/* It uses pointers to return multiple values. */
int gcdext(int a, int b, int *x, int *y)
{
  int _x, _y, gcd;
  
  if (a == 0) {
    *x = 0;
    *y = 1;
    return b;
  }

  gcd = gcdext(b % a, a, &_x, &_y);
 
  *x = _y - (b/a) * _x;
  *y = _x;

  /* we return the positive gcd */
  if (gcd < 0) {
    *x = -*x;
    *y = -*y;
    gcd = -gcd;
  }
  
  return gcd;
}
