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

int rational_lt(rational a, rational b)
{
  // we assume here that a and b are positive as they should be
  return a.z * b.n < b.z * a.n;
}

rational rational_sum(rational a, rational b)
{
  int d, dummy1, dummy2;
  rational c;

  c.z = a.z * b.n + b.z * a.n;
  c.n = a.n * b.n;

  d = gcdext(c.z, c.n, &dummy1, &dummy2);
  c.z /= d;
  c.n /= d;

  return c;
}
