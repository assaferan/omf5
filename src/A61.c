#include <stdio.h>
#include <time.h>

#include "arith.h"
#include "genus.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "neighbor.h"

int total(int p)
{
  int a[8] = {0};
  int num, e;

  clock_t cpuclock;
  double cputime;

  cpuclock = clock();
  
  for (num = 0; num < p; num++) {
    q61_nbs1(a, p, num, NULL);
  }

  cpuclock = clock() - cpuclock;
  cputime = cpuclock / CLOCKS_PER_SEC;

  // Since we know the specific eigenform, it is hard coded here
  e = -a[1]+a[2]+2*a[3]-2*a[5]+(3*a[4]-4*a[7])/6;
  
  printf("%4d %4d - %10d %10d %10d %10d %10d %10d %10d %10d - %10f\n",
    p, e, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], cputime);

  return 0;
}

int main(int argc, char* argv[])
{
  int p;
  
  /* printf("testing Tornaria in C...\n"); */
  if (argc > 2)
    return -1;

  p = 97;

  if (argc == 2)
    sscanf(argv[1], "%d", &p);
  
  total(p);
  
  return 0;
}
