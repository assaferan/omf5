#include <stdio.h>
#include <time.h>

#include <carat/typedef.h>

#include "arith.h"
#include "genus.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "neighbor.h"

int get_ev(hash_table* genus, int p)
{
  int* a;
  int num, e;
  /* int Q_coeffs[15] = {2,1,2,0,0,2,0,0,0,4,1,0,0,-1,6}; */
  /* matrix_TYP* Q; */

  a = (int*)malloc(genus->num_stored * sizeof(int));
  for (num = 0; num < genus->num_stored; num++)
    a[num] = 0;
  
  clock_t cpuclock;
  double cputime;

  /* cpuclock = clock(); */

  /* Q = init_sym_matrix(Q_coeffs); */
  /* genus = get_genus_reps(Q); */

  /* cpuclock = clock() - cpuclock; */
  /* cputime = cpuclock / CLOCKS_PER_SEC; */
  /* printf("computing genus took %f\n", cputime); */
  
  cpuclock = clock();
  
  for (num = 0; num < p; num++) {
    q61_nbs1(a, p, num, NULL, genus);
  }

  cpuclock = clock() - cpuclock;
  cputime = cpuclock / CLOCKS_PER_SEC;

  /* Since we know the specific eigenform, it is hard coded here */
  /* e = -a[1]+a[2]+2*a[3]-2*a[5]+(3*a[4]-4*a[7])/6; */
  /* e = -a[2]+a[0]+2*a[7]-2*a[5]+(3*a[3]-4*a[1])/6; */
  e = -a[2]+a[0]+2*a[6]-2*a[5]+(3*a[3]-4*a[1])/6;
  
  printf("%4d %4d - %10d %10d %10d %10d %10d %10d %10d %10d - %10f\n",
    p, e, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], cputime);

  /* free_hash(genus); */
  /* free_mat(Q); */
  
  /* if ((p == 97) && (e != 453)) */
  /*   return -1; */
  
  /* return 0; */

  free(a);
  
  return e;
}

int test_61()
{
  int i;
  int ps[8] = {2,3,5,7,11,13,17,19};
  int evs[8] = {-7,-3,3,-9,-4,-3,37,-75};
  hash_table* genus;
  clock_t cpuclock;
  double cputime;
  int Q_coeffs[15] = {2,1,2,0,0,2,0,0,0,4,1,0,0,-1,6};
  matrix_TYP* Q;

  cpuclock = clock();

  Q = init_sym_matrix(Q_coeffs);
  genus = get_genus_reps(Q);

  cpuclock = clock() - cpuclock;
  cputime = cpuclock / CLOCKS_PER_SEC;
  
  printf("computing genus took %f\n", cputime);
  
  for (i = 0; i < 8; i++) {
    if (get_ev(genus, ps[i]) != evs[i]) {
      free_hash(genus);
      free_mat(Q);
      return FALSE;
    }
  }

  free_hash(genus);
  free_mat(Q);
  
  return TRUE;
}

int main(int argc, char* argv[])
{
  int test_res;
  /* int p; */
  
  /* printf("testing Tornaria in C...\n"); */
  /* if (argc > 2) */
  /*   return -1; */

  /* p = 97; */

  /* if (argc == 2) */
  /*   sscanf(argv[1], "%d", &p); */
  
  /* return total(p); */

  if (argc > 1)
    return -1;

  test_res = test_61();

  if (test_res == TRUE)
    return 0;
  else
    return -1;
  
}
