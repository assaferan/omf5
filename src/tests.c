#include <time.h>

#include "genus.h"
#include "hash.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "tests.h"

int test_61()
{
  int i;
#ifdef DEBUG
  int j;
#endif // DEBUG
  int ps[25] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
  int test_evs[25] = {-7, -3, 3, -9, -4, -3, 37, -75, 10, 212, -6, -88, -3, 547, -147, -108, -45,
		      145, -632, -650, 859, -978, 931, -571, 453};
  hash_table* genus;
  clock_t cpuclock_0, cpuclock_1, cpudiff;
  double cputime;
  int Q_coeffs[15] = {2,1,2,0,0,2,0,0,0,4,1,0,0,-1,6};
  matrix_TYP* Q, *hecke;
  eigenvalues* evs;
  nf_elem_t ev;

  cpuclock_0 = clock();

  Q = init_sym_matrix(Q_coeffs);
  genus = get_genus_reps(Q);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  printf("computing genus took %f\n", cputime);
  
  hecke = hecke_matrix(genus, 2);
  evs = get_eigenvalues(hecke);
  free_mat(hecke);

  cpuclock_0 = clock();
  cpudiff = cpuclock_0 - cpuclock_1;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  printf("computing eigenvectors took %f\n", cputime);
  
#ifdef DEBUG
  printf("eigenvectors are:\n");
  for (i = 0; i < evs->num; i++) {
    for (j = 0; j < evs->dim; j++) {
	nf_elem_print_pretty(evs->eigenvecs[i][j], evs->nfs[i], "a");
	printf(" ");
      }
      
      printf("\n");
    }
#endif // DEBUG
  
  for (i = 0; i < 25; i++) {
    get_hecke_ev(ev, genus, evs, ps[i], 0);
    if (!nf_elem_equal_si(ev, test_evs[i], evs->nfs[0])) {
      free_eigenvalues(evs);
      free_hash(genus);
      free_mat(Q);
      return FALSE;
    }

#ifdef DEBUG
    hecke = hecke_matrix(genus, ps[i]);
    print_mat(hecke);

    printf("hecke eigenvalues are:\n");
    for (j = 0; j < evs->num; j++) {
      get_hecke_ev(ev, genus, evs, ps[i], j);
      nf_elem_print_pretty(ev, evs->nfs[j], "a");
      printf(" ");
      nf_elem_clear(ev, evs->nfs[j]);
    }
    printf("\n");
    free_mat(hecke);
#endif // DEBUG
   
  }

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  printf("computing eigenvalues took %f\n", cputime);

  free_eigenvalues(evs);

  free_hash(genus);
  free_mat(Q);
  
  return TRUE;
}
