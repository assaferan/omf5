#include <time.h>

#include "genus.h"
#include "hash.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "tests.h"

/* Run with test_evs = NULL to just print all the eigenvalues.
   Run with num_evs = 0 to print all the eigenvectors          */

int test(int* Q_coeffs, int* ps, int* test_evs, int num_evs, int form_idx)
{
  int i;
  int j;

  hash_table* genus;
  clock_t cpuclock_0, cpuclock_1, cpudiff;
  double cputime;

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
  
  // #ifdef DEBUG
  if (num_evs == 0) {
    printf("eigenvectors are:\n");
    for (i = 0; i < evs->num; i++) {
      for (j = 0; j < evs->dim; j++) {
	nf_elem_print_pretty(evs->eigenvecs[i][j], evs->nfs[i], "a");
	printf(" ");
      }
      
      printf("\n");
    }
  }
  // #endif // DEBUG

  printf("hecke eigenvalues are:\n");
  for (i = 0; i < num_evs; i++) {
    get_hecke_ev(ev, genus, evs, ps[i], form_idx);
    nf_elem_print_pretty(ev, evs->nfs[form_idx], "a");
    printf(" ");
    fflush(stdout); //make sure to print every time it computes an eigenvalue
    nf_elem_clear(ev, evs->nfs[form_idx]);
    if (test_evs != NULL) {
      if (!nf_elem_equal_si(ev, test_evs[i], evs->nfs[form_idx])) {
	free_eigenvalues(evs);
	free_hash(genus);
	free_mat(Q);
	return FALSE;
      }
    }

#ifdef DEBUG
    hecke = hecke_matrix(genus, ps[i]);
    print_mat(hecke);

    printf("all hecke eigenvalues are:\n");
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

int test_61()
{
  int ps[25] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
  int test_evs[25] = {-7, -3, 3, -9, -4, -3, 37, -75, 10, 212, -6, -88, -3, 547, -147, -108, -45,
		      145, -632, -650, 859, -978, 931, -571, 453};
  int Q_coeffs[15] = {2,1,2,0,0,2,0,0,0,4,1,0,0,-1,6};

  return test(Q_coeffs, ps, test_evs, 25, 0);
}

int test_69()
{
  int ps[25] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
  int test_evs[25] = {-6, -5, 8, -12, 2, -35, 50, 22, -156, -85, -63, 152, -119, -332, 369, 500, -240, 88, 250, 597, -389, 234, -656, -342, -1342};
  int Q_coeffs[15] = {2,0,2,0,0,2,1,0,0,2,0,0,1,0,12};

  return test(Q_coeffs, ps, test_evs, 25, 0);
}

int compute_eigenvectors(int* Q_coeffs)
{
  return test(Q_coeffs, NULL, NULL, 0, 0);
}

int compute_eigenvalues(int* Q_coeffs, int form_idx, int p)
{
  return test(Q_coeffs, &p, NULL, 1, form_idx);
}

int compute_eigenvalues_up_to(int* Q_coeffs, int form_idx, int prec)
{
  int* ps;
  int* sieve;
  int num_ps, p, i;

  num_ps = 0;
  sieve = malloc((prec+2) * sizeof(int));
  ps = malloc(prec * sizeof(int));
  // since prec is small we construct the first primes using aristothenes sieve
  // replace by fmpz_nextprime if improving.
  for (i = 0; i <= prec+1; i++)
    sieve[i] = i;
  num_ps = 0;
  p = 2;
  while(p <= prec) {
    ps[num_ps] = p;
    num_ps++;
    for (i = p; i <= prec; i+= p)
      sieve[i] = -1;
    while(sieve[p] == -1) p++;
  }

  free(sieve);
  
  int res = test(Q_coeffs, ps, NULL, num_ps, form_idx);

  free(ps);
  
  return res;
}

