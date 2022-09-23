#include <time.h>

#include "hecke.h"
#include "arith.h"
#include "genus.h"
#include "hash.h"
#include "matrix_tools.h"
#include "tests.h"

/* Run with test_evs = NULL to just print all the eigenvalues.
   Run with num_evs = 0 to print all the eigenvectors          */

STATUS test(const int* Q_coeffs, int* ps, int* test_evs, int num_evs, int form_idx, const char* inp_type)
{
  int i;
  int j;

  hash_table* slow_genus;
  hash_table* genus;
  clock_t cpuclock_0, cpuclock_1, cpudiff;
  double cputime;
  fmpq_t trace;

  matrix_TYP* Q; //, *hecke;
#ifdef DEBUG_LEVEL_FULL
  matrix_TYP* hecke;
#endif // DEBUG_LEVEL_FULL
  eigenvalues* evs;
  nf_elem_t ev;
  slong deg;

  fmpq_init(trace);
  
  cpuclock_0 = clock();

  Q = init_sym_matrix(Q_coeffs, inp_type);
  slow_genus = get_genus_reps(Q);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  printf("computing genus took %f\n", cputime);

  genus = recalibrate_hash(slow_genus);

  cpuclock_0 = clock();
  cpudiff = cpuclock_0 - cpuclock_1;
  cputime = cpudiff / CLOCKS_PER_SEC;

  printf("recalibrating genus took %f\n", cputime);
  
  // hecke = hecke_matrix(genus, 2);
  // evs = get_eigenvalues(hecke);
  evs = hecke_eigenforms(genus);
  //  free_mat(hecke);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  printf("computing eigenvectors took %f\n", cputime);
  
  // #ifdef DEBUG
  if (num_evs == 0) {
    printf("eigenvectors are:\n");
    for (i = 0; i < evs->num; i++) {
      deg = fmpq_poly_degree(evs->nfs[i]->pol);
      if (deg < 10) {
	  for (j = 0; j < evs->dim; j++) {
	    nf_elem_print_pretty(evs->eigenvecs[i][j], evs->nfs[i], "a");
	    printf(" ");
	  }
      }
      else
	printf("a vector of length %d ", evs->dim);
      printf("over ");
      if (deg < 10) {
	nf_print(evs->nfs[i]);
      }
      else {
	printf("a number field of degree %ld", deg);
      }
      printf("\n");
    }
  }
  // #endif // DEBUG

  printf("traces of hecke eigenvalues are:\n");
  for (i = 0; i < num_evs; i++) {
    get_hecke_ev(ev, genus, evs, ps[i], form_idx);
    nf_elem_trace(trace, ev, evs->nfs[form_idx]);
    fmpq_print(trace);
    // nf_elem_print_pretty(ev, evs->nfs[form_idx], "a");
    printf(" ");
    fflush(stdout); //make sure to print every time it computes an eigenvalue
    nf_elem_clear(ev, evs->nfs[form_idx]);
    if (test_evs != NULL) {
      if (!nf_elem_equal_si(ev, test_evs[i], evs->nfs[form_idx])) {
	free_eigenvalues(evs);
	free_hash(genus);
	free_mat(Q);
	return FAIL;
      }
    }

#ifdef DEBUG_LEVEL_FULL
    hecke = hecke_matrix(genus, ps[i]);
    print_mat(hecke);

    printf("traces of all hecke eigenvalues are:\n");
    for (j = 0; j < evs->num; j++) {
      get_hecke_ev(ev, genus, evs, ps[i], j);
      nf_elem_trace(trace, ev, evs->nfs[j]);
      // nf_elem_print_pretty(ev, evs->nfs[j], "a");
      fmpq_print(trace);
      printf(" ");
      nf_elem_clear(ev, evs->nfs[j]);
    }
    printf("\n");
    free_mat(hecke);
#endif // DEBUG_LEVEL_FULL
   
  }

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  printf("computing eigenvalues took %f\n", cputime);

  fmpq_clear(trace);
  free_eigenvalues(evs);

  free_hash(genus);
  free_mat(Q);
  
  return SUCCESS;
}

STATUS test_61()
{
  int ps[25] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
  int test_evs[25] = {-7, -3, 3, -9, -4, -3, 37, -75, 10, 212, -6, -88, -3, 547, -147, -108, -45,
		      145, -632, -650, 859, -978, 931, -571, 453};
  int Q_coeffs[15] = {2,1,2,0,0,2,0,0,0,4,1,0,0,-1,6};

  return test(Q_coeffs, ps, test_evs, 25, 0, "A");
}

STATUS test_69()
{
  int ps[25] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97};
  int test_evs[25] = {-6, -4, 8, -12, 2, -35, 50, 22, -155, -85, -63, 152, -119, -332, 369, 500, -240, 88, 250, 597, -389, 234, -656, -342, -1342};
  int Q_coeffs[15] = {2,0,2,0,0,2,1,0,0,2,0,0,1,0,12};

  return test(Q_coeffs, ps, test_evs, 25, 0, "A");
}

STATUS test_greedy(const int* Q_coeffs, const int* red_Q_coeffs)
{
  matrix_TYP *Q, *s, *red_Q;
  STATUS ret;

  ret = SUCCESS;
  
  s = init_mat(5, 5, "1");
  Q = init_sym_matrix(Q_coeffs, "A");
  greedy(Q, s, 5, 5);
  red_Q = init_sym_matrix(red_Q_coeffs, "A");
  
#ifdef DEBUG
  print_mat(Q);
#endif // DEBUG

  for (int i = 0; i < 5; i++) {
    for (int j = 0; j < 5; j++) {
      if (Q->array.SZ[i][j] != red_Q->array.SZ[i][j]) {
	ret = FAIL;
	break;
      }
    }
    if (ret == FAIL)
      break;
  }
  
  free_mat(s);
  free_mat(Q);
  free_mat(red_Q);

  return ret;
}

STATUS test_greedy_overflow()
{
  STATUS ret;
  
  int Q1[15] = {8,2,26,2,5,292,3,13,-115,1956,298,69,-60,88,11166};
  int Q2[15] = {12,6,16,-4,-4,842,-5,-5,175,23596,-125,772,-5402,-1517,88494};
  int red_Q[15] = {2,0,2,1,1,2,0,1,0,6,1,1,1,1,8};
  int red_Q2[15] = {2,1,2,1,0,2,1,0,1,6,1,1,1,0,8};

  ret = test_greedy(Q1, red_Q) << 1;

  ret |= test_greedy(Q2, red_Q2);

  return ret;
}

STATUS compute_eigenvectors(const int* Q_coeffs, const char* inp_type)
{
  return test(Q_coeffs, NULL, NULL, 0, 0, inp_type);
}

STATUS compute_eigenvalues(const int* Q_coeffs, int form_idx, int p, const char* inp_type)
{
  return test(Q_coeffs, &p, NULL, 1, form_idx, inp_type);
}

STATUS compute_eigenvalues_up_to(const int* Q_coeffs, int form_idx, int prec, const char* inp_type)
{
  int* ps;
  int num_ps;
  STATUS res;

  num_ps = primes_up_to(&ps, prec);
  
  res = test(Q_coeffs, ps, NULL, num_ps, form_idx, inp_type);

  free(ps);
  
  return res;
}
