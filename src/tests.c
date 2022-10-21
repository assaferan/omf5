#include <assert.h>
#include <time.h>

#include "eigenvalues.h"
#include "hecke.h"
#include "arith.h"
#include "genus.h"
#include "matrix_tools.h"
#include "tests.h"

void print_eigenvectors(const eigenvalues_t evs)
{
  int i, j;
  slong deg;
  
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
      printf("a vector of length %ld ", evs->dim);
    printf("over ");
    if (deg < 10) {
      nf_print(evs->nfs[i]);
    }
    else {
      printf("a number field of degree %ld", deg);
    }
    printf("\n");
  }

  return;
}

STATUS test_eigenvalues(const genus_t genus, const eigenvalues_t evs,
			int num_evs, int form_idx, const int* ps,
			const int* test_evs, int k, int c)
{
  int i;
  nf_elem_t ev;
  fmpq_t trace;
#ifdef DEBUG_LEVEL_FULL
  int j;
  matrix_TYP* hecke;
#endif // DEBUG_LEVEL_FULL

  fmpq_init(trace);
  
  printf("traces of hecke eigenvalues of T_{p^%d} are:\n", k);
  for (i = 0; i < num_evs; i++) {
    if (k == 2)
      get_hecke_ev_nbr_data_all_conductors(ev, genus, evs, ps[i], k, form_idx, c);
    else {
#ifdef NBR_DATA
      get_hecke_ev_nbr_data_all_conductors(ev, genus, evs, ps[i], k, form_idx, c);
#else
      get_hecke_ev(ev, genus, evs, ps[i], form_idx);
#endif
    }
    nf_elem_trace(trace, ev, evs->nfs[form_idx]);
    fmpq_print(trace);
    // nf_elem_print_pretty(ev, evs->nfs[form_idx], "a");
    printf(" ");
    fflush(stdout); //make sure to print every time it computes an eigenvalue
    nf_elem_clear(ev, evs->nfs[form_idx]);
    if (test_evs != NULL) {
      // if (!nf_elem_equal_si(ev, test_evs[i], evs->nfs[form_idx])) {
      if (!fmpq_equal_si(trace, test_evs[i])) {
	fmpq_clear(trace);
	return FAIL;
      }
    }

#ifdef DEBUG_LEVEL_FULL
    if (k == 1) {
      hecke = hecke_matrix(genus, ps[i]);
      print_mat(hecke);
    }
    printf("traces of all hecke eigenvalues are:\n");
    for (j = 0; j < evs->num; j++) {
      if (k == 1) {
#ifdef NBR_DATA
	get_hecke_ev_nbr_data(ev, genus, evs, ps[i], k, j);
#else
	get_hecke_ev(ev, genus, evs, ps[i], j);
#endif // NBR_DATA
      }
      else
	get_hecke_ev_nbr_data(ev, genus, evs, ps[i], k, j);
      nf_elem_trace(trace, ev, evs->nfs[j]);
      // nf_elem_print_pretty(ev, evs->nfs[j], "a");
      fmpq_print(trace);
      printf(" ");
      nf_elem_clear(ev, evs->nfs[j]);
    }
    printf("\n");
    if (k == 1)
      free_mat(hecke);
#endif // DEBUG_LEVEL_FULL
   
  }

  printf("\n");

  fmpq_clear(trace);
  return SUCCESS;
}

/* Run with test_evs = NULL to just print all the eigenvalues.
   Run with num_evs = 0 to print all the eigenvectors          */

STATUS test(const example_t ex)
{
  int form_idx, k;
  slong c, c2;

  genus_t genus;
  clock_t cpuclock_0, cpuclock_1, cpudiff;
  double cputime;

  matrix_TYP* Q;
  eigenvalues_t* evs;
  const int* test_evs = NULL;
  
  cpuclock_0 = clock();

  Q = init_sym_matrix(ex->Q_coeffs, ex->format);
  genus_init(genus, Q);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  printf("computing genus took %f\n", cputime);

  if (ex->num_conductors != 0)
    assert(genus->num_conductors == ex->num_conductors);
  
  if (ex->dims != NULL)
    for (c = 0; c < ex->num_conductors; c++)
      assert(genus->dims[c] == ex->dims[c]);

  evs = hecke_eigenforms_all_conductors(genus);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  printf("computing eigenvectors took %f\n", cputime);

  if (ex->num_forms != NULL)
    for (c = 0; c < ex->num_conductors; c++)
      assert(evs[c]->num == ex->num_forms[c]);
  else
    for (c = 0; c < genus->num_conductors; c++)
      print_eigenvectors(evs[c]);
  
  for (c = 0; c < genus->num_conductors; c++) {
    eigenvalues_set_lifts(evs[c], 2, c, genus);
    for (form_idx = 0; form_idx < evs[c]->num; form_idx++)
      if (!(evs[c]->is_lift[form_idx])) {
	for (k = 0; k < 2; k++) {
	  if (ex->num_conductors != 0)
	    test_evs = ex->test_evs[c][form_idx][k];
	  if (test_eigenvalues(genus, evs[c], ex->num_ps[k], form_idx,
			       ex->ps[k], test_evs, k+1, c) == FAIL) {
	    for (c2 = 0; c2 < genus->num_conductors; c2++)
	      eigenvalues_clear(evs[c2]);
	    free(evs);
	    genus_clear(genus);
	    free_mat(Q);
	    return FAIL;
	  }
	}
      }
  }
  
  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  printf("computing eigenvalues took %f\n", cputime);

  for (c = 0; c < genus->num_conductors; c++)
    eigenvalues_clear(evs[c]);
  free(evs);

  genus_clear(genus);
  free_mat(Q);
  
  return SUCCESS;
}

void example_init(example_t ex, const int* form, const char* format, const int* dims,
		  const int* num_forms, const int* num_ps, const int** ps, 
		  const int**** test_evs, int num_conductors)
{
  int c,i,j,k;
  
  for (i = 0; i < 15; i++)
    ex->Q_coeffs[i] = form[i];

  ex->format = format;

  ex->num_conductors = num_conductors;

  if (num_conductors == 0) {
    ex->dims = NULL;
    ex->num_forms = NULL;
  } else {
    ex->dims = (int*)malloc(num_conductors*sizeof(int));
    for (c = 0; c < num_conductors; c++)
      ex->dims[c] = dims[c];
    
    ex->num_forms = (int*)malloc(num_conductors*sizeof(int));
    for (c = 0; c < num_conductors; c++)
      ex->num_forms[c] = num_forms[c];
  }

  for (k = 0; k < 2; k++)
    ex->num_ps[k] = num_ps[k];

  for (k = 0; k < 2; k++) {
    if (ex->num_ps[k] == 0)
      ex->ps[k] = NULL;
    else {
      ex->ps[k] = (int*)malloc(ex->num_ps[k] * sizeof(int));
      for (i = 0; i < ex->num_ps[k]; i++)
	ex->ps[k][i] = ps[k][i];
    }
  }

  if (num_conductors != 0) {
    ex->test_evs = (int****)malloc(num_conductors * sizeof(int***));
    for (c = 0; c < num_conductors; c++) {
      if (ex->num_forms[c] == 0)
	ex->test_evs[c] = NULL;
      else {
	ex->test_evs[c] = (int***)malloc(ex->num_forms[c] * sizeof(int**));
	for (i = 0; i < ex->num_forms[c]; i++) {
	  ex->test_evs[c][i] = (int**)malloc(2 * sizeof(int**));
	  for (k = 0; k < 2; k++) {
	    ex->test_evs[c][i][k] = (int*)malloc(ex->num_ps[k] * sizeof(int));
	    for (j = 0; j < ex->num_ps[k]; j++)
	      ex->test_evs[c][i][k][j] = test_evs[c][i][k][j];
	  }
	}
      }
    }
  }
  
  return;
}

void example_clear(example_t ex)
{
  int c, i, k;
  
  for (c = 0; c < ex->num_conductors; c++) {
    if (ex->test_evs[c] != NULL) {
      for (i = 0; i < ex->num_forms[c]; i++) {
	for (k = 0; k < 2; k++)
	  free(ex->test_evs[c][i][k]);
	free(ex->test_evs[c][i]);
      }
      free(ex->test_evs[c]);
    }
  }
  if (ex->num_conductors != 0)
    free(ex->test_evs);

  for (k = 0; k < 2; k++)
    if (ex->ps[k] != NULL)
      free(ex->ps[k]);

  if (ex->num_forms != NULL)
    free(ex->num_forms);
  
  if (ex-> dims != NULL) 
    free(ex->dims);

  return;
}

void example_init_61(example_t ex)
{
  int form[15] = {2,1,2,0,0,2,0,0,0,4,1,0,0,-1,6};
  char* format = "A";
  int dims[2] = {8,0};
  int num_forms[2] = {3,0};
  int num_ps[2] = {25,4};
  const int* ps[2] = {
    (int[]){2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97},
    (int []){2, 3, 5, 7}
  };

  const int*** test_evs[2] = {
    (const int** []) {
	
      (const int* []) {
	(int[]){-7, -3, 3, -9, -4, -3, 37, -75, 10, 212, -6, -88, -3, 547, -147, -108, -45, 145, -632, -650, 859, -978, 931, -571, 453},
	(int[]){7,-9,-9,-42}
      },
	
      (const int* []) {
	(int []){15, 40, 156, 400, 1464, 2380, 5220, 7240, 12720, 25260, 30784, 52060, 70644, 81400, 106080, 151740, 208920, 227043, 305320, 363024, 394420, 499360, 578760, 712980, 922180},
	(int []){30, 120, 780, 2800} 
      },
	
      (const int* []) {
	(int []){29,59,164,309,612,1038,1786,2287,3235,4757,5717,8337,9491,11542,13464,16843,19388,23058,27780,29755,34353,38362,40541,47818,61861},
	(int []){-3,-4,48,72}
      }
    },
      
    (const int** []) {
    }
    
  };

  example_init(ex,form,format,dims,num_forms,num_ps,ps,test_evs,2);
  
  return;
}

void example_init_69(example_t ex)
{
  int form[15] = {2,0,2,0,0,2,1,0,0,2,0,0,1,0,12};
  char* format = "A";
  int dims[4] = {8,0,0,2};
  int num_forms[4] = {4,0,0,1};
  int num_ps[2] = {25,4};
  const int* ps[2] = {
    (int[]){2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97},
    (int []){2, 3, 5, 7}
  };

  const int*** test_evs[] = {
    (const int** []) {
	
      (const int* []) {
	(int[]){-6, -4, 8, -12, 2, -35, 50, 22, -155, -85, -63, 152, -119, -332, 369, 500, -240, 88, 250, 597, -389, 234, -656, -342, -1342},
	(int[]){5, 1, -32, -40}
      },
	
      (const int* []) {
	(int []){15,31,156, 400, 1464, 2380, 5220, 7240, 13249, 25260, 30784, 52060, 70644, 81400, 106080, 151740, 208920, 230764, 305320, 363024, 394420, 499360, 578760, 712980, 922180},
	(int []){30, 1, 780, 2800} 
      },
	
      (const int* []) {
	(int []){8, 30, 34, 102, 204, 340, 462, 714, 1058, 1524, 2308, 2952, 3080, 3658, 4632, 6214, 7048, 6656, 8238, 10712, 11456, 11442, 13844, 16122, 20636},
	(int []){-6, 2, -108, 16}
      },
	
      (const int* []) {
	(int []){26, 19, 134, 240, 536, 839, 1322, 1616, 2116, 3501, 3775, 5794, 6763, 7570, 8347, 11218, 13020, 15882, 18712, 20047, 23117, 24442, 28030, 34382, 39086},
	(int []){18, 4, 180, 320}
      }
    },
      
    (const int** []) {},
    (const int** []) {},
    (const int** []) {
      (const int* []) {
	(int []){10, -20, 68, 84, 192, 304, 708, 612, -1152, 1536, 1816, 2928, 3448, 3364, 4464, 5756, 7040, 6800, 8124, 10000, 11624, 14412, 15424, 16764, 18752},
	(int []){0, 0, 96, -128}
      }
    }
    
  };

  example_init(ex,form,format,dims,num_forms,num_ps,ps,test_evs,4);
  
  return;
}

STATUS test_61()
{
  example_t ex_61;
  STATUS ret;
  
  example_init_61(ex_61);

  ret = test(ex_61);

  example_clear(ex_61);

  return ret;
}

STATUS test_69()
{
  example_t ex_69;
  STATUS ret;
  
  example_init_69(ex_69);

  ret = test(ex_69);

  example_clear(ex_69);

  return ret;
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
  example_t dummy;
  STATUS ret;
  int num_ps[2] = {0,0};
  
  example_init(dummy, Q_coeffs, inp_type, NULL, NULL, num_ps, NULL, NULL, 0);
  ret = test(dummy);
  example_clear(dummy);

  return ret;
}

STATUS compute_eigenvalues(const int* Q_coeffs, int p, const char* inp_type)
{
  example_t dummy;
  STATUS ret;
  int num_ps[2] = {1,0};
  const int* ps[2] = {
    (int []){p},
    (int []){}
  };
  
  example_init(dummy, Q_coeffs, inp_type, NULL, NULL, num_ps, ps, NULL, 0);
  ret = test(dummy);
  example_clear(dummy);

  return ret;
}

// !! TODO - for now we leave the form_idx as a dummy argument, fix that later
// together with argument parsing

STATUS compute_eigenvalues_up_to(const int* Q_coeffs, int form_idx,
				 int prec, const char* inp_type)
{
  example_t dummy;
  STATUS ret;
  int num_ps[2];
  int* ps[2];

  num_ps[0] = primes_up_to(&(ps[0]), prec);
  num_ps[1] = primes_up_to(&(ps[1]), floor(sqrt(prec)));

  example_init(dummy, Q_coeffs, inp_type, NULL, NULL, num_ps, (const int**)ps, NULL, 0);
  ret = test(dummy);
  example_clear(dummy);

  free(ps[0]);
  free(ps[1]);
  
  return ret;
}
