#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>

#include "arith.h"
#include "eigenvalues.h"
#include "genus.h"
#include "hecke.h"
#include "io.h"
#include "matrix_tools.h"
#include "tests.h"
#include "weight.h"

void print_conductors(const genus_t genus)
{
  slong c;
  
  printf("The possible conductors are: \n");
  for (c = 0; c < genus->num_conductors; c++)
    printf("%ld ", genus->conductors[c]);
  printf("\n");
  printf("The corresponding dimensions are: ");
  for (c = 0; c < genus->num_conductors; c++)
    printf("%ld ", genus->dims[c]);
  printf("\n");
  return;
}

// Now this also checks initialization from a list of matrices
void compute_genus(genus_t genus, const int* Q_coeffs, const char* format)
{
  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;

  square_matrix_t Q;
  genus_t old_genus;

  cpuclock_0 = clock();

  square_matrix_init_set_symm(Q, Q_coeffs, format);
  genus_init_square_matrix(old_genus, Q, -1);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  fprintf(stderr, "computing genus took %f sec\n", cputime);

  square_matrix_clear(Q);

  cpuclock_0 = clock();
  genus_init_set_square_matrix_vec(genus, old_genus->genus_reps->keys, old_genus->dims[0]);
  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  fprintf(stderr, "recomputing genus took %f sec\n", cputime);

  genus_clear(old_genus);

  print_conductors(genus);
  
  return;
}

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
    if (c == 0)
      if (k == 1)
	get_hecke_ev(ev, genus, evs, ps[i], form_idx);
      else
	get_hecke_ev_nbr_data(ev, genus, evs, ps[i], k, form_idx);
    else
      if (k == 1)
	get_hecke_ev_all_conductors(ev, genus, evs, ps[i], form_idx, c);
      else
	get_hecke_ev_nbr_data_all_conductors(ev, genus, evs, ps[i], k, form_idx, c);
    
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
      if (c == 0)
	if (k == 1)
	  get_hecke_ev(ev, genus, evs, ps[i], j);
	else
	  get_hecke_ev_nbr_data(ev, genus, evs, ps[i], k, j);
      else
	if (k == 1)
	  get_hecke_ev_all_conductors(ev, genus, evs, ps[i], j, c);
	else
	  get_hecke_ev_nbr_data_all_conductors(ev, genus, evs, ps[i], k, j, c);
      
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

STATUS test_genus(genus_t genus, const example_genus_t ex)
{
  slong c;
  
  compute_genus(genus, ex->Q_coeffs, ex->format);

  if (ex->num_conductors != 0)
    assert(genus->num_conductors == ex->num_conductors);
  
  if (ex->dims != NULL)
    for (c = 0; c < ex->num_conductors; c++)
      assert(genus->dims[c] == ex->dims[c]);
  
  return SUCCESS;
}

STATUS test_evs(const genus_t genus, const example_evs_t ex)
{
  int form_idx, k;
  slong c, c2, num_conductors;

  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;

  eigenvalues_t* evs;
  bool has_spinor = false;
  const int* test_evs = NULL;

  num_conductors = genus->num_conductors;
  for (c = 1; c < num_conductors; c++)
    if (genus->dims[c] != 0)
      has_spinor = true;

  cpuclock_0 = clock();
  
  if (has_spinor)
    evs = hecke_eigenforms_all_conductors(genus);
  else
    // !!TODO - change it to be faster when there is no spinor
    evs = hecke_eigenforms_all_conductors(genus);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  fprintf(stderr, "computing eigenvectors took %f sec\n", cputime);

  if (ex->num_forms != NULL)
    for (c = 0; c < num_conductors; c++)
      assert(evs[c]->num == ex->num_forms[c]);
  else
    for (c = 0; c < num_conductors; c++) {
      printf("For conductor %ld, ", genus->conductors[c]);
      print_eigenvectors(evs[c]);
    }

  has_spinor = false;
  for (c = 0; c < num_conductors; c++) {
    eigenvalues_set_lifts(evs[c], 2, c, genus);
    if (c != 0)
      for (form_idx = 0; form_idx < evs[c]->num; form_idx++)
	if (!(evs[c]->is_lift[form_idx]))
	  has_spinor = true;
  }
  
  for (c = 0; c < num_conductors; c++) {
    printf("For conductor %ld:\n", genus->conductors[c]);
    for (form_idx = 0; form_idx < evs[c]->num; form_idx++)
      if (!(evs[c]->is_lift[form_idx])) {
	for (k = 0; k < 2; k++) {
	  if (ex->test_evs != NULL)
	    test_evs = ex->test_evs[c][form_idx][k];
	  cpuclock_0 = clock();
	  if (test_eigenvalues(genus, evs[c], ex->num_ps[k], form_idx,
			       ex->ps[k], test_evs, k+1, c) == FAIL) {
	    for (c2 = 0; c2 < num_conductors; c2++)
	      eigenvalues_clear(evs[c2]);
	    free(evs);
	    return FAIL;
	  }
	  cpuclock_1 = clock();
	  cpudiff = cpuclock_1 - cpuclock_0;
	  cputime = cpudiff / CLOCKS_PER_SEC;
	  // printf("cpudiff = %f, clocks_per_sec = %d\n", cpudiff, CLOCKS_PER_SEC);
	  fprintf(stderr, "computing eigenvalues for k = %d, took %f sec\n", k+1, cputime);
	}
      }
  }

  for (c = 0; c < num_conductors; c++)
    eigenvalues_clear(evs[c]);
  free(evs);
 
  return SUCCESS;
}

void example_genus_init(example_genus_t ex, const int* form, const char* format,
			const int* dims, int num_conductors)
{
  slong c;
  int i;
  
  for (i = 0; i < 15; i++)
    ex->Q_coeffs[i] = form[i];

  ex->format = format;

  ex->num_conductors = num_conductors;

  if (num_conductors == 0) {
    ex->dims = NULL;
  }
  else {
    ex->dims = (int*)malloc(num_conductors*sizeof(int));
    for (c = 0; c < num_conductors; c++)
      ex->dims[c] = dims[c];
  }
  
  return;
}

void example_genus_clear(example_genus_t ex)
{
  if (ex-> dims != NULL) 
    free(ex->dims);
  
  return;
}

void example_evs_init(example_evs_t ex, slong num_conductors,
		      const int* num_forms, const int* num_ps, const int** ps, 
		      const int**** test_evs)
{
  int c,i,j,k;

  if (num_forms != NULL) {
    ex->num_forms = (int*)malloc(num_conductors*sizeof(int));
    for (c = 0; c < num_conductors; c++)
      ex->num_forms[c] = num_forms[c];
  }
  else
    ex->num_forms = NULL;

  if (num_ps != NULL) {
    for (k = 0; k < 2; k++)
      ex->num_ps[k] = num_ps[k];
  }

  for (k = 0; k < 2; k++) {
    if (ex->num_ps[k] == 0)
      ex->ps[k] = NULL;
    else {
      ex->ps[k] = (int*)malloc(ex->num_ps[k] * sizeof(int));
      for (i = 0; i < ex->num_ps[k]; i++)
	ex->ps[k][i] = ps[k][i];
    }
  }

  if (test_evs != NULL) {
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
  else
    ex->test_evs = NULL;
  
  return;
}

void example_evs_clear(example_evs_t ex, int num_conductors)
{
  int c, i, k;

  if (ex->test_evs != NULL) {
    for (c = 0; c < num_conductors; c++) {
      if (ex->test_evs[c] != NULL) {
	for (i = 0; i < ex->num_forms[c]; i++) {
	  for (k = 0; k < 2; k++)
	    free(ex->test_evs[c][i][k]);
	  free(ex->test_evs[c][i]);
	}
	free(ex->test_evs[c]);
      }
    }
  }
  if (num_conductors != 0)
    free(ex->test_evs);

  for (k = 0; k < 2; k++)
    if (ex->ps[k] != NULL)
      free(ex->ps[k]);

  if (ex->num_forms != NULL)
    free(ex->num_forms);
  
  return;
}

void example_genus_init_61(example_genus_t ex)
{
  int form[15] = {2,1,2,0,0,2,0,0,0,4,1,0,0,-1,6};
  char* format = "A";
  int dims[2] = {8,0};

  example_genus_init(ex,form,format,dims,2);
  
  return;
}

void example_evs_init_61(example_evs_t ex)
{
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

  example_evs_init(ex,2,num_forms,num_ps,ps,test_evs);
  
  return;
}

void example_genus_init_69( example_genus_t ex)
{
  int form[15] = {2,0,2,0,0,2,1,0,0,2,0,0,1,0,12};
  char* format = "A";
  int dims[4] = {8,0,0,2};

  example_genus_init(ex,form,format,dims,4);
  
  return;
}

void example_evs_init_69(example_evs_t ex)
{
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

  example_evs_init(ex,4,num_forms,num_ps,ps,test_evs);
  
  return;
}

STATUS test_61()
{
  example_genus_t ex_61_genus;
  example_evs_t ex_61_evs;
  genus_t genus;
  STATUS ret;

  example_genus_init_61(ex_61_genus);
  ret = test_genus(genus, ex_61_genus);
  example_genus_clear(ex_61_genus);

  if (ret == FAIL)
    return ret;
  
  example_evs_init_61(ex_61_evs);
  ret = test_evs(genus, ex_61_evs);
  example_evs_clear(ex_61_evs, genus->num_conductors);

  genus_clear(genus);

  return ret;
}

STATUS test_69()
{
  example_genus_t ex_69_genus;
  example_evs_t ex_69_evs;
  genus_t genus;
  STATUS ret;

  example_genus_init_69(ex_69_genus);
  ret = test_genus(genus, ex_69_genus);
  example_genus_clear(ex_69_genus);

  if (ret == FAIL)
    return ret;
  
  example_evs_init_69(ex_69_evs);
  ret = test_evs(genus, ex_69_evs);
  example_evs_clear(ex_69_evs, genus->num_conductors);

  genus_clear(genus);
  
  return ret;
}

STATUS test_greedy(const int* Q_coeffs, const int* red_Q_coeffs)
{
  square_matrix_t Q, red_Q;
  isometry_t s;

  STATUS ret;

  isometry_init(s);
  square_matrix_init_set_symm(Q, Q_coeffs, "A");
  greedy(Q, s, QF_RANK);
  square_matrix_init_set_symm(red_Q, red_Q_coeffs, "A");
  
#ifdef DEBUG
  square_matrix_print(Q);
#endif // DEBUG

  ret = (square_matrix_is_equal(Q, red_Q) ? SUCCESS : FAIL);
  
  isometry_clear(s);
  square_matrix_clear(Q);
  square_matrix_clear(red_Q);

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

STATUS compute_eigenvectors(const genus_t genus)
{
  example_evs_t dummy;
  STATUS ret;
  int num_ps[2] = {0,0};

  example_evs_init(dummy, genus->num_conductors, NULL, num_ps, NULL, NULL);
  
  ret = test_evs(genus, dummy);
  example_evs_clear(dummy, genus->num_conductors);

  return ret;
}

STATUS compute_eigenvalues(const genus_t genus, int p)
{
  example_evs_t dummy;
  STATUS ret;
  int num_ps[2] = {1,0};
  const int* ps[2] = {
    (int []){p},
    (int []){}
  };

  example_evs_init(dummy, genus->num_conductors, NULL, num_ps, ps, NULL);
  
  ret = test_evs(genus, dummy);
  example_evs_clear(dummy, genus->num_conductors);

  return ret;
}

STATUS compute_eigenvalues_up_to(const genus_t genus, int prec)
{
  example_evs_t dummy;
  STATUS ret;
  int num_ps[2];
  int* ps[2];

  /*
  num_ps[0] = primes_up_to(&(ps[0]), prec);
  num_ps[1] = primes_up_to(&(ps[1]), floor(sqrt(prec)));
  */
  
  num_ps[0] = primes_up_to_prime_to(&(ps[0]), prec, fmpz_get_si(genus->disc));
  num_ps[1] = primes_up_to_prime_to(&(ps[1]), floor(sqrt(prec)), fmpz_get_si(genus->disc));
  
  example_evs_init(dummy, genus->num_conductors, NULL, num_ps, (const int**)ps, NULL);
  ret = test_evs(genus, dummy);
  example_evs_clear(dummy, genus->num_conductors);

  free(ps[0]);
  free(ps[1]);
  
  return ret;
}

STATUS compute_hecke_col_all_conds(const genus_t genus, slong p, int gen_idx)
{
  int** hecke;
  slong c, i;

  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;  

  cpuclock_0 = clock();
  
  hecke = hecke_col_all_conds_sparse(p, gen_idx, genus);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;
  
  // printf("{%ld : {", p);
  printf("{");
  for (c = 0; c < genus->num_conductors; c++) {
    printf("%ld : ", genus->conductors[c]);
    printf("[");
    for (i = 0; i < genus->dims[c]; i++) {
      printf("%d", hecke[c][i]);
      if (i != genus->dims[c]-1)
	printf(",");
    }
    printf("]");
    if (c != genus->num_conductors-1)
      printf(",");
  }
  printf("}");
  //printf("} }");

  fprintf(stderr, "computing hecke took %f sec\n", cputime);
  
  for (c = 0; c < genus->num_conductors; c++)
    free(hecke[c]);

  free(hecke);

  return SUCCESS;
}

STATUS compute_hecke_col(const genus_t genus, slong p, slong c)
{
  int* hecke;
  int c_idx, i, gen_idx;

  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;  

  for (c_idx = 0; (c_idx < genus->num_conductors) && (genus->conductors[c_idx] != c);) c_idx++;

  gen_idx = 0;
  while (genus->lut_positions[c_idx][gen_idx] == -1)
    gen_idx++;
  
  // right now we only handle differently the trivial spinor norm case
  if (c_idx != 0)
    return compute_hecke_col_all_conds(genus, p, gen_idx);

  hecke = (int*)malloc(genus->dims[0] * sizeof(int));

  cpuclock_0 = clock();

  hecke_col(hecke, p, gen_idx, genus);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;

  // We stop printing out the p, p will be part of the filename (externally)
  //   printf("{%ld : {", p);
  
  printf("{%ld : ", genus->conductors[c]);
  printf("[");
  for (i = 0; i < genus->dims[c]; i++) {
    printf("%d", hecke[i]);
    if (i != genus->dims[c]-1)
      printf(",");
  }
  printf("]");

  printf("}");

  fprintf(stderr, "computing hecke took %f sec\n", cputime);
  
  free(hecke);

  return SUCCESS;
}

STATUS compute_hecke_matrix(const genus_t genus, slong p, slong c)
{
  matrix_TYP* hecke;
  int c_idx, gen_idx;

  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;  

  for (c_idx = 0; (c_idx < genus->num_conductors) && (genus->conductors[c_idx] != c);) c_idx++;

  gen_idx = 0;
  while (genus->lut_positions[c_idx][gen_idx] == -1)
    gen_idx++;
  
  // right now we only handle differently the trivial spinor norm case
  if (c_idx != 0)
    return compute_hecke_matrix_all_conds(genus, p);

  cpuclock_0 = clock();

  hecke = hecke_matrix(genus, p);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;

  //printf("{%ld : {", p);
  printf("{%ld : ", genus->conductors[0]);
  print_mat_dense(hecke);
  printf("}");
  //  printf("} }");

  fprintf(stderr, "computing hecke took %f sec\n", cputime);
  
  free_mat(hecke);

  return SUCCESS;
}

STATUS compute_hecke_matrix_all_conds(const genus_t genus, slong p)
{
  matrix_TYP** hecke;
  slong c;
 
  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;
  
  cpuclock_0 = clock();
 
  hecke = hecke_matrices_all_conductors(genus,p);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;

  // printf("{%ld : {", p);
  printf("{");
  for (c = 0; c < genus->num_conductors; c++) {
    printf("%ld : ", genus->conductors[c]);
    print_mat_dense(hecke[c]);
    if (c != genus->num_conductors-1)
      printf(",");
  }
  printf("}");
  //  printf("} }");

  fprintf(stderr, "computing hecke matrices took %f sec\n", cputime);
  
  for (c = 0; c < genus->num_conductors; c++)
    free_mat(hecke[c]);

  free(hecke);
 
  return SUCCESS;
  
}

STATUS compute_first_hecke_matrix_all_conds(const genus_t genus)
{
  slong p;
  fmpz_t prime;

  fmpz_init(prime);
  fmpz_set_ui(prime, 1);
  
  do {
    fmpz_nextprime(prime, prime, true);
  }
  while (fmpz_divisible(genus->disc,prime));

  p = fmpz_get_si(prime);

  compute_hecke_matrix_all_conds(genus, p);
 
  fmpz_clear(prime);
 
  return SUCCESS;
}

STATUS compute_first_large_hecke(const genus_t genus)
{
  square_matrix_t** hecke;
  slong i,j,p;
  fmpz_t prime;

  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;
  
  cpuclock_0 = clock();

  fmpz_init(prime);
  fmpz_set_ui(prime, 1);
  
  do {
    fmpz_nextprime(prime, prime, true);
  }
  while (fmpz_divisible(genus->disc,prime));

  p = fmpz_get_si(prime);
  hecke = hecke_matrices_isometries(genus,p);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;

  printf("%ld, ", p);
  printf("[");
  for (i = 0; i < genus->dims[0]; i++) {
    printf("[");
    for (j = 0; j < genus->dims[0]; j++) {
      square_matrix_print(hecke[i][j]);
      if (j != genus->dims[0]-1)
	printf(",");
    }
    printf("]");
    if (i != genus->dims[0]-1)
      printf(",");
  }
  printf("]");

  fprintf(stderr, "computing hecke matrices took %f sec\n", cputime);
  
  for (i = 0; i < genus->dims[0]; i++) {
    for (j = 0; j < genus->dims[0]; j++) {
      square_matrix_clear(hecke[i][j]);
    }
    free(hecke[i]);
  }

  free(hecke);
  fmpz_clear(prime);
 
  return SUCCESS;
}
