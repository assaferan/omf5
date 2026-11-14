#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "arith.h"
#include "eigenvalues.h"
#include "genus.h"
#include "hecke.h"
#include "io.h"
#include "mass.h"
#include "matrix_tools.h"
#include "tests.h"
#include "weight.h"

void print_conductors(const genus_t genus)
{
  slong c;
  
  fprintf(stderr, "The possible conductors are: \n");
  for (c = 0; c < genus->num_conductors; c++)
    fprintf(stderr, "%ld ", genus->conductors[c]);
  fprintf(stderr, "\n");
  fprintf(stderr, "The corresponding dimensions are: ");
  for (c = 0; c < genus->num_conductors; c++)
    fprintf(stderr, "%ld ", genus->dims[c]);
  fprintf(stderr, "\n");
  return;
}

// Now this also checks initialization from a list of matrices
void compute_genus(genus_t genus, const Z64* Q_coeffs, const char* format)
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
  char* nf_str;
  
  fprintf(stderr, "eigenvectors are:\n");
  for (i = 0; i < evs->num; i++) {
    /*
    if (evs->lift_type[i] == O)
      continue;
    */
    deg = fmpq_poly_degree(evs->nfs[i]->pol);
    if (deg < 10) {
      for (j = 0; j < evs->dim; j++) {
	nf_str = nf_elem_get_str_pretty(evs->eigenvecs[i][j], "a", evs->nfs[i]);
	fprintf(stderr, "%s ", nf_str);
	free(nf_str);
      }
      fprintf(stderr, "over Number field with defining polynomial ");
      fmpq_poly_fprint_pretty(stderr, evs->nfs[i]->pol, "x");
    }
    else
      fprintf(stderr, "a vector of length %ld over a number field of degree %ld",
	      evs->dim, deg);
    fprintf(stderr, "\n");
  }

  return;
}

STATUS test_eigenvalues(const genus_t genus, const eigenvalues_t evs,
			int num_evs, int form_idx, const int* ps,
			const int* test_evs, int k, int c, bool ev_print)
{
  int i;
  nf_elem_t ev;
  fmpq_t trace;
#ifdef DEBUG_LEVEL_FULL
  int j;
  matrix_TYP* hecke;
#endif // DEBUG_LEVEL_FULL

  fmpq_init(trace);

  fprintf(stderr, "traces of hecke eigenvalues of T_{p^%d} are:\n", k);

  if (ev_print) {
    if (k == 1)
      printf("'lambda_p' : [");
    if (k == 2)
      printf("'lambda_p_square' : [");
  }
  
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

    if (ev_print) {
      printf("'");
      nf_elem_print_pretty(ev, evs->nfs[form_idx], "a");
      printf("',");
      fflush(stdout);
    }
    
    nf_elem_trace(trace, ev, evs->nfs[form_idx]);
    fmpq_fprint(stderr, trace);
    fprintf(stderr, " ");
    fflush(stderr); // make sure to print every time it computes an eigenvalue
    
    nf_elem_clear(ev, evs->nfs[form_idx]);
    if (test_evs != NULL) {
      // if (!nf_elem_equal_si(ev, test_evs[i], evs->nfs[form_idx])) {
      if (!fmpq_equal_si(trace, test_evs[i])) {
	fprintf(stderr, "\nFAIL: trace of T_{%d^%d} for form %d is ", ps[i], k, form_idx);
	fmpq_fprint(stderr, trace);
	fprintf(stderr, ", expected %d\n", test_evs[i]);
	fmpq_clear(trace);
	return FAIL;
      }
    }

#ifdef DEBUG_LEVEL_FULL
    if (k == 1) {
      hecke = hecke_matrix(genus, ps[i]);
      print_mat(hecke);
    }
    fprintf(stderr, "traces of all hecke eigenvalues are:\n");
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
      fmpq_fprint(stderr, trace);
      fprintf(stderr, " ");
      nf_elem_clear(ev, evs->nfs[j]);
    }
    fprintf(stderr, "\n");
    if (k == 1)
      free_mat(hecke);
#endif // DEBUG_LEVEL_FULL
   
  }

  if (ev_print)
     printf("], ");

  fprintf(stderr, "\n");

  fmpq_clear(trace);
  return SUCCESS;
}

/* Run with test_evs = NULL to just print all the eigenvalues.
   Run with num_evs = 0 to print all the eigenvectors          */

STATUS test_genus(genus_t genus, const example_genus_t ex)
{
  slong c;
  
  compute_genus(genus, ex->Q_coeffs, ex->format);

  if ((ex->num_conductors != 0) && (genus->num_conductors != ex->num_conductors)) {
    fprintf(stderr, "FAIL: the genus has %ld conductors, expected %d\n",
	    genus->num_conductors, ex->num_conductors);
    return FAIL;
  }
  
  if (ex->dims != NULL)
    for (c = 0; c < ex->num_conductors; c++)
      if (genus->dims[c] != ex->dims[c]) {
	fprintf(stderr, "FAIL: the space of conductor %ld has dimension %ld, expected %d\n",
		genus->conductors[c], (slong)genus->dims[c], ex->dims[c]);
	return FAIL;
      }
  
  return SUCCESS;
}

STATUS test_evs(const genus_t genus, const example_evs_t ex, bool nonlifts, int num_idxs, const int* idxs)
{
  int form_idx, k, i;
  slong c, c2, num_conductors;
  int total_idx;
  slong num_lift_ps;

  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;

  eigenvalues_t* evs;
  bool has_spinor = false;
  const int* test_evs = NULL;

  char lift_types[5] = {'G', 'Y', 'P', 'F', 'O'};

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

  if (ex->num_forms != NULL) {
    for (c = 0; c < num_conductors; c++)
      if (evs[c]->num != ex->num_forms[c]) {
	fprintf(stderr, "FAIL: conductor %ld has %ld eigenforms, expected %d\n",
		genus->conductors[c], (slong)evs[c]->num, ex->num_forms[c]);
	for (c2 = 0; c2 < num_conductors; c2++)
	  eigenvalues_clear(evs[c2]);
	free(evs);
	return FAIL;
      }
  }
  else
    for (c = 0; c < num_conductors; c++) {
      fprintf(stderr, "For conductor %ld, ", genus->conductors[c]);
      print_eigenvectors(evs[c]);
    }

  // The lifts are detected at the primes p not dividing the discriminant
  // for which T_{p^2} is computed, and at two primes at least.
  // These are the first such primes, as eigenvalues_set_lifts expects.
  num_lift_ps = 0;
  for (i = 0; i < ex->num_ps[1]; i++)
    if (fmpz_fdiv_ui(genus->disc, ex->ps[1][i]) != 0)
      num_lift_ps++;
  if (num_lift_ps < 2)
    num_lift_ps = 2;

  // When indices are given, only the forms at these indices are computed.
  // They are the ones which are given the type G below.
  if (num_idxs > 0)
    nonlifts = true;

  has_spinor = false;
  total_idx = 0;
  for (c = 0; c < num_conductors; c++) {
    if (num_idxs == 0)
      // No indices were given, so detect the lifts.
      eigenvalues_set_lifts(evs[c], num_lift_ps, c, genus);
    else
      // The forms at the given indices are taken to be the nonlifts,
      // and all the others are skipped. The indices can be in any order.
      for (form_idx = 0; form_idx < evs[c]->num; form_idx++, total_idx++) {
	evs[c]->lift_type[form_idx] = O;
	for (i = 0; i < num_idxs; i++)
	  if (idxs[i] == total_idx)
	    evs[c]->lift_type[form_idx] = G;
      }
    if (c != 0)
      for (form_idx = 0; form_idx < evs[c]->num; form_idx++)
	if ((!nonlifts) || (evs[c]->lift_type[form_idx] == G))
	  has_spinor = true;
  }

  // Here total_idx is the number of forms, unless no indices were given.
  for (i = 0; i < num_idxs; i++)
    if ((idxs[i] < 0) || (idxs[i] >= total_idx))
      fprintf(stderr, "Warning: there is no form at index %d, the indices of the %d forms start from 0\n",
	      idxs[i], total_idx);

  // We report the forms whose eigenvalues are going to be computed
  total_idx = 0;
  for (c = 0; c < num_conductors; c++)
    for (form_idx = 0; form_idx < evs[c]->num; form_idx++, total_idx++)
      if ((!nonlifts) || (evs[c]->lift_type[form_idx] == G)) {
	if (evs[c]->lift_type[form_idx] == G)
	  fprintf(stderr, "Found nonlift");
	else if (evs[c]->lift_type[form_idx] == O)
	  // an old form
	  fprintf(stderr, "Found form of type O");
	else
	  fprintf(stderr, "Found lift of type %c", lift_types[evs[c]->lift_type[form_idx]]);
	fprintf(stderr, " at total_idx = %d", total_idx);
	fprintf(stderr, ", c = %ld, form_idx = %d\n", c, form_idx);
      }
    
  printf("{");
  for (c = 0; c < num_conductors; c++) {
    printf("%ld : [", genus->conductors[c]);
    fprintf(stderr, "For conductor %ld:\n", genus->conductors[c]);
    for (form_idx = 0; form_idx < evs[c]->num; form_idx++)
      if ((!nonlifts) || (evs[c]->lift_type[form_idx] == G)) {
	printf("{ 'aut_rep_type' : '%c', ", lift_types[evs[c]->lift_type[form_idx]]); 
	printf("'field_poly' : '");
	fmpq_poly_print(evs[c]->nfs[form_idx]->pol);
	printf("', ");
	for (k = 0; k < 2; k++) {
	  if (ex->test_evs != NULL)
	    test_evs = ex->test_evs[c][form_idx][k];
	  cpuclock_0 = clock();
	  if (test_eigenvalues(genus, evs[c], ex->num_ps[k], form_idx,
			       ex->ps[k], test_evs, k+1, c, true) == FAIL) {
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
	printf("}, ");
      }
    printf("], ");
  }
  printf("}");

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
	(int[]){5, 0, -32, -40}
      },
	
      (const int* []) {
	(int []){15,31,156, 400, 1464, 2380, 5220, 7240, 13249, 25260, 30784, 52060, 70644, 81400, 106080, 151740, 208920, 230764, 305320, 363024, 394420, 499360, 578760, 712980, 922180},
	(int []){30, 0, 780, 2800} 
      },
	
      (const int* []) {
	(int []){8, 30, 34, 102, 204, 340, 462, 714, 1058, 1524, 2308, 2952, 3080, 3658, 4632, 6214, 7048, 6656, 8238, 10712, 11456, 11442, 13844, 16122, 20636},
	(int []){-6, 0, -108, 16}
      },
	
      (const int* []) {
	(int []){26, 19, 134, 240, 536, 839, 1322, 1616, 2116, 3501, 3775, 5794, 6763, 7570, 8347, 11218, 13020, 15882, 18712, 20047, 23117, 24442, 28030, 34382, 39086},
	(int []){18, 0, 180, 320}
      }
    },
      
    (const int** []) {},
    (const int** []) {},
    (const int** []) {
      (const int* []) {
	(int []){10, 18, 68, 84, 192, 304, 708, 612, 1150, 1536, 1816, 2928, 3448, 3364, 4464, 5756, 7040, 6800, 8124, 10000, 11624, 14412, 15424, 16764, 18752},
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

  if (ret == FAIL) {
    genus_clear(genus);
    return ret;
  }
  
  example_evs_init_61(ex_61_evs);
  // Check all the forms, lifts included.
  ret = test_evs(genus, ex_61_evs, false, 0, NULL);
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

  if (ret == FAIL) {
    genus_clear(genus);
    return ret;
  }
  
  example_evs_init_69(ex_69_evs);
  // Check all the forms, lifts included.
  ret = test_evs(genus, ex_69_evs, false, 0, NULL);
  example_evs_clear(ex_69_evs, genus->num_conductors);

  genus_clear(genus);
  
  return ret;
}

STATUS test_greedy(const Z64* Q_coeffs, const Z64* red_Q_coeffs)
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
  
  Z64 Q1[15] = {8,2,26,2,5,292,3,13,-115,1956,298,69,-60,88,11166};
  Z64 Q2[15] = {12,6,16,-4,-4,842,-5,-5,175,23596,-125,772,-5402,-1517,88494};
  Z64 red_Q[15] = {2,0,2,1,1,2,0,1,0,6,1,1,1,1,8};
  Z64 red_Q2[15] = {2,1,2,1,0,2,1,0,1,6,1,1,1,0,8};

  ret = test_greedy(Q1, red_Q) << 1;

  ret |= test_greedy(Q2, red_Q2);

  return ret;
}

// Tests whether x is found to be a square in Q[a] / (f), for an integer x,
// where the polynomial f is given as for fmpq_poly_set_str.
STATUS test_is_square_of(const char* f_str, slong x_int, bool is_square)
{
  STATUS ret;
  fmpq_poly_t f;
  nf_t K;
  nf_elem_t x;

  fmpq_poly_init(f);
  fmpq_poly_set_str(f, f_str);
  nf_init(K, f);
  nf_elem_init(x, K);
  nf_elem_set_si(x, x_int, K);

  ret = (nf_elem_is_square_fast(x, K) == is_square) ? SUCCESS : FAIL;

  if (ret == FAIL)
    fprintf(stderr, "FAIL: %ld was found %sto be a square in the field defined by %s\n",
	    x_int, is_square ? "not " : "", f_str);

  nf_elem_clear(x, K);
  nf_clear(K);
  fmpq_poly_clear(f);

  return ret;
}

STATUS test_is_square()
{
  STATUS ret;

  // a^2 - 3 has a repeated factor modulo 3, where Z[a] is maximal.
  // Modulo 3 is the only place where 2 is seen not to be a square,
  // among the primes which are tested.
  ret = test_is_square_of("3  -3 0 1", 2, false) << 1;

  // a^2 + 9 has a repeated factor modulo 3, where Z[a] is not maximal.
  // Here -4 = (2a/3)^2 is not a square modulo 3.
  ret |= test_is_square_of("3  9 0 1", -4, true);

  return ret;
}

// The entry which follows the coefficients, to see that it is not written to.
#define PARSE_GUARD 123456789

STATUS test_parse_matrix_of(const char* mat_str, bool is_valid)
{
  STATUS ret;
  Z64 Q_coeffs[16];
  int i;
  bool parsed;

  Q_coeffs[15] = PARSE_GUARD;
  parsed = parse_matrix(Q_coeffs, mat_str);

  ret = SUCCESS;
  if (parsed != is_valid) {
    fprintf(stderr, "FAIL: the matrix %s was found %sto be valid\n",
	    mat_str, is_valid ? "not " : "");
    ret = FAIL;
  }
  if (Q_coeffs[15] != PARSE_GUARD) {
    fprintf(stderr, "FAIL: more than 15 coefficients were read from %s\n",
	    mat_str);
    ret = FAIL;
  }
  // The coefficients of the valid matrices below are -7,...,7.
  if (parsed && is_valid)
    for (i = 0; i < 15; i++)
      if (Q_coeffs[i] != i - 7) {
	fprintf(stderr, "FAIL: coefficient %d of %s was read as %" PRId64 "\n",
		i, mat_str, Q_coeffs[i]);
	ret = FAIL;
      }

  return ret;
}

STATUS test_parse_matrix()
{
  STATUS ret;
  Z64 Q_coeffs[15];

  ret = test_parse_matrix_of("-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7", true);
  // brackets and spaces are allowed
  ret |= test_parse_matrix_of("[-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7]", true);
  ret |= test_parse_matrix_of(" -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7 ", true);
  ret |= test_parse_matrix_of(" [ -7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7 ] ", true);
  ret |= test_parse_matrix_of("[-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7", false);
  ret |= test_parse_matrix_of("-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7]", false);
  // too many coefficients, which used to be written past the end of the array
  ret |= test_parse_matrix_of("-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8", false);
  ret |= test_parse_matrix_of("-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,", false);
  // too few coefficients
  ret |= test_parse_matrix_of("-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6", false);
  ret |= test_parse_matrix_of("", false);
  ret |= test_parse_matrix_of("-7,,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7", false);
  // coefficients which are not integers, which used to be read as integers
  ret |= test_parse_matrix_of("-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,x", false);
  ret |= test_parse_matrix_of("-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7x", false);
  ret |= test_parse_matrix_of("-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,99999999999999999999", false);

  // a coefficient which does not fit in 32 bits, which used to be truncated
  if ((!parse_matrix(Q_coeffs, "1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1099511627776")) ||
      (Q_coeffs[14] != -INT64_C(1099511627776))) {
    fprintf(stderr, "FAIL: the coefficient -2^40 was not read as a 64-bit integer\n");
    ret = FAIL;
  }

  return ret;
}

// Checks the mass formula for the lattice whose gram matrix is given by Q_coeffs
STATUS test_mass_of(const Z64* Q_coeffs, slong mass_num, slong mass_denom)
{
  square_matrix_t Q;
  fmpq_t mass, expected;
  STATUS ret;

  fmpq_init(mass);
  fmpq_init(expected);
  square_matrix_init_set_symm(Q, Q_coeffs, "A");

  get_mass(mass, Q);
  fmpq_set_si(expected, mass_num, mass_denom);

  ret = (fmpq_equal(mass, expected) ? SUCCESS : FAIL);

  square_matrix_clear(Q);
  fmpq_clear(expected);
  fmpq_clear(mass);

  return ret;
}

// The mass of each of these lattices is the sum of the reciprocals of the
// orders of the automorphism groups of the lattices in its genus, as computed
// by Magma. Together they cover lattices that are odd or even, maximal at 2 or
// not, with one, two, three and five Jordan constituents at 2, and with
// constituents of type I of even rank.
STATUS test_mass()
{
  STATUS ret;

  // Z^5, A_5 and D_5
  Z64 Z5[15] = {1,0,1,0,0,1,0,0,0,1,0,0,0,0,1};
  Z64 A5[15] = {2,-1,2,0,-1,2,0,0,-1,2,0,0,0,-1,2};
  Z64 D5[15] = {2,-1,2,0,-1,2,0,0,-1,2,0,0,-1,0,2};
  // lattices of discriminants 12, 48, 768 and 385
  Z64 Q12[15] = {2,-1,2,0,0,2,0,0,1,2,-1,0,-1,-1,4};
  Z64 Q48[15] = {2,0,2,0,-1,2,1,0,0,2,1,-1,1,1,12};
  Z64 Q768[15] = {2,1,2,0,0,2,0,0,-1,2,1,0,-1,0,172};
  Z64 Q385[15] = {2,1,2,1,1,2,1,0,1,8,-1,-1,-1,-4,30};
  // diag(1,2,4,8,16), A_2 + 2 A_2 + <4> and diag(1,3,3,9,15)
  Z64 diag2[15] = {1,0,2,0,0,4,0,0,0,8,0,0,0,0,16};
  Z64 A2A2[15] = {2,-1,2,0,0,4,0,0,-2,4,0,0,0,0,4};
  Z64 diag3[15] = {1,0,3,0,0,3,0,0,0,9,0,0,0,0,15};
  // diag(1,1,4,4,20), diag(1,3,3,3,4) and [[4,3],[3,4]] + diag(5,3,5),
  // where the constituents of type I are free, and of ranks 2 to 5
  Z64 free2[15] = {1,0,1,0,0,4,0,0,0,4,0,0,0,0,20};
  Z64 free4[15] = {1,0,3,0,0,3,0,0,0,3,0,0,0,0,4};
  Z64 free5[15] = {4,3,4,0,0,5,0,0,0,3,0,0,0,0,5};

  ret = test_mass_of(Z5, 1, 3840);
  ret |= test_mass_of(A5, 1, 1440);
  ret |= test_mass_of(D5, 1, 3840);
  ret |= test_mass_of(Q12, 1, 144);
  ret |= test_mass_of(Q48, 1, 9);
  ret |= test_mass_of(Q768, 256, 9);
  ret |= test_mass_of(Q385, 3, 1);
  ret |= test_mass_of(diag2, 1, 4);
  ret |= test_mass_of(A2A2, 5, 576);
  ret |= test_mass_of(diag3, 65, 576);
  ret |= test_mass_of(free2, 13, 384);
  ret |= test_mass_of(free4, 5, 192);
  ret |= test_mass_of(free5, 65, 16);

  return ret;
}

STATUS compute_eigenvectors(const genus_t genus, bool nonlifts, int num_idxs, const int* idxs)
{
  example_evs_t dummy;
  STATUS ret;
  int num_ps[2] = {0,0};

  example_evs_init(dummy, genus->num_conductors, NULL, num_ps, NULL, NULL);
  
  ret = test_evs(genus, dummy, nonlifts, num_idxs, idxs);
  example_evs_clear(dummy, genus->num_conductors);

  return ret;
}

STATUS compute_eigenvalues(const genus_t genus, int p, bool nonlifts, int num_idxs, const int* idxs)
{
  example_evs_t dummy;
  STATUS ret;
  int num_ps[2] = {1,0};
  const int* ps[2] = {
    (int []){p},
    (int []){}
  };

  example_evs_init(dummy, genus->num_conductors, NULL, num_ps, ps, NULL);
  
  ret = test_evs(genus, dummy, nonlifts, num_idxs, idxs);
  example_evs_clear(dummy, genus->num_conductors);

  return ret;
}

STATUS compute_eigenvalues_up_to(const genus_t genus, int prec, bool nonlifts, int num_idxs, const int* idxs)
{
  example_evs_t dummy;
  STATUS ret;
  int num_ps[2];
  int* ps[2];
  
  num_ps[0] = primes_up_to_prime_to(&(ps[0]), prec, fmpz_get_si(genus->disc));
  num_ps[1] = primes_up_to_prime_to(&(ps[1]), floor(sqrt(prec)), fmpz_get_si(genus->disc));

  // This doens't work in general at the moment - figure out why and how to fix it
  
  // num_ps[0] = primes_up_to_prime_to(&(ps[0]), prec, 1);
  // num_ps[1] = primes_up_to_prime_to(&(ps[1]), floor(sqrt(prec)), 1);
  
  example_evs_init(dummy, genus->num_conductors, NULL, num_ps, (const int**)ps, NULL);
  ret = test_evs(genus, dummy, nonlifts, num_idxs, idxs);
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

// The index of the conductor c among the conductors of the genus.
// Exits with an error if c is not one of them.
static int conductor_index(const genus_t genus, slong c)
{
  int c_idx;

  for (c_idx = 0; (c_idx < genus->num_conductors) && (genus->conductors[c_idx] != c);) c_idx++;

  if (c_idx == genus->num_conductors) {
    fprintf(stderr, "Error: -cond=%ld is not one of the possible conductors\n", c);
    exit(EXIT_FAILURE);
  }

  return c_idx;
}

STATUS compute_hecke_col(const genus_t genus, slong p, slong c)
{
  int* hecke;
  int c_idx, i, gen_idx;

  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;  

  c_idx = conductor_index(genus, c);

  // The column is that of the first lattice which contributes to the space of
  // this conductor, or of the first lattice if the space has dimension 0.
  gen_idx = 0;
  if (genus->dims[c_idx] > 0)
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
  
  printf("{%ld : ", genus->conductors[c_idx]);
  printf("[");
  for (i = 0; i < genus->dims[c_idx]; i++) {
    printf("%d", hecke[i]);
    if (i != genus->dims[c_idx]-1)
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
  int c_idx;

  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;  

  c_idx = conductor_index(genus, c);

  // right now we only handle differently the trivial spinor norm case
  if (c_idx != 0)
    return compute_hecke_matrix_all_conds(genus, p);

  cpuclock_0 = clock();

  hecke = hecke_matrix(genus, p);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;

  //printf("{%ld : {", p);
  printf("{%ld : ", genus->conductors[c_idx]);
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
