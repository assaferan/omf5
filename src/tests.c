#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include <flint/nmod_mat.h>

#include "arith.h"
#include "eigenvalues.h"
#include "genus.h"
#include "hecke.h"
#include "io.h"
#include "mass.h"
#include "matrix_tools.h"
#include "spinor.h"
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

// The matrices of T_{p^k} on all the conductors, acting on row vectors.
static fmpq_mat_t* hecke_fmpq_mats_all_conductors(const genus_t genus, int p, int k)
{
  matrix_TYP** hecke;
  fmpq_mat_t* hecke_fmpq;
  slong c;

  // These are found by the functions which find the eigenvalues at the
  // conductors other than 1. An eigenvalue which is printed is found again,
  // from one column, and at the conductor 1 by other functions (hecke_col,
  // hecke_col_nbr_data): what is checked on these matrices is not that column.
  if (k == 1)
    hecke = hecke_matrices_all_conductors(genus, p);
  else
    hecke = hecke_matrices_nbr_data_all_conductors(genus, p, k);

  hecke_fmpq = (fmpq_mat_t*)malloc(genus->num_conductors * sizeof(fmpq_mat_t));
  for (c = 0; c < genus->num_conductors; c++) {
    fmpq_mat_init_set_matrix_TYP(hecke_fmpq[c], hecke[c]);
    // we transpose to match with our conventions
    fmpq_mat_transpose(hecke_fmpq[c], hecke_fmpq[c]);
    free_mat(hecke[c]);
  }
  free(hecke);

  return hecke_fmpq;
}

// bad_hecke has an entry for every prime in ps. It is NULL until the matrices
// of T_{p^k} on all the conductors are needed, and they are kept in it then,
// for the next forms. The caller clears them.
STATUS test_eigenvalues(const genus_t genus, const eigenvalues_t evs,
			int num_evs, int form_idx, const int* ps,
			const int* test_evs, int k, int c, bool ev_print,
			fmpq_mat_t** bad_hecke)
{
  int i;
  nf_elem_t ev;
  fmpq_t trace;
  bool* skip;
#ifdef DEBUG_LEVEL_FULL
  int j;
  matrix_TYP* hecke;
#endif // DEBUG_LEVEL_FULL

  fmpq_init(trace);

  // When no eigenvector was found for the form, the vector we have is an
  // eigenvector of the last Hecke operator of the decomposition, in a space
  // which the operators of the decomposition, at primes not dividing the
  // discriminant, do not split. It need not be an eigenvector of the operator
  // at a prime dividing the discriminant, and the number found from it then
  // depends on the vector, which depends on the basis of the lattice. It does
  // not when the operator acts on the space of the form as an element of the
  // field of the form: every vector which could have been kept is then an
  // eigenvector, with that eigenvalue. So the eigenvalue is given at such a
  // prime, and at no other prime dividing the discriminant. This takes the
  // matrix of the operator, and not one column of it.
  skip = (bool*)malloc((num_evs + 1) * sizeof(bool));
  for (i = 0; i < num_evs; i++) {
    skip[i] = false;
    if (evs->is_old[form_idx] && (fmpz_fdiv_ui(genus->disc, ps[i]) == 0)) {
      if (bad_hecke[i] == NULL)
	bad_hecke[i] = hecke_fmpq_mats_all_conductors(genus, ps[i], k);
      skip[i] = !ev_is_scalar_on_old_space(evs, form_idx, bad_hecke[i][c]);
      if (skip[i])
	fprintf(stderr, "Warning: no eigenvalue of T_{%d^%d} is given for the form %d "
		"of conductor %ld, for which no eigenvector was found: "
		"%d divides the discriminant, and T_{%d^%d} is not a scalar "
		"on the space of the form\n",
		ps[i], k, form_idx, (slong)genus->conductors[c], ps[i], ps[i], k);
    }
  }

  fprintf(stderr, "traces of hecke eigenvalues of T_{p^%d} are:\n", k);

  if (ev_print) {
    if (k == 1)
      printf("'lambda_p' : [");
    if (k == 2)
      printf("'lambda_p_square' : [");
  }
  
  for (i = 0; i < num_evs; i++) {
    // See the warning above.
    if (skip[i])
      continue;

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
	free(skip);
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
  free(skip);
  return SUCCESS;
}

// Clear the matrices which test_eigenvalues kept, for the primes ps[k] of ex
static void bad_hecke_clear(fmpq_mat_t** bad_hecke[2], const genus_t genus, const example_evs_t ex)
{
  int k, i;
  slong c;

  for (k = 0; k < 2; k++) {
    for (i = 0; i < ex->num_ps[k]; i++)
      if (bad_hecke[k][i] != NULL) {
	for (c = 0; c < genus->num_conductors; c++)
	  fmpq_mat_clear(bad_hecke[k][i][c]);
	free(bad_hecke[k][i]);
      }
    free(bad_hecke[k]);
  }

  return;
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
  // The matrices at the primes dividing the discriminant, see test_eigenvalues
  fmpq_mat_t** bad_hecke[2];

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

  for (k = 0; k < 2; k++)
    bad_hecke[k] = (fmpq_mat_t**)calloc(ex->num_ps[k] + 1, sizeof(fmpq_mat_t*));
    
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
			       ex->ps[k], test_evs, k+1, c, true,
			       bad_hecke[k]) == FAIL) {
	    bad_hecke_clear(bad_hecke, genus, ex);
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

  bad_hecke_clear(bad_hecke, genus, ex);
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

// Tests the two functions which find whether x is a square in Q[a] / (f).
// The polynomial f, and x as a polynomial in a, are given as for
// fmpq_poly_set_str. nf_elem_is_square has to find what is_square says.
// nf_elem_is_square_fast has to find what is_square_fast says: it can
// find x to be a square when it is not, and never the other way.
STATUS test_is_square_of(const char* f_str, const char* x_str, bool is_square, bool is_square_fast)
{
  STATUS ret;
  fmpq_poly_t f, x_poly;
  nf_t K;
  nf_elem_t x;

  fmpq_poly_init(f);
  fmpq_poly_init(x_poly);
  fmpq_poly_set_str(f, f_str);
  fmpq_poly_set_str(x_poly, x_str);
  nf_init(K, f);
  nf_elem_init(x, K);
  nf_elem_set_fmpq_poly(x, x_poly, K);

  ret = SUCCESS;

  if (nf_elem_is_square(x, K) != is_square) {
    ret = FAIL;
    fprintf(stderr, "FAIL: %s was found %sto be a square in the field defined by %s\n",
	    x_str, is_square ? "not " : "", f_str);
  }

  if (nf_elem_is_square_fast(x, K) != is_square_fast) {
    ret = FAIL;
    fprintf(stderr, "FAIL: %s was found %sto be a square modulo small primes in the field defined by %s\n",
	    x_str, is_square_fast ? "not " : "", f_str);
  }

  nf_elem_clear(x, K);
  nf_clear(K);
  fmpq_poly_clear(x_poly);
  fmpq_poly_clear(f);

  return ret;
}

STATUS test_is_square()
{
  STATUS ret;

  // a^2 - 3 has a repeated factor modulo 3, where Z[a] is maximal.
  // Modulo 3 is the only place where 2 is seen not to be a square,
  // among the primes which are tested.
  ret = test_is_square_of("3  -3 0 1", "1  2", false, false);

  // a^2 + 9 has a repeated factor modulo 3, where Z[a] is not maximal.
  // Here -4 = (2a/3)^2 is not a square modulo 3.
  ret |= test_is_square_of("3  9 0 1", "1  -4", true, true);

  // 0 is a square.
  ret |= test_is_square_of("2  0 1", "1  0", true, true);
  ret |= test_is_square_of("3  -2 0 1", "1  0", true, true);

  // In the field of the rational numbers.
  ret |= test_is_square_of("2  0 1", "1  9/4", true, true);
  ret |= test_is_square_of("2  0 1", "1  2", false, false);
  // No prime is tested for -1.
  ret |= test_is_square_of("2  0 1", "1  -1", false, true);

  // In Q[a] / (a^2 - 2), 2 is a square, and it is not one in the field
  // which it generates. 3 + 2a = (1 + a)^2, and 1 + a has the norm -1.
  ret |= test_is_square_of("3  -2 0 1", "1  2", true, true);
  ret |= test_is_square_of("3  -2 0 1", "2  3 2", true, true);
  ret |= test_is_square_of("3  -2 0 1", "2  1 1", false, true);
  ret |= test_is_square_of("3  -2 0 1", "1  3", false, false);

  // In Q[a] / (a^4 + 1), -1 = (a^2)^2 and 2 = (a + a^(-1))^2.
  // 3 is not a square, and it is one modulo 3, 5 and 7, the odd primes
  // which are tested: a^4 + 1 has two factors of degree 2 modulo each.
  // It is not a square modulo 17.
  ret |= test_is_square_of("5  1 0 0 0 1", "1  -1", true, true);
  ret |= test_is_square_of("5  1 0 0 0 1", "1  2", true, true);
  ret |= test_is_square_of("5  1 0 0 0 1", "1  3", false, true);

  // In Q[a] / (a^3 - 2), a is not a square, and no prime is tested for it.
  ret |= test_is_square_of("4  -2 0 0 1", "3  0 0 1", true, true);
  ret |= test_is_square_of("4  -2 0 0 1", "2  0 1", false, true);

  // Elements with a denominator. In Q[a] / (a^2 - 2), 9/4 + a = (1/2 + a)^2,
  // and (1 + a)/2 has the norm -1/4.
  ret |= test_is_square_of("3  -2 0 1", "2  9/4 1", true, true);
  ret |= test_is_square_of("3  -2 0 1", "2  1/2 1/2", false, false);

  // a^2 - 5 has a repeated factor modulo 2, where Z[a] is not maximal:
  // (3 + a)/2 is the square of (1 + a)/2, which is not in Z[a],
  // and which has the norm -1.
  ret |= test_is_square_of("3  -5 0 1", "2  3/2 1/2", true, true);
  ret |= test_is_square_of("3  -5 0 1", "2  1/2 1/2", false, false);
  ret |= test_is_square_of("3  -5 0 1", "1  2", false, false);

  // A polynomial which is not monic. In Q[a] / (2a^2 - 1), 2 = (2a)^2,
  // and a has the norm -1/2.
  ret |= test_is_square_of("3  -1 0 2", "1  2", true, true);
  ret |= test_is_square_of("3  -1 0 2", "2  0 1", false, true);

  // In Q[a] / (a^5 - a - 1), 1 + a^2 is not a square.
  ret |= test_is_square_of("6  -1 -1 0 0 0 1", "5  1 0 2 0 1", true, true);
  ret |= test_is_square_of("6  -1 -1 0 0 0 1", "3  1 0 1", false, true);

  // In Q[a] / (a^6 - 2), a^2 is a square, and it is not one in the field
  // of degree 3 which it generates. 2 = (a^3)^2, and a^3 is not a square.
  ret |= test_is_square_of("7  -2 0 0 0 0 0 1", "3  0 0 1", true, true);
  ret |= test_is_square_of("7  -2 0 0 0 0 0 1", "1  2", true, true);
  ret |= test_is_square_of("7  -2 0 0 0 0 0 1", "4  0 0 0 1", false, true);
  ret |= test_is_square_of("7  -2 0 0 0 0 0 1", "1  3", false, false);

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

// Whether s1 / denom1 and s2 / denom2 are the same matrix. The products are
// found with big integers: a product which has wrapped around is not the same
// as one which has not, though the two can agree modulo 2^64.
static bool rational_matrix_is_equal(const square_matrix_t s1, Z64 denom1,
				     const square_matrix_t s2, Z64 denom2)
{
  fmpz_t lhs, rhs;
  bool is_equal;
  int i, j;

  fmpz_init(lhs);
  fmpz_init(rhs);
  is_equal = true;
  for (i = 0; is_equal && (i < QF_RANK); i++)
    for (j = 0; is_equal && (j < QF_RANK); j++) {
      fmpz_set_si(lhs, s1[i][j]);
      fmpz_mul_si(lhs, lhs, denom2);
      fmpz_set_si(rhs, s2[i][j]);
      fmpz_mul_si(rhs, rhs, denom1);
      is_equal = fmpz_equal(lhs, rhs);
    }
  fmpz_clear(lhs);
  fmpz_clear(rhs);

  return is_equal;
}

static bool isometry_is_same(const isometry_t isom1, const isometry_t isom2)
{
  return rational_matrix_is_equal(isom1->s, isom1->denom, isom2->s, isom2->denom) &&
    rational_matrix_is_equal(isom1->s_inv, isom1->inv_denom, isom2->s_inv, isom2->inv_denom);
}

// Sets isom to num / denom, written as the fraction (c * num) / (c * denom),
// and its inverse as a fraction with c times its numerators and denominator.
static void isometry_init_set_scaled(isometry_t isom, const square_matrix_t num, Z64 denom, Z64 c)
{
  isometry_init_set_square_matrix(isom, num, denom);
  square_matrix_mul_scalar(isom->s, isom->s, c);
  isom->denom *= c;
  square_matrix_mul_scalar(isom->s_inv, isom->s_inv, c);
  isom->inv_denom *= c;

  return;
}

// Checks the products of two isometries which are written as fractions with
// the common factor c, against their products when c = 1.
STATUS test_isometry_mul_of(Z64 c)
{
  // three times an isometry of denominator 3, and a matrix of determinant 1
  square_matrix_t M = {{3,0,0,-3,3},{2,3,-2,1,1},{-1,0,4,1,-2},{1,0,-1,2,-1},{-1,0,1,1,1}};
  square_matrix_t N = {{1,2,3,0,0},{0,1,2,0,0},{0,0,1,5,0},{0,0,0,1,7},{0,0,0,0,1}};
  isometry_t a, b, a1, b1, ab, bb, prod;
  bool is_same;

  isometry_init_set_scaled(a, M, 3, c);
  isometry_init_set_scaled(b, N, 1, c);
  isometry_init_set_scaled(a1, M, 3, 1);
  isometry_init_set_scaled(b1, N, 1, 1);
  isometry_init(ab);
  isometry_init(bb);
  isometry_mul(ab, a1, b1);
  isometry_mul(bb, b1, b1);

  isometry_init(prod);
  isometry_mul(prod, a, b);
  is_same = isometry_is_same(prod, ab);
  isometry_mul(prod, b, b);
  is_same = is_same && isometry_is_same(prod, bb);
  isometry_clear(prod);

  // in place
  isometry_init_set(prod, a);
  isometry_muleq_right(prod, b);
  is_same = is_same && isometry_is_same(prod, ab);
  isometry_clear(prod);

  isometry_init_set(prod, b);
  isometry_muleq_left(prod, a);
  is_same = is_same && isometry_is_same(prod, ab);
  isometry_clear(prod);

  // by a matrix
  isometry_init_set(prod, b);
  isometry_mul_mateq_left(prod, a->s, a->denom);
  is_same = is_same && isometry_is_same(prod, ab);
  isometry_clear(prod);

  isometry_init_set(prod, a);
  isometry_mul_mateq_right(prod, b->s, b->denom);
  is_same = is_same && isometry_is_same(prod, ab);
  isometry_clear(prod);

  if (!is_same)
    fprintf(stderr, "FAIL: a product of isometries written with the common factor %" PRId64 " is wrong\n", c);

  isometry_clear(a);
  isometry_clear(b);
  isometry_clear(a1);
  isometry_clear(b1);
  isometry_clear(ab);
  isometry_clear(bb);

  return (is_same ? SUCCESS : FAIL);
}

// The numerators and the denominators of isometries are multiplied when the
// isometries are composed. The denominators used to be of type int, and
// neither product was checked.
STATUS test_isometry_mul()
{
  STATUS ret;

  ret = test_isometry_mul_of(1);
  // the denominators of the products do not fit in an int
  ret |= test_isometry_mul_of(65537);
  // The products do not fit in 64 bits before the fractions are reduced. With
  // 2^31 the product of the denominators of the second matrix with itself
  // does fit, and the numerators do not.
  ret |= test_isometry_mul_of(INT64_C(2147483648));
  ret |= test_isometry_mul_of(INT64_C(4294967311));

  return ret;
}

// Checks the vector found to span the radical modulo 4 of the lattice whose
// Gram matrix is given by Q_coeffs, for each of the 120 orders of the vectors
// of the basis. The radical modulo 2 of these lattices is a line.
STATUS test_radical_mod_4_of(const Z64* Q_coeffs)
{
  square_matrix_t Q;
  fmpz_mat_t q, perm_q;
  spinor_t spinor;
  nmod_mat_t perm_q_4, rad_q;
  int perm[QF_RANK];
  int code, digits, i, j;
  slong prime_idx;
  bool is_perm, is_odd;
  STATUS ret;

  ret = SUCCESS;

  square_matrix_init_set_symm(Q, Q_coeffs, "A");
  fmpz_mat_init_set_square_matrix(q, Q);
  fmpz_mat_init(perm_q, QF_RANK, QF_RANK);
  nmod_mat_init(rad_q, 1, QF_RANK, 4);

  // The orders are the numbers with QF_RANK different digits in base QF_RANK.
  for (code = 0; code < 3125; code++) {
    digits = code;
    for (i = 0; i < QF_RANK; i++) {
      perm[i] = digits % QF_RANK;
      digits /= QF_RANK;
    }
    is_perm = true;
    for (i = 0; i < QF_RANK; i++)
      for (j = 0; j < i; j++)
	if (perm[i] == perm[j])
	  is_perm = false;
    if (!is_perm)
      continue;

    for (i = 0; i < QF_RANK; i++)
      for (j = 0; j < QF_RANK; j++)
	fmpz_set(fmpz_mat_entry(perm_q, i, j), fmpz_mat_entry(q, perm[i], perm[j]));

    spinor_init_fmpz_mat(spinor, perm_q);

    for (prime_idx = 0; prime_idx < spinor->num_primes; prime_idx++)
      if (spinor->primes[prime_idx].n == 4)
	break;

    if ((prime_idx == spinor->num_primes) || (nmod_mat_nrows(spinor->rads[prime_idx]) != 1))
      ret = FAIL;
    else {
      // The vector is not 0 modulo 2, and its products with the vectors of
      // the basis are 0 modulo 4.
      is_odd = false;
      for (i = 0; i < QF_RANK; i++)
	if (nmod_mat_entry(spinor->rads[prime_idx], 0, i) % 2 != 0)
	  is_odd = true;
      nmod_mat_init_set_fmpz_mat(perm_q_4, perm_q, 4);
      nmod_mat_mul(rad_q, spinor->rads[prime_idx], perm_q_4);
      if ((!is_odd) || (!nmod_mat_is_zero(rad_q)))
	ret = FAIL;
      nmod_mat_clear(perm_q_4);
    }

    spinor_clear(spinor);
  }

  nmod_mat_clear(rad_q);
  fmpz_mat_clear(perm_q);
  fmpz_mat_clear(q);
  square_matrix_clear(Q);

  return ret;
}

// The vector found for the radical modulo 4 could be wrong when the first
// column of the Gram matrix is 0 modulo 2, as it is for the last two lattices
// here, and was wrong for both. The second lattice of discriminant 124 is the
// first one in another basis.
STATUS test_radical_mod_4()
{
  STATUS ret;

  // D_5, of discriminant 2, and a lattice of discriminant 12
  Z64 D5[15] = {2,-1,2,0,-1,2,0,0,-1,2,0,0,-1,0,2};
  Z64 Q12[15] = {2,-1,2,0,0,2,0,0,1,2,-1,0,-1,-1,4};
  // lattices of discriminants 124, 124 and 116; the first vector of the basis
  // of the last two spans the radical modulo 2
  Z64 Q124[15] = {2,0,2,1,1,4,0,0,0,4,-1,0,-1,-1,6};
  Z64 Q124_even[15] = {8,-2,2,0,1,4,-6,2,1,6,0,0,-1,-1,6};
  Z64 Q116_even[15] = {8,2,2,4,1,6,0,0,1,2,-2,0,-1,1,4};

  ret = test_radical_mod_4_of(D5);
  ret |= test_radical_mod_4_of(Q12);
  ret |= test_radical_mod_4_of(Q124);
  ret |= test_radical_mod_4_of(Q124_even);
  ret |= test_radical_mod_4_of(Q116_even);

  return ret;
}

// Whether the two matrices of the given dimension commute.
static bool hecke_matrices_commute(const matrix_TYP* a, const matrix_TYP* b, slong dim)
{
  slong i, j, k, ab, ba;

  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) {
      ab = ba = 0;
      for (k = 0; k < dim; k++) {
	ab += (slong)(a->array.SZ[i][k]) * b->array.SZ[k][j];
	ba += (slong)(b->array.SZ[i][k]) * a->array.SZ[k][j];
      }
      if (ab != ba)
	return false;
    }

  return true;
}

// Checks the Hecke operator T_p when k = 1, and T_{p^2} when k = 2, at the prime
// p for the lattice given by Q_coeffs in the format GG: on each of the spaces,
// in the order of the conductors, it commutes with the operator at the prime q,
// which does not divide the discriminant, it is self-adjoint for the weights
// given by the numbers of automorphisms, and it has the given trace. The
// operator T_{p^2} is also checked to commute with T_p.
static STATUS test_hecke_bad_prime_power_of(const Z64* Q_coeffs, int p, int k, int q,
					    slong num_conds, const slong* traces)
{
  genus_t genus;
  matrix_TYP** hecke_p;
  matrix_TYP** hecke_q;
  matrix_TYP** hecke_p1;
  slong c, i, j, dim, trace;
  bool is_ok, self_adjoint;
  const char* name = (k == 1) ? "" : " squared";

  compute_genus(genus, Q_coeffs, "GG");

  if (genus->num_conductors != num_conds) {
    fprintf(stderr, "FAIL: there are %ld conductors and not %ld\n", genus->num_conductors, num_conds);
    genus_clear(genus);
    return FAIL;
  }

  hecke_q = hecke_matrices_all_conductors(genus, q);
  hecke_p1 = hecke_matrices_all_conductors(genus, p);
  hecke_p = (k == 1) ? hecke_p1 : hecke_matrices_nbr_data_all_conductors(genus, p, k);

  is_ok = true;
  for (c = 0; c < genus->num_conductors; c++) {
    dim = genus->dims[c];
    trace = 0;
    for (i = 0; i < dim; i++)
      trace += hecke_p[c]->array.SZ[i][i];
    if (trace != traces[c]) {
      fprintf(stderr, "FAIL: the trace of the Hecke operator at %d%s for the conductor %ld is %ld and not %ld\n",
	      p, name, genus->conductors[c], trace, traces[c]);
      is_ok = false;
    }
    if (!hecke_matrices_commute(hecke_p[c], hecke_q[c], dim)) {
      fprintf(stderr, "FAIL: the Hecke operators at %d%s and %d for the conductor %ld do not commute\n",
	      p, name, q, genus->conductors[c]);
      is_ok = false;
    }
    if ((k != 1) && !hecke_matrices_commute(hecke_p[c], hecke_p1[c], dim)) {
      fprintf(stderr, "FAIL: the Hecke operators at %d%s and %d for the conductor %ld do not commute\n",
	      p, name, p, genus->conductors[c]);
      is_ok = false;
    }
    self_adjoint = true;
    for (i = 0; i < dim; i++)
      for (j = 0; j < dim; j++)
	if ((slong)(hecke_p[c]->array.SZ[i][j]) * genus->num_auts[c][j] !=
	    (slong)(hecke_p[c]->array.SZ[j][i]) * genus->num_auts[c][i])
	  self_adjoint = false;
    if (!self_adjoint) {
      fprintf(stderr, "FAIL: the Hecke operator at %d%s for the conductor %ld is not self-adjoint\n",
	      p, name, genus->conductors[c]);
      is_ok = false;
    }
  }

  for (c = 0; c < genus->num_conductors; c++) {
    free_mat(hecke_p1[c]);
    free_mat(hecke_q[c]);
    if (k != 1)
      free_mat(hecke_p[c]);
  }
  free(hecke_p1);
  free(hecke_q);
  if (k != 1)
    free(hecke_p);

  genus_clear(genus);

  return (is_ok ? SUCCESS : FAIL);
}

STATUS test_hecke_bad_prime_of(const Z64* Q_coeffs, int p, int q, slong num_conds, const slong* traces)
{
  return test_hecke_bad_prime_power_of(Q_coeffs, p, 1, q, num_conds, traces);
}

static STATUS test_hecke_bad_prime_square_of(const Z64* Q_coeffs, int p, int q, slong num_conds, const slong* traces)
{
  return test_hecke_bad_prime_power_of(Q_coeffs, p, 2, q, num_conds, traces);
}

// The neighbors at a prime p dividing the discriminant are images of the lattice
// under isometries which are not integral at p, and the character of a conductor
// divisible by p was not found correctly at them, for example when p = 2,
// or when p^2 divides the discriminant. The traces are those of the operators
// HeckeOperator(M, p) of the package ModFrmAlg of Magma, for the spaces
// M := OrthogonalModularForms(L : d := d, Special := true).
// From discriminant 500 on, the radical of the lattice modulo a prime is not a
// line. Magma does not find the spaces of such lattices, and the traces are
// those found here.
STATUS test_hecke_bad_primes()
{
  // discriminant 66, with the conductors 1,2,3,6,11,22,33,66
  Z64 Q_66[15] = {1,1,-1,1,1,1,-1,1,1,1,-1,-1,3,3,4};
  slong traces_66_2[8] = {20,0,0,0,0,-5,3,0};
  slong traces_66_3[8] = {61,0,0,-5,0,6,-20,0};
  // discriminant 75, with the conductors 1,3,5,15
  Z64 Q_75[15] = {1,0,-1,-1,0,1,-1,0,0,2,0,1,2,0,2};
  slong traces_75_3[4] = {57,0,0,-20};
  slong traces_75_5[4] = {320,0,0,60};
  // discriminant 76, with the conductors 1,2,19,38
  Z64 Q_76[15] = {1,-1,0,-1,0,1,0,1,0,1,0,1,3,2,3};
  slong traces_76_2[4] = {37,0,0,7};
  // discriminant 40, with the conductors 1,2,5,10
  Z64 Q_40[15] = {1,-1,0,-1,1,1,0,1,0,1,-1,0,2,1,3};
  slong traces_40_2[4] = {24,0,0,6};
  // discriminant 54, with the conductors 1,2,3,6
  Z64 Q_54[15] = {1,0,0,0,-1,1,0,1,-1,1,0,1,2,-2,3};
  slong traces_54_2[4] = {17,0,0,-5};
  slong traces_54_3[4] = {105,0,0,12};
  // discriminant 80, with the conductors 1,2,5,10
  Z64 Q_80[15] = {1,0,1,-1,0,1,-1,0,1,2,-2,-1,2,-1,3};
  slong traces_80_2[4] = {40,0,0,16};
  // discriminant 162, with the conductors 1,2,3,6
  Z64 Q_162[15] = {1,1,0,0,-1,2,1,1,1,2,-1,2,2,1,3};
  slong traces_162_3[4] = {231,0,0,51};
  // Discriminant 500, with the conductors 1,2,5,10. The radical modulo 5 is
  // a plane, and the space of conductor 5 is of dimension 3.
  Z64 Q_500[15] = {1,0,0,-1,0,2,-2,-1,-1,3,-2,3,3,0,4};
  slong traces_500_5[4] = {390,0,20,0};
  slong traces_500_3[4] = {105,0,1,0};
  // Discriminant 104, with the conductors 1,13. The radical modulo 2 is of
  // dimension 3, and there is no character at 2.
  Z64 Q_104[15] = {1,0,0,-1,-1,1,0,1,1,1,-1,1,2,1,6};
  slong traces_104_2[2] = {31,0};
  // Discriminant 125, with the conductors 1,5. The radical modulo 5 is of
  // dimension 3, and no vector outside it is isotropic modulo 5: there are
  // no 5-neighbors, and the Hecke operator at 5 is 0. One of the two lattices
  // of the genus has no isotropic vector modulo 5 whose first coordinate is 1.
  Z64 Q_125[15] = {1,-1,-1,1,0,1,1,0,0,2,-2,0,2,0,5};
  slong traces_125_5[2] = {0,0};
  // Discriminant 324, with the conductors 1,2,3,6. The same holds modulo 3,
  // and the neighbors at 3 are used to time the isometry tests when the genus
  // is found.
  Z64 Q_324[15] = {1,-1,0,0,0,1,0,0,0,3,3,-3,4,-4,4};
  slong traces_324_3[4] = {0,0,0,0};
  // Discriminant 2625, with the conductors 1,3,5,15,7,21,35,105. The same
  // holds modulo 5, and the spaces of the conductors 5,15,21,35,105 are not 0.
  Z64 Q_2625[15] = {1,1,0,0,0,1,0,0,0,5,0,0,5,5,10};
  slong traces_2625_5[8] = {0,0,0,0,0,0,0,0};
  slong traces_2625_3[8] = {113,0,-2,10,0,-50,-8,9};
  STATUS ret;

  // the prime divides the discriminant once
  ret = test_hecke_bad_prime_of(Q_66, 2, 5, 8, traces_66_2);
  ret |= test_hecke_bad_prime_of(Q_66, 3, 5, 8, traces_66_3);
  ret |= test_hecke_bad_prime_of(Q_75, 3, 2, 4, traces_75_3);
  ret |= test_hecke_bad_prime_of(Q_54, 2, 5, 4, traces_54_2);
  // to the power 2
  ret |= test_hecke_bad_prime_of(Q_75, 5, 2, 4, traces_75_5);
  ret |= test_hecke_bad_prime_of(Q_76, 2, 3, 4, traces_76_2);
  // to the power 3
  ret |= test_hecke_bad_prime_of(Q_40, 2, 3, 4, traces_40_2);
  ret |= test_hecke_bad_prime_of(Q_54, 3, 5, 4, traces_54_3);
  // to the power 4
  ret |= test_hecke_bad_prime_of(Q_80, 2, 3, 4, traces_80_2);
  ret |= test_hecke_bad_prime_of(Q_162, 3, 5, 4, traces_162_3);
  // The radical modulo 5 is not a line: the character at 5 is a determinant,
  // at the isometries to the 5-neighbors and at those which are integral at 5.
  ret |= test_hecke_bad_prime_of(Q_500, 5, 3, 4, traces_500_5);
  ret |= test_hecke_bad_prime_of(Q_500, 3, 7, 4, traces_500_3);
  // the radical modulo 2 is not a line
  ret |= test_hecke_bad_prime_of(Q_104, 2, 3, 2, traces_104_2);
  // there are no neighbors at the prime
  ret |= test_hecke_bad_prime_of(Q_125, 5, 3, 2, traces_125_5);
  ret |= test_hecke_bad_prime_of(Q_324, 3, 5, 4, traces_324_3);
  ret |= test_hecke_bad_prime_of(Q_2625, 5, 2, 8, traces_2625_5);
  ret |= test_hecke_bad_prime_of(Q_2625, 3, 2, 8, traces_2625_3);

  return ret;
}

// The operator T_{p^2} at a prime p dividing the discriminant is the sum over
// the neighbors N of the lattice L with pL in N and pN in L which are found
// from the planes isotropic modulo p that do not meet the radical. The
// isometries to them were given no character at p when p^2 divides the
// discriminant, and the number of the planes was not found when there are none.
// The traces of the first six lattices are those of HeckeOperator(M, p, 2) in
// the package ModFrmAlg of Magma, for the spaces
// M := OrthogonalModularForms(L : d := d, Special := true), and the others are
// those found here.
STATUS test_hecke_square_bad_primes()
{
  // discriminant 33, with the conductors 1,3,11,33
  Z64 Q_33[15] = {1,1,1,1,-1,1,1,1,-1,1,1,-1,2,0,3};
  slong traces_33_3[4] = {192,0,0,24};
  slong traces_33_11[4] = {0,0,0,0};
  // discriminant 66, with the conductors 1,2,3,6,11,22,33,66
  Z64 Q_66[15] = {1,1,-1,1,1,1,-1,1,1,1,-1,-1,3,3,4};
  slong traces_66_2[8] = {0,0,0,0,0,0,0,0};
  // discriminant 124, with the conductors 1,2,31,62
  Z64 Q_124[15] = {1,0,1,0,-1,1,1,0,0,2,0,-1,2,-1,3};
  slong traces_124_2[4] = {48,0,0,24};
  // discriminant 153, with the conductors 1,3,17,51
  Z64 Q_153[15] = {1,0,0,-1,-1,1,1,1,0,2,2,-1,3,1,3};
  slong traces_153_3[4] = {240,0,0,48};
  // discriminant 136, with the conductors 1,2,17,34
  Z64 Q_136[15] = {1,0,1,-1,0,1,1,0,-1,2,1,-1,3,-2,3};
  slong traces_136_2[4] = {52,0,0,16};
  // discriminant 189, with the conductors 1,3,7,21
  Z64 Q_189[15] = {1,0,0,1,0,1,0,0,-1,1,0,0,2,0,7};
  slong traces_189_3[4] = {270,0,0,90};
  // Discriminant 500, with the conductors 1,2,5,10. The radical modulo 5 is
  // a plane, and its complement has no plane which is isotropic modulo 5.
  Z64 Q_500[15] = {1,0,0,-1,0,2,-2,-1,-1,3,-2,3,3,0,4};
  slong traces_500_5[4] = {0,0,0,0};
  // Discriminant 192, with the conductors 1,3. The form is 0 modulo 2, and
  // there are no neighbors at 2. Counting them used to end in a division by 0.
  Z64 Q_192[15] = {2,2,2,-2,2,2,0,-2,0,2,-2,2,2,0,4};
  slong traces_192_2[2] = {0,0};
  // Discriminant 2625, with the conductors 1,3,5,15,7,21,35,105. The radical
  // modulo 5 is of dimension 3, and no vector outside it is isotropic modulo 5.
  // The spaces of the conductors 5,15,21,35,105 are not 0.
  Z64 Q_2625[15] = {1,1,0,0,0,1,0,0,0,5,0,0,5,5,10};
  slong traces_2625_5[8] = {0,0,0,0,0,0,0,0};
  STATUS ret;

  // The prime divides the discriminant once. The operator is 0 when the
  // lattice has no plane which is isotropic modulo the prime.
  ret = test_hecke_bad_prime_square_of(Q_33, 3, 2, 4, traces_33_3);
  ret |= test_hecke_bad_prime_square_of(Q_33, 11, 2, 4, traces_33_11);
  ret |= test_hecke_bad_prime_square_of(Q_66, 2, 5, 8, traces_66_2);
  // to the power 2
  ret |= test_hecke_bad_prime_square_of(Q_124, 2, 3, 4, traces_124_2);
  ret |= test_hecke_bad_prime_square_of(Q_153, 3, 2, 4, traces_153_3);
  // to the power 3
  ret |= test_hecke_bad_prime_square_of(Q_136, 2, 3, 4, traces_136_2);
  ret |= test_hecke_bad_prime_square_of(Q_189, 3, 2, 4, traces_189_3);
  // the radical modulo the prime is not a line
  ret |= test_hecke_bad_prime_square_of(Q_500, 5, 3, 4, traces_500_5);
  ret |= test_hecke_bad_prime_square_of(Q_2625, 5, 2, 8, traces_2625_5);
  // there are no neighbors at the prime
  ret |= test_hecke_bad_prime_of(Q_192, 2, 5, 2, traces_192_2);
  ret |= test_hecke_bad_prime_square_of(Q_192, 2, 5, 2, traces_192_2);

  return ret;
}

// Checks ev_is_scalar_on_old_space for the lattice given by Q_coeffs in the
// format GG, which has one form of conductor 1 for which no eigenvector is
// found: the operator T_{p^k} acts on the space of this form as an element of
// its field for the primes p = ps[i] and the powers k = ks[i] with
// is_scalar[i], and for no other. The identity is such an element, and a
// matrix which has one entry which is not 0 is not, since it is not 0 on the
// space and is not invertible on it.
static STATUS test_old_space_scalars_of(const Z64* Q_coeffs, int num, const int* ps,
					const int* ks, const bool* is_scalar)
{
  genus_t genus;
  eigenvalues_t* evs;
  fmpq_mat_t* hecke;
  fmpq_mat_t T;
  slong c, i, j, form_idx, num_old;
  bool is_ok;

  compute_genus(genus, Q_coeffs, "GG");
  evs = hecke_eigenforms_all_conductors(genus);

  num_old = 0;
  form_idx = -1;
  for (i = 0; i < evs[0]->num; i++)
    if (evs[0]->is_old[i]) {
      num_old++;
      form_idx = i;
    }

  is_ok = (num_old == 1);
  if (!is_ok)
    fprintf(stderr, "FAIL: there are %ld forms of conductor 1 for which no eigenvector was found, and not 1\n",
	    num_old);

  for (i = 0; is_ok && (i < num); i++) {
    hecke = hecke_fmpq_mats_all_conductors(genus, ps[i], ks[i]);
    if (ev_is_scalar_on_old_space(evs[0], form_idx, hecke[0]) != is_scalar[i]) {
      fprintf(stderr, "FAIL: the Hecke operator T_{%d^%d} %s a scalar on the space of the form\n",
	      ps[i], ks[i], is_scalar[i] ? "is not" : "is");
      is_ok = false;
    }
    for (c = 0; c < genus->num_conductors; c++)
      fmpq_mat_clear(hecke[c]);
    free(hecke);
  }

  if (is_ok) {
    fmpq_mat_init(T, evs[0]->dim, evs[0]->dim);
    fmpq_mat_one(T);
    if (!ev_is_scalar_on_old_space(evs[0], form_idx, T)) {
      fprintf(stderr, "FAIL: the identity is not a scalar on the space of the form\n");
      is_ok = false;
    }
    // the entry is in the row of an entry of the vector which is not 0
    j = 0;
    while (nf_elem_is_zero(evs[0]->eigenvecs[form_idx][j], evs[0]->nfs[form_idx]))
      j++;
    fmpq_mat_zero(T);
    fmpq_set_si(fmpq_mat_entry(T, j, 0), 1, 1);
    if (ev_is_scalar_on_old_space(evs[0], form_idx, T)) {
      fprintf(stderr, "FAIL: a matrix with one entry is a scalar on the space of the form\n");
      is_ok = false;
    }
    fmpq_mat_clear(T);
  }

  for (c = 0; c < genus->num_conductors; c++)
    eigenvalues_clear(evs[c]);
  free(evs);
  genus_clear(genus);

  return (is_ok ? SUCCESS : FAIL);
}

// The lattice of discriminant 124 = 2^2 * 31 has a form of conductor 1 for
// which no eigenvector is found, and whose field is quadratic. The operators
// T_3, at a prime of the decomposition, and T_31 are scalars on its space, and
// T_2 and T_{2^2} are not. The second matrix is that of the same lattice in
// another basis.
STATUS test_old_space_scalars()
{
  Z64 Q_124[15] = {1,0,1,0,-1,1,1,0,0,2,0,-1,2,-1,3};
  Z64 Q_124_other[15] = {4,-2,0,-6,0,1,1,2,0,2,1,-1,3,-1,3};
  int ps[4] = {3, 31, 2, 2};
  int ks[4] = {1, 1, 1, 2};
  bool is_scalar[4] = {true, true, false, false};
  STATUS ret;

  ret = test_old_space_scalars_of(Q_124, 4, ps, ks, is_scalar);
  ret |= test_old_space_scalars_of(Q_124_other, 4, ps, ks, is_scalar);

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

// The eigenvalues of T_p, and also those of T_{p^2} if square is true.
// The prime p may divide the discriminant.
STATUS compute_eigenvalues(const genus_t genus, int p, bool square, bool nonlifts, int num_idxs, const int* idxs)
{
  example_evs_t dummy;
  STATUS ret;
  int num_ps[2] = {1, square ? 1 : 0};
  const int* ps[2] = {
    (int []){p},
    (int []){p}
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

// The matrices of T_{p^2}, for the conductor 1 if c = 1 is given and for all
// the conductors otherwise. The prime p may divide the discriminant.
STATUS compute_hecke_square_matrix(const genus_t genus, slong p, bool has_c, slong c)
{
  matrix_TYP** hecke;
  slong c_idx, num_printed;

  clock_t cpuclock_0, cpuclock_1;
  double cputime, cpudiff;

  // as for T_p, only the conductor 1 is printed on its own
  num_printed = genus->num_conductors;
  if (has_c && (conductor_index(genus, c) == 0))
    num_printed = 1;

  cpuclock_0 = clock();

  hecke = hecke_matrices_nbr_data_all_conductors(genus, p, 2);

  cpuclock_1 = clock();
  cpudiff = cpuclock_1 - cpuclock_0;
  cputime = cpudiff / CLOCKS_PER_SEC;

  printf("{");
  for (c_idx = 0; c_idx < num_printed; c_idx++) {
    printf("%ld : ", genus->conductors[c_idx]);
    print_mat_dense(hecke[c_idx]);
    if (c_idx != num_printed-1)
      printf(",");
  }
  printf("}");

  fprintf(stderr, "computing hecke matrices took %f sec\n", cputime);

  for (c_idx = 0; c_idx < genus->num_conductors; c_idx++)
    free_mat(hecke[c_idx]);

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
