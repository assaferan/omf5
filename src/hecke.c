#include <assert.h>

#include <carat/matrix.h>

#include <flint/fmpz_mat.h>

#include "nbr_data.h"
#include "neighbor.h"

#include "arith.h"
#include "hecke.h"
#include "matrix_tools.h"
#include "spinor.h"

#include "typedefs.h"

int process_isotropic_vector_nbr_data(nbr_data_t nbr_man, int* T, const hash_table* genus, slong spinor_prime,
				      double* theta_time, double* isom_time, double* total_time, int* num_isom)
{
  int i;
  clock_t cputime;
  matrix_TYP* nbr;

  fmpz_mat_t nbr_fmpz, nbr_isom;
  fmpz_mat_init(nbr_fmpz, N, N);
  fmpz_mat_init(nbr_isom, N, N);
  
  nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
  matrix_TYP_init_set_fmpz_mat(&nbr, nbr_fmpz);

  cputime = clock();
  i = indexof(genus, nbr, 0, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

#ifdef DEBUG
  if ((i < 0) || (i > genus->num_stored)) {
    printf("Error! Couldn't find element in genus!\n");
    exit(-1);
  }
#endif // DEBUG

  if (spinor_prime == 0) {
    T[i]++;
  }
  else {
    T[i] += fmpz_mat_spinor_norm(nbr_man->q, nbr_isom, fq_nmod_ctx_prime(nbr_man->GF), spinor_prime);
  }

  fmpz_mat_clear(nbr_fmpz);
  fmpz_mat_clear(nbr_isom);

  free_mat(nbr);
  
  return 0;
}
  
int process_isotropic_vector(neighbor_manager* nbr_man, int* T, const hash_table* genus,
			     double* theta_time, double* isom_time, double* total_time, int* num_isom)

{

  int i;
  clock_t cputime;
  matrix_TYP* nbr;

  nbr = build_nb(nbr_man);

  cputime = clock();
  i = indexof(genus, nbr, 0, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

#ifdef DEBUG
  if ((i < 0) || (i > genus->num_stored)) {
    printf("Error! Couldn't find element in genus!\n");
    exit(-1);
  }
#endif // DEBUG
  
  T[i]++;

  free_mat(nbr);
  
  return 0;
}

// !! TODO - use i to cut the parameters to chunks, need to convert from a number to a vector in the parameter space
int process_neighbour_chunk(int* T, int p, int i, int gen_idx, const hash_table* genus,
			    double* theta_time, double* isom_time, double* total_time, int* num_isom)
{
  matrix_TYP *Q;
  neighbor_manager nbr_man;
  int lc;

  lc = 0;

  Q = genus->keys[gen_idx];

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  print_mat(Q);
#endif // DEBUG_LEVEL_FULL

  init_nbr_process(&nbr_man, Q, p, i);

  while (!(has_ended(&nbr_man))) {
    process_isotropic_vector(&nbr_man, T, genus, theta_time, isom_time, total_time, num_isom);
    advance_nbr_process(&nbr_man);
    lc++;
  }

  free_nbr_process(&nbr_man);
 
  return lc;
}

int process_neighbour_chunk_nbr_data(int* T, int p, int k, int gen_idx, const hash_table* genus, slong spinor_prime,
				     double* theta_time, double* isom_time, double* total_time, int* num_isom)
{
  matrix_TYP *Q;
  nbr_data_t nbr_man;
  int lc;

  lc = 0;

  Q = genus->keys[gen_idx];

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  print_mat(Q);
#endif // DEBUG_LEVEL_FULL

  nbr_data_init(nbr_man, Q, p, k);

  while (!(nbr_data_has_ended(nbr_man))) {
    process_isotropic_vector_nbr_data(nbr_man, T, genus, spinor_prime, theta_time, isom_time, total_time, num_isom);
    nbr_data_get_next_neighbor(nbr_man);
    lc++;
  }

  nbr_data_clear(nbr_man);
 
  return lc;
}

// assumes T is initialized to zeros

void hecke_col_nbr_data(int* T, int p, int k, int gen_idx, const hash_table* genus, slong spinor_prime)
{
  int num_isom, lc;
  double theta_time, isom_time, total_time;
  
  num_isom = lc = 0;
  theta_time = isom_time = total_time = 0;

  lc += process_neighbour_chunk_nbr_data(T, p, k, gen_idx, genus, spinor_prime,
					 &theta_time, &isom_time, &total_time, &num_isom);

#ifdef DEBUG
  printf("theta_time = %f, isom_time = %f, total_time = %f, num_isom = %d / %d \n", theta_time/lc, isom_time/lc, total_time, num_isom, lc);
#endif // DEBUG

  return;
}

void hecke_col(int* T, int p, int gen_idx, const hash_table* genus)
{
  int num;
  int num_isom, lc;
  double theta_time, isom_time, total_time;
  num_isom = lc = 0;
  theta_time = isom_time = total_time = 0;

  for (num = 0; num < p; num++) {
    lc += process_neighbour_chunk(T, p, num, gen_idx, genus, &theta_time, &isom_time, &total_time, &num_isom);
  }

#ifdef DEBUG
  printf("theta_time = %f, isom_time = %f, total_time = %f, num_isom = %d / %d \n", theta_time/lc, isom_time/lc, total_time, num_isom, lc);
#endif // DEBUG

  return;
}

matrix_TYP* hecke_matrix(const hash_table* genus, int p)
{
  matrix_TYP* hecke;
  int gen_idx;

  hecke = init_mat(genus->num_stored, genus->num_stored, "");

  for (gen_idx = 0; gen_idx < genus->num_stored; gen_idx++)
    hecke_col(hecke->array.SZ[gen_idx], p, gen_idx, genus);
  
  return hecke;
}

void get_hecke_ev_nbr_data(nf_elem_t e, const hash_table* genus, eigenvalues* evs,
			   int p, int k, int ev_idx, slong spinor_prime)
{
  int* a;
  int num, i, pivot;
  nf_elem_t prod;
  fmpz_mat_t q;
  fmpz_t disc;

  nf_elem_init(e, evs->nfs[ev_idx]);
  nf_elem_init(prod, evs->nfs[ev_idx]);
  a = (int*)malloc(genus->num_stored * sizeof(int));
  for (num = 0; num < genus->num_stored; num++)
    a[num] = 0;

  for (pivot = 0; pivot < evs->dim;) {
    if (!(nf_elem_is_zero(evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]))) {
      break;
    }
    pivot++;
  }
#ifdef DEBUG
  clock_t cpuclock;
  double cputime;
  
  cpuclock = clock();
#endif // DEBUG
  
  hecke_col_nbr_data(a, p, k, pivot, genus, spinor_prime);

#ifdef DEBUG
  cpuclock = clock() - cpuclock;
  cputime = cpuclock / CLOCKS_PER_SEC;
#endif // DEBUG
  
  nf_elem_zero(e, evs->nfs[ev_idx]);
  for (i = 0; i < evs->dim; i++) {
    nf_elem_set_si(prod, a[i], evs->nfs[ev_idx]);
    nf_elem_mul(prod, prod, evs->eigenvecs[ev_idx][i], evs->nfs[ev_idx]);
    nf_elem_add(e, e, prod, evs->nfs[ev_idx]);
  }
  nf_elem_div(e, e, evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]);

  // handling the case p divides the discriminant
  if (genus->num_stored > 0) {
    fmpz_mat_init_set_matrix_TYP(q, genus->keys[0]);
    fmpz_init(disc);
    fmpz_mat_det(disc, q);
    fmpz_divexact_si(disc,disc,2);
    if (fmpz_get_si(disc) % p == 0) {
      nf_elem_add_si(e, e, 1, evs->nfs[ev_idx]);
    }
    fmpz_mat_clear(q);
    fmpz_clear(disc);
  }

#ifdef DEBUG
  printf("%4d ", p);
  nf_elem_print_pretty(e, evs->nfs[ev_idx], "a");
  printf(" - ");
  for (num = 0; num < genus->num_stored; num++)
    printf("%10d ", a[num]);
  
  printf("- %10f\n", cputime);
#endif // DEBUG

  nf_elem_clear(prod, evs->nfs[ev_idx]);
  free(a);
  
  return;
}

void get_hecke_ev(nf_elem_t e, const hash_table* genus, eigenvalues* evs, int p, int ev_idx)
{
  int* a;
  int num, i, pivot;
  nf_elem_t prod;

  nf_elem_init(e, evs->nfs[ev_idx]);
  nf_elem_init(prod, evs->nfs[ev_idx]);
  a = (int*)malloc(genus->num_stored * sizeof(int));
  for (num = 0; num < genus->num_stored; num++)
    a[num] = 0;

  for (pivot = 0; pivot < evs->dim;) {
    if (!(nf_elem_is_zero(evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]))) {
      break;
    }
    pivot++;
  }
  clock_t cpuclock;
  double cputime;
  
  cpuclock = clock();
  
  hecke_col(a, p, pivot, genus);

  cpuclock = clock() - cpuclock;
  cputime = cpuclock / CLOCKS_PER_SEC;

  nf_elem_zero(e, evs->nfs[ev_idx]);
  for (i = 0; i < evs->dim; i++) {
    nf_elem_set_si(prod, a[i], evs->nfs[ev_idx]);
    nf_elem_mul(prod, prod, evs->eigenvecs[ev_idx][i], evs->nfs[ev_idx]);
    nf_elem_add(e, e, prod, evs->nfs[ev_idx]);
  }
  nf_elem_div(e, e, evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]);

#ifdef DEBUG
  printf("%4d ", p);
  nf_elem_print_pretty(e, evs->nfs[ev_idx], "a");
  printf(" - ");
  for (num = 0; num < genus->num_stored; num++)
    printf("%10d ", a[num]);
  
  printf("- %10f\n", cputime);
#endif // DEBUG

  nf_elem_clear(prod, evs->nfs[ev_idx]);
  free(a);
  
  return;
}

void get_hecke_fmpq_mat(fmpq_mat_t hecke_fmpq_mat, const hash_table* genus, int p)
{
  matrix_TYP* hecke_mat;
   
  hecke_mat = hecke_matrix(genus, p);
  fmpq_mat_init_set_matrix_TYP(hecke_fmpq_mat, hecke_mat);
  // we transpose to match with our conventions
  fmpq_mat_transpose(hecke_fmpq_mat, hecke_fmpq_mat);

  free_mat(hecke_mat);
  return;
}

// TODO - get rid of the unnecessary recursion
bool decomposition_finite_subspace(decomposition* decomp, const hash_table* genus, const fmpq_mat_t basis_V,
				   const int* ps, slong idx, slong num_ps)
{
  slong i, dim_V, next_idx;
  fmpq_mat_t T, fT, W;
  fmpq_poly_t f;
  fmpq_poly_factor_t fac;
  bool is_complete, is_complete_W;
  fmpq_mat_t hecke_fmpq_mat;

  dim_V = fmpq_mat_nrows(basis_V);
  if (dim_V == 0) {
    return true;
  }

  if (idx >= num_ps) {
    (decomp->num)++;
    decomp->bases = (fmpq_mat_t*)realloc(decomp->bases, (decomp->num)*sizeof(fmpq_mat_t));
    fmpq_mat_init(decomp->bases[decomp->num-1], dim_V, dim_V);
    fmpq_mat_one(decomp->bases[decomp->num-1]);
    return false;
  }

  fmpq_mat_init(T, dim_V, dim_V);
  fmpq_mat_init(fT, dim_V, dim_V);
  fmpq_poly_init(f);

  get_hecke_fmpq_mat(hecke_fmpq_mat, genus, ps[idx]);
  restrict_mat(T, hecke_fmpq_mat, basis_V);
  fmpq_mat_charpoly(f, T);
  fmpq_poly_factor(fac, f);

  is_complete = true;
  for (i = 0; i < fac->num; i++) {
    fmpq_poly_evaluate_fmpq_mat(fT, &(fac->p[i]), T);
    kernel_on(W, fT, basis_V);
    assert(W->rows != 0);
    // optimally, we will already compute the eigenvectors at this point.
    if (fac->exp[i] == 1) { // test for irreduciblity
      (decomp->num)++;
      decomp->bases = (fmpq_mat_t*)realloc(decomp->bases, (decomp->num)*sizeof(fmpq_mat_t));
      fmpq_mat_init_set(decomp->bases[decomp->num-1], W);
      is_complete_W = true;
    }
    else {
      next_idx = (fmpq_mat_nrows(W) == dim_V) ? idx+1 : 0;
      is_complete_W = decomposition_finite_subspace(decomp, genus, W, ps, next_idx, num_ps);
    }
    is_complete = (is_complete) && (is_complete_W);
  }

  fmpq_mat_clear(hecke_fmpq_mat);
  fmpq_mat_clear(W);
  fmpq_mat_clear(fT);
  fmpq_mat_clear(T);
  fmpq_poly_factor_free(fac);
  fmpq_poly_clear(f);
  
  return is_complete;
}

bool decomposition_finite(decomposition* decomp, const hash_table* genus, const int* ps, slong num_ps)
{
  fmpq_mat_t basis_M;
  slong dim;
  bool is_complete;
  
  decomp->num = 0;
  decomp->bases = NULL;

  if (num_ps == 0)
    dim = 0;
  else
    dim = genus->num_stored;

  fmpq_mat_init(basis_M, dim, dim);
  fmpq_mat_one(basis_M);

  is_complete = decomposition_finite_subspace(decomp, genus, basis_M, ps, 0, num_ps);
  
  fmpq_mat_clear(basis_M);
  
  return is_complete;
}

decomposition* decompose(const hash_table* genus)
{
  slong bound, num_ps;
  int* ps;
  decomposition* decomp;
  bool is_complete;

  is_complete = false;
  decomp = (decomposition*)malloc(sizeof(decomposition));
  bound = 10;

  while (!(is_complete)) {
    num_ps = primes_up_to(&ps, bound);
    is_complete = decomposition_finite(decomp, genus, ps, num_ps);
    free(ps);
    bound *= 2;
  }

  return decomp;
}

eigenvalues* hecke_eigenforms(const hash_table* genus)
{
  fmpq_mat_t T;
  slong i;
  eigenvalues* evs;
  decomposition* D;
  
  evs = (eigenvalues*)malloc(sizeof(eigenvalues));
  D = decompose(genus);

  eigenvalues_init(&evs, D->num, genus->num_stored);

  get_hecke_fmpq_mat(T, genus, 2);
  
  for (i = 0; i < D->num; i++) {
    get_eigenvector(evs->eigenvecs[i], evs->nfs[i], T, D->bases[i]);
    nf_elem_init(evs->eigenvals[i], evs->nfs[i]);
    nf_elem_gen(evs->eigenvals[i], evs->nfs[i]);
  }
  
  return evs;
}
