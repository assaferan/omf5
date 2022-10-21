#include <assert.h>

#include <carat/matrix.h>

#include <flint/fmpz_mat.h>

#include "nbr_data.h"
#include "neighbor.h"

#include "arith.h"
#include "decomposition.h"
#include "hecke.h"
#include "matrix_tools.h"

#include "typedefs.h"

#define likely(x)   __builtin_expect(!!(x), 1)
#define unlikely(x) __builtin_expect(!!(x), 0)

int process_isotropic_vector_nbr_data(nbr_data_t nbr_man, int* T, const genus_t genus,
				      double* theta_time, double* isom_time, double* total_time, int* num_isom)
{
  int i;
  clock_t cputime;
  matrix_TYP* nbr;
  fmpz_mat_t nbr_fmpz, nbr_isom;
  slong n;

  assert(genus->genus_reps->num_stored > 0);
  n = genus->genus_reps->keys[0]->rows;
  assert(n == genus->genus_reps->keys[0]->cols);
  
  fmpz_mat_init(nbr_fmpz, n, n);
  fmpz_mat_init(nbr_isom, n, n);
  
  nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
  matrix_TYP_init_set_fmpz_mat(&nbr, nbr_fmpz);

  cputime = clock();
  i = hash_table_indexof(genus->genus_reps, nbr, 0, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

#ifdef DEBUG
  if ((i < 0) || (i > genus->genus_reps->num_stored)) {
    printf("Error! Couldn't find element in genus!\n");
    exit(-1);
  }
#endif // DEBUG

  // TODO - add here the spinor norm
  T[i]++;

  fmpz_mat_clear(nbr_fmpz);
  fmpz_mat_clear(nbr_isom);

  free_mat(nbr);
  
  return 0;
}

int process_isotropic_vector_nbr_data_all_conductors(nbr_data_t nbr_man, W64* spin_vals,
						     const genus_t genus,
						     double* theta_time, double* isom_time,
						     double* total_time, int* num_isom, int gen_idx)
{
  int i;
  clock_t cputime;
  matrix_TYP* nbr;
  isometry_t s_nbr, s_inv;
  matrix_TYP* hash_isom;
  fmpz_mat_t nbr_fmpz, nbr_isom;
  slong n;

  assert(genus->genus_reps->num_stored > 0);
  n = genus->genus_reps->keys[0]->rows;
  assert(n == genus->genus_reps->keys[0]->cols);
  
  fmpz_mat_init(nbr_fmpz, n, n);
  fmpz_mat_init(nbr_isom, n, n);
  
  nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
  matrix_TYP_init_set_fmpz_mat(&nbr, nbr_fmpz);
  isometry_init_fmpz_mat(s_nbr, nbr_isom, fmpz_get_si(fq_nmod_ctx_prime(nbr_man->GF)));
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
  
  cputime = clock();
  i = hash_table_index_and_isom(genus->genus_reps, nbr, &hash_isom, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

#ifdef DEBUG
  if ((i < 0) || (i > genus->genus_reps->num_stored)) {
    printf("Error! Couldn't find element in genus!\n");
    exit(-1);
  }
#endif // DEBUG

  // TODO - compute only the upper half, and complete using hermitian property
  // Determine which subspaces this representative contributes.

  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
  assert(isometry_is_isom(genus->isoms[gen_idx], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[gen_idx]));
  isometry_mul(s_nbr, genus->isoms[gen_idx], s_nbr);
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], nbr));
  assert(is_isometry(hash_isom, nbr, genus->genus_reps->keys[i], 1));
  isometry_mul_mat_right(s_nbr, s_nbr, hash_isom);
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], genus->genus_reps->keys[i]));
  assert(isometry_is_isom(genus->isoms[i], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[i]));
  isometry_inv(s_inv, genus->isoms[i]);
  assert(isometry_is_isom(s_inv, genus->genus_reps->keys[i],
			  genus->genus_reps->keys[0]));
  isometry_mul(s_nbr, s_nbr, s_inv);
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0],
			  genus->genus_reps->keys[0]));
  
  *spin_vals = ((i << (genus->spinor->num_primes)) | spinor_norm_isom(genus->spinor, s_nbr));

  isometry_clear(s_nbr);
  isometry_clear(s_inv);
  
  fmpz_mat_clear(nbr_fmpz);
  fmpz_mat_clear(nbr_isom);

  free_mat(nbr);
  
  return 0;
}
  
int process_isotropic_vector(neighbor_manager* nbr_man, int* T, const genus_t genus,
			     double* theta_time, double* isom_time, double* total_time, int* num_isom)

{

  int i;
  clock_t cputime;
  matrix_TYP* nbr;

  nbr = build_nb(nbr_man);

  cputime = clock();
  i = hash_table_indexof(genus->genus_reps, nbr, 0, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

#ifdef DEBUG
  if ((i < 0) || (i > genus->genus_reps->num_stored)) {
    printf("Error! Couldn't find element in genus!\n");
    exit(-1);
  }
#endif // DEBUG
  
  T[i]++;

  free_mat(nbr);
  
  return 0;
}



// !! TODO - use i to cut the parameters to chunks, need to convert from a number to a vector in the parameter space
int process_neighbour_chunk(int* T, int p, int i, int gen_idx, const genus_t genus,
			    double* theta_time, double* isom_time, double* total_time, int* num_isom)
{
  matrix_TYP *Q;
  neighbor_manager nbr_man;
  int lc;
  
  lc = 0;

  Q = genus->genus_reps->keys[gen_idx];

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

int process_neighbour_chunk_nbr_data_all_conductors(W64* spin_vals, int p, int k, int gen_idx,
						    const genus_t genus,
						    double* theta_time, double* isom_time,
						    double* total_time,
						    int* num_isom)
{
  matrix_TYP *Q;
  nbr_data_t nbr_man;
  int lc;

  lc = 0;

  Q = genus->genus_reps->keys[gen_idx];

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  print_mat(Q);
#endif // DEBUG_LEVEL_FULL

  assert(k > 0);
  
  nbr_data_init(nbr_man, Q, p, k);

  while (!(nbr_data_has_ended(nbr_man))) {
    process_isotropic_vector_nbr_data_all_conductors(nbr_man, &(spin_vals[lc]), genus,
						     theta_time, isom_time, total_time,
						     num_isom, gen_idx);
    nbr_data_get_next_neighbor(nbr_man);
    lc++;
  }

  nbr_data_clear(nbr_man);
 
  return lc;
}

int process_neighbour_chunk_nbr_data(int* T, int p, int k, int gen_idx, const genus_t genus,
				     double* theta_time, double* isom_time, double* total_time, int* num_isom)
{
  matrix_TYP *Q;
  nbr_data_t nbr_man;
  int lc;

  lc = 0;

  Q = genus->genus_reps->keys[gen_idx];

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  print_mat(Q);
#endif // DEBUG_LEVEL_FULL

  nbr_data_init(nbr_man, Q, p, k);

  while (!(nbr_data_has_ended(nbr_man))) {
    process_isotropic_vector_nbr_data(nbr_man, T, genus, theta_time, isom_time, total_time, num_isom);
    nbr_data_get_next_neighbor(nbr_man);
    lc++;
  }

  nbr_data_clear(nbr_man);
 
  return lc;
}

// assumes T is initialized to zeros

void hecke_col_nbr_data(int* T, int p, int k, int gen_idx, const genus_t genus)
{
  int num_isom, lc;
  double theta_time, isom_time, total_time;
  
  num_isom = lc = 0;
  theta_time = isom_time = total_time = 0;

  lc += process_neighbour_chunk_nbr_data(T, p, k, gen_idx, genus, 
					 &theta_time, &isom_time, &total_time, &num_isom);

#ifdef DEBUG
  printf("theta_time = %f, isom_time = %f, total_time = %f, num_isom = %d / %d \n", theta_time/lc, isom_time/lc, total_time, num_isom, lc);
#endif // DEBUG

  return;
}

void hecke_col(int* T, int p, int gen_idx, const genus_t genus)
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

void hecke_col_nbr_data_all_conductors(W64* spin_vals, int p, int k, int gen_idx, const genus_t genus)
{
  int num_isom, lc;
  double theta_time, isom_time, total_time;
  
  num_isom = lc = 0;
  theta_time = isom_time = total_time = 0;

  lc += process_neighbour_chunk_nbr_data_all_conductors(&(spin_vals[lc]), p, k, gen_idx, genus,
							&theta_time, &isom_time, &total_time, &num_isom);

#ifdef DEBUG
  printf("theta_time = %f, isom_time = %f, total_time = %f, num_isom = %d / %d \n", theta_time/lc, isom_time/lc, total_time, num_isom, lc);
#endif // DEBUG

  return;
}

int char_vals[256] = {
  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,
  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,
  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,
  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1,
  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1,
  -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,
  1, -1, -1,  1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,
  1, -1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1,  1, -1, -1,  1,
  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,  1, -1, -1,  1,  1, -1,
  -1,  1, -1,  1,  1, -1,  1, -1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1,
  1, -1, -1,  1
};
 
int char_val(W64 x)
{
  if (x <= 0xff) return char_vals[x];
  int value = 1;
  while (x) {
    value *= char_vals[x&0xff];
    x >>= 8;
  }
  return value;
}

matrix_TYP** hecke_matrices_all_conductors(const genus_t genus, int p, int k)
{
  matrix_TYP** hecke;
  int gen_idx, rpos;
  slong c, npos, num_nbrs, spin_idx;
  int* row;
  slong* lut;
  W64* spin_vals;
  W64 r;
  int* hecke_ptr;

  hecke = (matrix_TYP**)malloc(genus->num_conductors * sizeof(matrix_TYP*));
  hecke_ptr = (int*)malloc(genus->num_conductors * sizeof(int));
  // hecke_ptr = (slong**)malloc(genus->num_conductors * sizeof(slong*));
  for (c = 0; c < genus -> num_conductors; c++) {
    // hecke_ptr[c] = (slong*)malloc(genus->dims[c] * genus->dims[c] * sizeof(slong));
    hecke[c] = init_mat(genus->dims[c], genus->dims[c], "");
    hecke_ptr[c] = 0;
  }

  num_nbrs = p*p*p+p*p+p+1;
  if (k == 2) {
    num_nbrs *= p;
  }
  assert((k == 1) || (k == 2));
  spin_vals = (W64*)malloc(num_nbrs*sizeof(W64));
  
  for (gen_idx = 0; gen_idx < genus->genus_reps->num_stored; gen_idx++) {
    hecke_col_nbr_data_all_conductors(spin_vals, p, k, gen_idx, genus);
    for (c = 0; c < genus->num_conductors; c++) {
      lut = genus->lut_positions[c];
      npos = lut[gen_idx];
      if (unlikely(npos==-1)) continue;
      row = hecke[c]->array.SZ[hecke_ptr[c]];
      for (spin_idx = 0; spin_idx < num_nbrs; spin_idx++) {
	r = spin_vals[spin_idx] >> (genus->spinor->num_primes);
	rpos = lut[r];
	if (unlikely(rpos==-1)) continue;
	row[rpos] += char_val(spin_vals[spin_idx] &c);
      }
      hecke_ptr[c]++;
    }
    
  }
  free(hecke_ptr);
  free(spin_vals);
  return hecke;
}

matrix_TYP* hecke_matrix(const genus_t genus, int p)
{
  matrix_TYP* hecke;
  int gen_idx;

  hecke = init_mat(genus->genus_reps->num_stored, genus->genus_reps->num_stored, "");

  for (gen_idx = 0; gen_idx < genus->genus_reps->num_stored; gen_idx++)
    hecke_col(hecke->array.SZ[gen_idx], p, gen_idx, genus);
  
  return hecke;
}

void get_hecke_ev_nbr_data(nf_elem_t e, const genus_t genus, const eigenvalues_t evs,
			   int p, int k, int ev_idx)
{
  int* a;
  int num, i, pivot;
  nf_elem_t prod;
  fmpz_mat_t q;
  fmpz_t disc;

  nf_elem_init(e, evs->nfs[ev_idx]);
  nf_elem_init(prod, evs->nfs[ev_idx]);
  a = (int*)malloc(evs->dim * sizeof(int));
  for (num = 0; num < evs->dim; num++)
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
  
  hecke_col_nbr_data(a, p, k, pivot, genus);

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
  if (genus->genus_reps->num_stored > 0) {
    fmpz_mat_init_set_matrix_TYP(q, genus->genus_reps->keys[0]);
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
  for (num = 0; num < genus->genus_reps->num_stored; num++)
    printf("%10d ", a[num]);
  
  printf("- %10f\n", cputime);
#endif // DEBUG

  nf_elem_clear(prod, evs->nfs[ev_idx]);
  free(a);
  
  return;
}


void get_hecke_ev_nbr_data_all_conductors(nf_elem_t e, const genus_t genus,
					  const eigenvalues_t evs,
					  int p, int k, int ev_idx, slong ev_cond)
{
  int** hecke;
  int num, i, pivot, num_nbrs,spin_idx;
  slong c, npos, rpos, gen_idx;
  nf_elem_t prod;
  fmpz_mat_t q;
  fmpz_t disc;
  slong* lut;
  int* row;
  W64* spin_vals;
  W64 r;
#ifdef DEBUG
  clock_t cpuclock;
  double cputime;
  
  cpuclock = clock();
#endif // DEBUG
  
  nf_elem_init(e, evs->nfs[ev_idx]);
  nf_elem_init(prod, evs->nfs[ev_idx]);

  for (pivot = 0; pivot < evs->dim;) {
    if (!(nf_elem_is_zero(evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]))) {
      break;
    }
    pivot++;
  }

  for (gen_idx = 0; gen_idx < genus->dims[0];) {
    if (pivot == genus->lut_positions[ev_cond][gen_idx])
      break;
    gen_idx++;
  }
  
  hecke = (int**)malloc(genus->num_conductors * sizeof(int*));
  
  for (c = 0; c < genus->num_conductors; c++) {
    hecke[c] = (int*)malloc(genus->dims[c] * sizeof(int));
    for (i = 0; i < genus->dims[c]; i++)
      hecke[c][i] = 0;
  }
  
  num_nbrs = p*p*p+p*p+p+1;
  if (k == 2) {
    num_nbrs *= p;
  }
  assert((k == 1) || (k == 2));
  spin_vals = (W64*)malloc(num_nbrs*sizeof(W64));
  hecke_col_nbr_data_all_conductors(spin_vals, p, k, gen_idx, genus);

  for (c = 0; c < genus->num_conductors; c++) {
    lut = genus->lut_positions[c];
    npos = lut[gen_idx];
    if (unlikely(npos==-1)) continue;
    row = hecke[c];
    for (spin_idx = 0; spin_idx < num_nbrs; spin_idx++) {
      r = spin_vals[spin_idx] >> (genus->spinor->num_primes);
      rpos = lut[r];
      if (unlikely(rpos==-1)) continue;
      row[rpos] += char_val(spin_vals[spin_idx] & c);
    }
  }

  free(spin_vals);

#ifdef DEBUG
  cpuclock = clock() - cpuclock;
  cputime = cpuclock / CLOCKS_PER_SEC;
#endif // DEBUG
  
  nf_elem_zero(e, evs->nfs[ev_idx]);
  for (i = 0; i < evs->dim; i++) {
    nf_elem_set_si(prod, hecke[ev_cond][i], evs->nfs[ev_idx]);
    nf_elem_mul(prod, prod, evs->eigenvecs[ev_idx][i], evs->nfs[ev_idx]);
    nf_elem_add(e, e, prod, evs->nfs[ev_idx]);
  }
  nf_elem_div(e, e, evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]);

  // handling the case p divides the discriminant
  if (genus->genus_reps->num_stored > 0) {
    fmpz_mat_init_set_matrix_TYP(q, genus->genus_reps->keys[0]);
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
  for (num = 0; num < genus->genus_reps->num_stored; num++)
    printf("%10d ", hecke[ev_cond][num]);
  
  printf("- %10f\n", cputime);
#endif // DEBUG

  nf_elem_clear(prod, evs->nfs[ev_idx]);

  for (c = 0; c < genus->num_conductors; c++)
    free(hecke[c]);

  free(hecke);
  
  return;
}


void get_hecke_ev(nf_elem_t e, const genus_t genus, const eigenvalues_t evs, int p, int ev_idx)
{
  int* a;
  int num, i, pivot;
  nf_elem_t prod;

  nf_elem_init(e, evs->nfs[ev_idx]);
  nf_elem_init(prod, evs->nfs[ev_idx]);
  a = (int*)malloc(genus->genus_reps->num_stored * sizeof(int));
  for (num = 0; num < genus->genus_reps->num_stored; num++)
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
  
  hecke_col(a, p, pivot, genus);

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

#ifdef DEBUG
  printf("%4d ", p);
  nf_elem_print_pretty(e, evs->nfs[ev_idx], "a");
  printf(" - ");
  for (num = 0; num < genus->genus_reps->num_stored; num++)
    printf("%10d ", a[num]);
  
  printf("- %10f\n", cputime);
#endif // DEBUG

  nf_elem_clear(prod, evs->nfs[ev_idx]);
  free(a);
  
  return;
}

void get_hecke_fmpq_mat(fmpq_mat_t hecke_fmpq_mat, const genus_t genus, int p)
{
  matrix_TYP* hecke_mat;
   
  hecke_mat = hecke_matrix(genus, p);
  fmpq_mat_init_set_matrix_TYP(hecke_fmpq_mat, hecke_mat);
  // we transpose to match with our conventions
  fmpq_mat_transpose(hecke_fmpq_mat, hecke_fmpq_mat);

  free_mat(hecke_mat);
  return;
}
 
void get_hecke_fmpq_mat_all_conductors(fmpq_mat_t* hecke_fmpq_mat, const genus_t genus, int p, int k)
{
  matrix_TYP** hecke_matrices;
  slong c;

  assert(hecke_fmpq_mat != NULL);
  hecke_matrices = hecke_matrices_all_conductors(genus, p, k);
  for (c = 0; c < genus->num_conductors; c++) {
    fmpq_mat_init_set_matrix_TYP(hecke_fmpq_mat[c], hecke_matrices[c]);
    // we transpose to match with our conventions
    fmpq_mat_transpose(hecke_fmpq_mat[c], hecke_fmpq_mat[c]);
    free_mat(hecke_matrices[c]);
  }
  
  return;
}

void hecke_eigenforms(eigenvalues_t evs, const decomposition_t D, const genus_t genus, slong c)
{
  slong i, p_idx;
  bool ev_init;

  eigenvalues_init(evs, D->num[c], genus->dims[c]);
  
  for (i = 0; i < D->num[c]; i++) {
    p_idx = 0;
    ev_init = false;
    do {
      ev_init = get_eigenvector_on_subspace(evs->eigenvecs[i], evs->nfs[i],
					    D->hecke[p_idx][c], D->bases[c][i]);
      p_idx++;
    } while (!ev_init);
    nf_elem_init(evs->eigenvals[i], evs->nfs[i]);
    nf_elem_gen(evs->eigenvals[i], evs->nfs[i]);
  }
  
  return;
}

eigenvalues_t* hecke_eigenforms_all_conductors(const genus_t genus)
{
  eigenvalues_t* all_evs;
  decomposition_t D;
  slong c;

  decomposition_init(D, genus->num_conductors);
 
  all_evs = (eigenvalues_t*)malloc(genus->num_conductors * sizeof(eigenvalues_t));
  for (c = 0; c < genus->num_conductors; c++) {
#ifdef DEBUG
    printf("decomposing conductor %ld\n", genus->conductors[c]);
#endif // DEBUG
    decompose(D, genus, c);
#ifdef DEBUG
    printf("computing eigenvectors\n");
#endif // DEBUG
    hecke_eigenforms(all_evs[c], D, genus, c);
  }

  decomposition_clear(D);
  return all_evs;
}
