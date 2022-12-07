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
  square_matrix_t nbr;
  fmpz_mat_t nbr_fmpz, nbr_isom;

  assert(genus->genus_reps->num_stored > 0);
  
  fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
  fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
  
  nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
  square_matrix_set_fmpz_mat(nbr, nbr_fmpz);

  cputime = clock();
  i = hash_table_indexof(genus->genus_reps, nbr, 0, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

  assert( (i >= 0) && (i <= genus->dims[0]));

  T[i]++;

  fmpz_mat_clear(nbr_fmpz);
  fmpz_mat_clear(nbr_isom);

  square_matrix_clear(nbr);
  
  return 0;
}

int process_isotropic_vector_nbr_data_all_conductors(nbr_data_t nbr_man, W64* spin_vals,
						     const genus_t genus,
						     double* theta_time, double* isom_time,
						     double* total_time, int* num_isom, int gen_idx)
{
  int i;
  clock_t cputime;
  square_matrix_t nbr;
  isometry_t s_nbr, s_inv;
  isometry_t hash_isom;
  fmpz_mat_t nbr_fmpz, nbr_isom;

  assert(genus->genus_reps->num_stored > 0);
  
  fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
  fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
  
  nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
  square_matrix_set_fmpz_mat(nbr, nbr_fmpz);
  isometry_init_set_fmpz_mat(s_nbr, nbr_isom, fmpz_get_si(fq_nmod_ctx_prime(nbr_man->GF)));
  assert(gen_idx < genus->genus_reps->num_stored);
#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
#endif // DEBUG_LEVEL_FULL
  
  cputime = clock();
  i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

  assert( (0 <= i) && (i < genus->genus_reps->num_stored) );

  // TODO - compute only the upper half, and complete using hermitian property
  // Determine which subspaces this representative contributes.

  assert(gen_idx < genus->genus_reps->num_stored);
#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
  assert(isometry_is_isom(genus->isoms[gen_idx], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[gen_idx]));
#endif // DEBUG_LEVEL_FULL
  isometry_muleq_left(s_nbr, genus->isoms[gen_idx]);
#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], nbr));
  assert(isometry_is_isom(hash_isom, nbr, genus->genus_reps->keys[i]));
#endif // DEBUG_LEVEL_FULL
  isometry_muleq_right(s_nbr, hash_isom);
#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], genus->genus_reps->keys[i]));
  assert(isometry_is_isom(genus->isoms[i], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[i]));
#endif // DEBUG_LEVEL_FULL
  isometry_inv(s_inv, genus->isoms[i]);
#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_inv, genus->genus_reps->keys[i],
			  genus->genus_reps->keys[0]));
#endif // DEBUG_LEVEL_FULL
  isometry_muleq_right(s_nbr, s_inv);
#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0],
			  genus->genus_reps->keys[0]));
#endif // DEBUG_LEVEL_FULL
  
  *spin_vals = ((i << (genus->spinor->num_primes)) | spinor_norm_isom(genus->spinor, s_nbr));

  isometry_clear(s_nbr);
  isometry_clear(s_inv);
  isometry_clear(hash_isom);
  
  fmpz_mat_clear(nbr_fmpz);
  fmpz_mat_clear(nbr_isom);

  square_matrix_clear(nbr);
  
  return 0;
}

int process_isotropic_vector_all_conductors(neighbor_manager_t nbr_man, W64* spin_vals,
					    const genus_t genus,
					    double* theta_time, double* isom_time,
					    double* total_time, int* num_isom, int gen_idx)
{
  int i, j;
  clock_t cputime;
  square_matrix_t nbr /* , aut_tr */;
  isometry_t s_nbr, s_inv;
  isometry_t hash_isom;
  bool non_singular;
  int vec_cmp;
  slong orb_size, stab_size;
  isometry_t* orbit_isom;
  vector_t* orbit;
  bool found;
  W64 spinorm;
  vector_t g_vec;

  cputime = clock();
  
  assert(genus->genus_reps->num_stored > gen_idx);

  orbit = (vector_t*)malloc(nbr_man->num_auts * sizeof(vector_t));
  orbit_isom = (isometry_t*)malloc(nbr_man->num_auts * sizeof(isometry_t));
  isometry_init(orbit_isom[0]);

  // right now we grow the orbit. Could also compute on the fly and in the end divide by the stabilizer's size
  // weird that it isn't already reduced?
  // vector_mod_p(nbr_man->iso_vec, nbr_man->p);
  // this might be unnecessary, but we want to be on the safe side for now
  // normalize_mod_p(nbr_man->iso_vec, nbr_man->p);
  vector_set(orbit[0], nbr_man->iso_vec);
  // Applying the automorphism group to the isotropic vector
  stab_size = 0;
  orb_size = 1;
  for (i = 0; i < nbr_man->num_auts; i++) {
    square_matrix_mul_vec_left(g_vec, nbr_man->iso_vec, nbr_man->auts[i]);
    // vector_mod_p(g_vec, nbr_man->p);
    // normalize_mod_p(g_vec, nbr_man->p);
    // vec_cmp = vector_cmp(g_vec, nbr_man->iso_vec);
    vec_cmp = vector_cmp_unred(g_vec, nbr_man->iso_vec, nbr_man->p);
    if (vec_cmp < 0) {
      free(orbit);
      free(orbit_isom);
      return 0;
    }
    if (vec_cmp == 0)
      stab_size++;
    if (vec_cmp > 0) {
      found = false;
      for (j = 0; j < orb_size; j++)
	if (vector_cmp_unred(g_vec, orbit[j], nbr_man->p) == 0) {
	  found = true;
	  break;
	}
      if (!found) {
	vector_set(orbit[orb_size], g_vec);
	isometry_init_set_square_matrix(orbit_isom[orb_size], nbr_man->conj_auts[i], 1);
	orb_size++;
      }
    }
  }
  assert(orb_size == nbr_man->num_auts / stab_size);
  
  non_singular = nbr_process_build_nb_and_isom(nbr, s_nbr, nbr_man);

  if (!non_singular) {
    free(orbit);
    free(orbit_isom);
    return 0;
  }

#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
#endif // DEBUG_LEVEL_FULL
  

  i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, theta_time, isom_time, num_isom);

  assert( (i >= 0) && (i <= genus->dims[0]) );

  // TODO - compute only the upper half, and complete using hermitian property
  // Determine which subspaces this representative contributes.

  // TODO - reduce the number of actual computations required here

#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
  for (j = 0; j < orb_size; j++)
    assert(isometry_is_isom(orbit_isom[j], genus->genus_reps->keys[gen_idx], genus->genus_reps->keys[gen_idx]));
  assert(isometry_is_isom(genus->isoms[gen_idx], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[gen_idx]));
#endif // DEBUG_LEVEL_FULL
  
  isometry_muleq_left(s_nbr, genus->isoms[gen_idx]);

#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], nbr));
  for (j = 0; j < orb_size; j++)
     assert(isometry_is_isom(orbit_isom[j], genus->genus_reps->keys[0], genus->genus_reps->keys[gen_idx]));
#endif // DEBUG_LEVEL_FULL
  isometry_inv(s_inv, genus->isoms[gen_idx]);

#ifdef DEBUG_LEVEL_FULL
  for (j = 0; j < orb_size; j++)
    assert(isometry_is_isom(orbit_isom[j], genus->genus_reps->keys[0],
			    genus->genus_reps->keys[0]));
  
  assert(isometry_is_isom(hash_isom, nbr, genus->genus_reps->keys[i]));
#endif // DEBUG_LEVEL_FULL
  
  isometry_muleq_right(s_nbr, hash_isom);

#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], genus->genus_reps->keys[i]));
  assert(isometry_is_isom(genus->isoms[i], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[i]));
#endif // DEBUG_LEVEL_FULL
  isometry_inv(s_inv, genus->isoms[i]);

#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_inv, genus->genus_reps->keys[i],
			  genus->genus_reps->keys[0]));
#endif // DEBUG_LEVEL_FULL
  
  isometry_muleq_right(s_nbr, s_inv);

#ifdef DEBUG_LEVEL_FULL
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0],
			  genus->genus_reps->keys[0]));
#endif // DEBUG_LEVEL_FULL
  
  spinorm = spinor_norm_isom(genus->spinor, s_nbr);
  for (j = 0; j < orb_size; j++) {
    *(spin_vals+j) = ((i << (genus->spinor->num_primes)) | (spinorm ^ spinor_norm_isom(genus->spinor, orbit_isom[j])) );
  }

  isometry_clear(s_nbr);
  isometry_clear(s_inv);
  isometry_clear(hash_isom);

  square_matrix_clear(nbr);

  free(orbit);
  free(orbit_isom);

  (*total_time) += clock() - cputime;
  
  return orb_size;
}

int process_isotropic_vector(neighbor_manager_t nbr_man, int* T, const genus_t genus,
			     double* theta_time, double* isom_time, double* total_time, int* num_isom)

{

  int i;
  clock_t cputime;
  square_matrix_t nbr;
  bool non_singular;
  vector_t g_vec;
  int vec_cmp;
  slong orbit_size, stab_size;

  // weird that it isn't already reduced?
  //  vector_mod_p(nbr_man->iso_vec, nbr_man->p);
  // this might be unnecessary, but we want to be on the safe side for now
  // normalize_mod_p(nbr_man->iso_vec, nbr_man->p);
  // Applying the automorphism group to the isotropic vector
  stab_size = 0;
  for (i = 0; i < nbr_man->num_auts; i++) {
    square_matrix_mul_vec_left(g_vec, nbr_man->iso_vec, nbr_man->auts[i]);
    // vector_mod_p(g_vec, nbr_man->p);
    // normalize_mod_p(g_vec, nbr_man->p);
    vec_cmp = vector_cmp_unred(g_vec, nbr_man->iso_vec, nbr_man->p);
    if (vec_cmp < 0)
      return 0;
    if (vec_cmp == 0)
      stab_size++;
  }
  orbit_size = nbr_man->num_auts / stab_size;
  
  non_singular = nbr_process_build_nb(nbr, nbr_man);
  if (!non_singular)
    return 0;

  cputime = clock();
  i = hash_table_indexof(genus->genus_reps, nbr, 0, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

  assert((i >= 0) && (i <= genus->dims[0]));

  T[i] += orbit_size;
  
  return 1;
}



// !! TODO - use i to cut the parameters to chunks, need to convert from a number to a vector in the parameter space
int process_neighbour_chunk(int* T, int p, int i, int gen_idx, const genus_t genus,
			    double* theta_time, double* isom_time, double* total_time, int* num_isom)
{
  square_matrix_t Q;
  neighbor_manager_t nbr_man;
  int lc;
  
  lc = 0;

  assert(gen_idx < genus->genus_reps->num_stored);
  square_matrix_set(Q,genus->genus_reps->keys[gen_idx]);

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL

  nbr_process_init(nbr_man, Q, p, i, genus->isoms[gen_idx]);

  while (!(nbr_process_has_ended(nbr_man))) {
    lc += process_isotropic_vector(nbr_man, T, genus, theta_time, isom_time, total_time, num_isom);
    nbr_process_advance(nbr_man);
  }

  nbr_process_clear(nbr_man);
 
  return lc;
}

slong hecke_col_nbr_data_all_conductors(W64* spin_vals, int p, int k, int gen_idx, const genus_t genus)
{
  square_matrix_t Q;
  nbr_data_t nbr_man;
  int lc;
  int num_isom;
  double theta_time, isom_time, total_time;

  num_isom = lc = 0;
  theta_time = isom_time = total_time = 0;

  lc = 0;
  
  assert(gen_idx < genus->genus_reps->num_stored);
  square_matrix_set(Q,genus->genus_reps->keys[gen_idx]);

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL

  assert(k > 0);
  
  nbr_data_init(nbr_man, Q, p, k);

  while (!(nbr_data_has_ended(nbr_man))) {
    process_isotropic_vector_nbr_data_all_conductors(nbr_man, &(spin_vals[lc]), genus,
						     &theta_time, &isom_time, &total_time,
						     &num_isom, gen_idx);
    nbr_data_get_next_neighbor(nbr_man);
    lc++;
  }

#ifdef DEBUG
  printf("theta_time = %f, isom_time = %f, total_time = %f, num_isom = %d / %d \n",
	 theta_time/lc, isom_time/lc, total_time, num_isom, lc);
#endif // DEBUG
  
  nbr_data_clear(nbr_man);

  return lc;
}

slong hecke_col_all_conductors(W64* spin_vals, int p, int gen_idx, const genus_t genus)
{
  square_matrix_t Q;
  neighbor_manager_t nbr_man;
  int lc, i;
  int num_isom;
  double theta_time, isom_time, total_time;

  num_isom = lc = 0;
  theta_time = isom_time = total_time = 0;

  lc = 0;

  assert(gen_idx < genus->genus_reps->num_stored);
  square_matrix_set(Q,genus->genus_reps->keys[gen_idx]);

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL

  for (i = 0; i < p; i++) {
    nbr_process_init(nbr_man, Q, p, i, genus->isoms[gen_idx]);

    while (!(nbr_process_has_ended(nbr_man))) {
      lc += process_isotropic_vector_all_conductors(nbr_man, &(spin_vals[lc]), genus,
						    &theta_time, &isom_time, &total_time,
						    &num_isom, gen_idx);
      nbr_process_advance(nbr_man);
    }
    nbr_process_clear(nbr_man);
  }

#ifdef DEBUG
  printf("theta_time = %f, isom_time = %f, total_time = %f, num_isom = %d / %d \n",
	 theta_time/lc, isom_time/lc, total_time, num_isom, lc);
#endif // DEBUG

  return lc;
}

int process_neighbour_chunk_nbr_data(int* T, int p, int k, int gen_idx, const genus_t genus,
				     double* theta_time, double* isom_time, double* total_time, int* num_isom)
{
  square_matrix_t Q;
  nbr_data_t nbr_man;
  int lc;

  lc = 0;

  assert(gen_idx < genus->genus_reps->num_stored);
  square_matrix_set(Q,genus->genus_reps->keys[gen_idx]);

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  square_matrix_print(Q);
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
  slong i;
  int num;
  int num_isom, lc;
  double theta_time, isom_time, total_time;
  num_isom = lc = 0;
  theta_time = isom_time = total_time = 0;

  if (genus->dims[0] == 0)
    return;
  
  for (i = 0; i < genus->dims[0]; i++) {
    T[i] = 0;
  }
  
  for (num = 0; num < p; num++) {
    lc += process_neighbour_chunk(T, p, num, gen_idx, genus, &theta_time, &isom_time, &total_time, &num_isom);
  }

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

int** hecke_col_all_conds_sparse(int p, int col_idx, const genus_t genus)
{
  slong c, i;

  nbr_data_t nbr_man;
  
  W64* spin_vals;
  int** hecke;
  slong npos, rpos, num_nbrs;
  slong* lut;
  int* row;
  slong spin_idx;
  W64 r;
  bool only_trivial;

  hecke = (int**)malloc(genus->num_conductors * sizeof(int*));

  only_trivial = genus_is_trivial_cond(genus);

  for (c = 0; c < genus->num_conductors; c++) {
    hecke[c] = (int*)malloc(genus->dims[c] * sizeof(int));
    for (i = 0; i < genus->dims[c]; i++)
      hecke[c][i] = 0;
  }

  if (only_trivial) {
    hecke_col(hecke[0], p, col_idx, genus);
    return hecke;
  }

  assert(genus->genus_reps->num_stored > 0);
  nbr_data_init(nbr_man, genus->genus_reps->keys[0], p, 1);
  // in neighbor_manager we also go over the radical
  // Maybe fixed that
  num_nbrs = number_of_neighbors(nbr_man, false);
  nbr_data_clear(nbr_man);
  
  spin_vals = (W64*)malloc(num_nbrs*sizeof(W64));
  hecke_col_all_conductors(spin_vals, p, col_idx, genus);
  
  for (c = 0; c < genus->num_conductors; c++) {
    lut = genus->lut_positions[c];
    npos = lut[col_idx];
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
 
  return hecke;
}

matrix_TYP** hecke_matrices_nbr_data_all_conductors(const genus_t genus, int p, int k)
{
  matrix_TYP** hecke;
  nbr_data_t nbr_man;
  int gen_idx, rpos;
  slong c, npos, num_nbrs, spin_idx;
#ifdef DEBUG
  slong nnbrs;
#endif // DEBUG
  int* row;
  slong* lut;
  W64* spin_vals;
  W64 r;
  int* hecke_ptr;

  hecke = (matrix_TYP**)malloc(genus->num_conductors * sizeof(matrix_TYP*));
  hecke_ptr = (int*)malloc(genus->num_conductors * sizeof(int));

  for (c = 0; c < genus -> num_conductors; c++) {
    // hecke_ptr[c] = (slong*)malloc(genus->dims[c] * genus->dims[c] * sizeof(slong));
    hecke[c] = init_mat(genus->dims[c], genus->dims[c], "");
    hecke_ptr[c] = 0;
  }

  // just computing the number of neighbors, to collect all the spin values at once
  assert(genus->genus_reps->num_stored > 0);
  nbr_data_init(nbr_man, genus->genus_reps->keys[0], p, k);
  // nbr_data only goes over the isotropic vectors that are not in the radical
  num_nbrs = number_of_neighbors(nbr_man, false);
  spin_vals = (W64*)malloc(num_nbrs*sizeof(W64));
  nbr_data_clear(nbr_man);
  
  for (gen_idx = 0; gen_idx < genus->genus_reps->num_stored; gen_idx++) {
#ifdef DEBUG
    nnbrs = hecke_col_nbr_data_all_conductors(spin_vals, p, k, gen_idx, genus);
    assert(nnbrs == num_nbrs);
#else
    hecke_col_nbr_data_all_conductors(spin_vals, p, k, gen_idx, genus);
#endif // DEBUG
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

matrix_TYP** hecke_matrices_all_conductors(const genus_t genus, int p)
{
  matrix_TYP** hecke;
  nbr_data_t nbr_man;
  int gen_idx, rpos;
  slong c, npos, num_nbrs, spin_idx;
#ifdef DEBUG
  slong nnbrs;
#endif // DEBUG
  int* row;
  slong* lut;
  W64* spin_vals;
  W64 r;
  int* hecke_ptr;
  bool only_trivial;

  hecke = (matrix_TYP**)malloc(genus->num_conductors * sizeof(matrix_TYP*));

  only_trivial = genus_is_trivial_cond(genus);

  if (only_trivial)
    hecke[0] = hecke_matrix(genus, p);
  else
    hecke[0] = init_mat(genus->dims[0], genus->dims[0], "");

  for (c = 1; c < genus -> num_conductors; c++)
    hecke[c] = init_mat(genus->dims[c], genus->dims[c], "");

  if (only_trivial)
    return hecke;

  hecke_ptr = (int*)malloc(genus->num_conductors * sizeof(int));
  
  for (c = 0; c < genus -> num_conductors; c++)
    hecke_ptr[c] = 0;

  // just computing the number of neighbors, to collect all the spin values at once
  if (genus->dims[0] == 0)
    num_nbrs = 0;
  else {
    assert(genus->genus_reps->num_stored > 0);
    nbr_data_init(nbr_man, genus->genus_reps->keys[0], p, 1);
    num_nbrs = number_of_neighbors(nbr_man, false);
    nbr_data_clear(nbr_man);
  }
  
  spin_vals = (W64*)malloc(num_nbrs*sizeof(W64));
  
  for (gen_idx = 0; gen_idx < genus->genus_reps->num_stored; gen_idx++) {
#ifdef DEBUG
    nnbrs = hecke_col_all_conductors(spin_vals, p, gen_idx, genus);
    assert(nnbrs == num_nbrs);
#else
    hecke_col_all_conductors(spin_vals, p, gen_idx, genus);
#endif // DEBUG
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
  nf_elem_t e_new;
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

  pivot = ev_get_pivot(evs->eigenvecs[ev_idx], evs->nfs[ev_idx], genus, 0);
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
  // this line is leaky, so we do it cautiously
  nf_elem_init(e_new, evs->nfs[ev_idx]);
  nf_elem_div(e_new, e, evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]);
  nf_elem_set(e, e_new, evs->nfs[ev_idx]);

  nf_elem_clear(e_new, evs->nfs[ev_idx]);

  // handling the case p divides the discriminant
  if (genus->genus_reps->num_stored > 0) {
    fmpz_mat_init_set_square_matrix(q, genus->genus_reps->keys[0]);
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
  nf_elem_t e_new;
  nbr_data_t nbr_man;
  int** hecke;
  int i, pivot;
  slong num_nbrs, spin_idx;
#ifdef DEBUG
  slong nnbrs;
#endif // DEBUG
  slong c, npos, rpos, gen_idx;
  nf_elem_t prod;
  fmpz_mat_t q;
  fmpz_t disc;
  slong* lut;
  int* row;
  W64* spin_vals;
  W64 r;
#ifdef DEBUG
  int num;
  clock_t cpuclock;
  double cputime;
  
  cpuclock = clock();
#endif // DEBUG
  
  nf_elem_init(e, evs->nfs[ev_idx]);
  nf_elem_init(prod, evs->nfs[ev_idx]);

  pivot = ev_get_pivot(evs->eigenvecs[ev_idx], evs->nfs[ev_idx], genus, ev_cond);

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

  assert(genus->genus_reps->num_stored > 0);
  nbr_data_init(nbr_man, genus->genus_reps->keys[0], p, k);
  num_nbrs = number_of_neighbors(nbr_man, false);
  nbr_data_clear(nbr_man);

  spin_vals = (W64*)malloc(num_nbrs*sizeof(W64));
#ifdef DEBUG
  nnbrs = hecke_col_nbr_data_all_conductors(spin_vals, p, k, gen_idx, genus);
  assert(num_nbrs == nnbrs);
#else
  hecke_col_nbr_data_all_conductors(spin_vals, p, k, gen_idx, genus);
#endif // DEBUG
  
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
  // this line is leaky, so we do it cautiously
  nf_elem_init(e_new, evs->nfs[ev_idx]);
  nf_elem_div(e_new, e, evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]);
  nf_elem_set(e, e_new, evs->nfs[ev_idx]);

  nf_elem_clear(e_new, evs->nfs[ev_idx]);

  // handling the case p divides the discriminant
  if (genus->genus_reps->num_stored > 0) {
    fmpz_mat_init_set_square_matrix(q, genus->genus_reps->keys[0]);
    fmpz_init(disc);
    fmpz_mat_det(disc, q);
    fmpz_divexact_si(disc,disc,2);
    if (fmpz_get_si(disc) % p == 0) {
      nf_elem_add_si(e, e, (genus->conductors[ev_cond] % p == 0) ? -1 : 1, evs->nfs[ev_idx]);
    }
    fmpz_mat_clear(q);
    fmpz_clear(disc);
  }

#ifdef DEBUG
  printf("%4d ", p);
  nf_elem_print_pretty(e, evs->nfs[ev_idx], "a");
  printf(" - ");
  for (num = 0; num < genus->dims[ev_cond]; num++)
    printf("%10d ", hecke[ev_cond][num]);
  
  printf("- %10f\n", cputime);
#endif // DEBUG

  nf_elem_clear(prod, evs->nfs[ev_idx]);

  for (c = 0; c < genus->num_conductors; c++)
    free(hecke[c]);

  free(hecke);
  
  return;
}


void get_hecke_ev_all_conductors(nf_elem_t e, const genus_t genus,
				 const eigenvalues_t evs,
				 int p, int ev_idx, slong ev_cond)
{
  nf_elem_t e_new;
  nbr_data_t nbr_man;
  int** hecke;
  int i, pivot;
  slong num_nbrs, spin_idx;
#ifdef DEBUG
  slong nnbrs;
#endif // DEBUG
  slong c, npos, rpos, gen_idx;
  nf_elem_t prod;
  fmpz_mat_t q;
  fmpz_t disc;
  slong* lut;
  int* row;
  W64* spin_vals;
  W64 r;
#ifdef DEBUG
  int num;
  clock_t cpuclock;
  double cputime;
  
  cpuclock = clock();
#endif // DEBUG
  
  nf_elem_init(e, evs->nfs[ev_idx]);
  nf_elem_init(prod, evs->nfs[ev_idx]);

  pivot = ev_get_pivot(evs->eigenvecs[ev_idx], evs->nfs[ev_idx], genus, ev_cond);

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

  assert(genus->genus_reps->num_stored > 0);
  nbr_data_init(nbr_man, genus->genus_reps->keys[0], p, 1);
  num_nbrs = number_of_neighbors(nbr_man, false);
  nbr_data_clear(nbr_man);

  spin_vals = (W64*)malloc(num_nbrs*sizeof(W64));
#ifdef DEBUG
  nnbrs = hecke_col_all_conductors(spin_vals, p, gen_idx, genus);
  assert(num_nbrs == nnbrs);
#else
  hecke_col_all_conductors(spin_vals, p, gen_idx, genus);
#endif // DEBUG
  
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
  // this line is leaky, so we do it cautiously
  nf_elem_init(e_new, evs->nfs[ev_idx]);
  nf_elem_div(e_new, e, evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]);
  nf_elem_set(e, e_new, evs->nfs[ev_idx]);

  nf_elem_clear(e_new, evs->nfs[ev_idx]);

  // handling the case p divides the discriminant
  if (genus->genus_reps->num_stored > 0) {
    fmpz_mat_init_set_square_matrix(q, genus->genus_reps->keys[0]);
    fmpz_init(disc);
    fmpz_mat_det(disc, q);
    fmpz_divexact_si(disc,disc,2);
    if (fmpz_get_si(disc) % p == 0) {
      nf_elem_add_si(e, e, 1, evs->nfs[ev_idx]);
      nf_elem_scalar_mul_si(e, e, (genus->conductors[ev_cond] % p == 0) ? -1 : 1, evs->nfs[ev_idx]); 
    }
    fmpz_mat_clear(q);
    fmpz_clear(disc);
  }

#ifdef DEBUG
  printf("%4d ", p);
  nf_elem_print_pretty(e, evs->nfs[ev_idx], "a");
  printf(" - ");
  for (num = 0; num < genus->dims[ev_cond]; num++)
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
  nf_elem_t e_new;
  int* a;
  int num, i, pivot;
  nf_elem_t prod;
  fmpz_mat_t q;
  fmpz_t disc;

  nf_elem_init(e, evs->nfs[ev_idx]);
  nf_elem_init(prod, evs->nfs[ev_idx]);
  a = (int*)malloc(genus->genus_reps->num_stored * sizeof(int));
  for (num = 0; num < genus->genus_reps->num_stored; num++)
    a[num] = 0;

  pivot = ev_get_pivot(evs->eigenvecs[ev_idx], evs->nfs[ev_idx], genus, 0);
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
  // this line is leaky, so we do it cautiously
  nf_elem_init(e_new, evs->nfs[ev_idx]);
  nf_elem_div(e_new, e, evs->eigenvecs[ev_idx][pivot], evs->nfs[ev_idx]);
  nf_elem_set(e, e_new, evs->nfs[ev_idx]);

  nf_elem_clear(e_new, evs->nfs[ev_idx]);

#ifdef DEBUG
  printf("%4d ", p);
  nf_elem_print_pretty(e, evs->nfs[ev_idx], "a");
  printf(" - ");
  for (num = 0; num < genus->genus_reps->num_stored; num++)
    printf("%10d ", a[num]);
  
  printf("- %10f\n", cputime);
#endif // DEBUG

  if (genus->genus_reps->num_stored > 0) {
    fmpz_mat_init_set_square_matrix(q, genus->genus_reps->keys[0]);
    fmpz_init(disc);
    fmpz_mat_det(disc, q);
    fmpz_divexact_si(disc,disc,2);
    if (fmpz_get_si(disc) % p == 0) {
      nf_elem_add_si(e, e, 1, evs->nfs[ev_idx]);
    }
    fmpz_mat_clear(q);
    fmpz_clear(disc);
  }

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
  hecke_matrices = hecke_matrices_nbr_data_all_conductors(genus, p, k);
  for (c = 0; c < genus->num_conductors; c++) {
    fmpq_mat_init_set_matrix_TYP(hecke_fmpq_mat[c], hecke_matrices[c]);
    // we transpose to match with our conventions
    fmpq_mat_transpose(hecke_fmpq_mat[c], hecke_fmpq_mat[c]);
    free_mat(hecke_matrices[c]);
  }
  free(hecke_matrices);
  
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
    } while ((!ev_init) && (p_idx < D->num_primes));
    if (ev_init) {
      nf_elem_init(evs->eigenvals[i], evs->nfs[i]);
      nf_elem_gen(evs->eigenvals[i], evs->nfs[i]);
      evs->lift_type[i] = G; // default
    }
    else {
      nf_elem_init(evs->eigenvals[i], evs->nfs[i]);
      nf_elem_gen(evs->eigenvals[i], evs->nfs[i]);
      evs->lift_type[i] = O; // it's an old form
    }
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
    fprintf(stderr, "decomposing conductor %ld\n", genus->conductors[c]);
#endif // DEBUG
    decompose(D, genus, c);
#ifdef DEBUG
    fprintf(stderr, "computing eigenvectors\n");
#endif // DEBUG
    hecke_eigenforms(all_evs[c], D, genus, c);
  }

  decomposition_clear(D);
  return all_evs;
}
