/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : genus.c
 *
 * Description: Genus of a lattice.
 *              Used also to record relevant data about
 *              the space of orthogonal modular forms
 *              supported on this genus with different
 *              conductors.
 *
 **********************************************************
 */

// System dependencies

#include <assert.h>

// Required packages dependencies

#include <carat/matrix.h>

#include <flint/fmpq.h>
#include <flint/fmpz.h>

// Self dependencies

#include "aut_grp.h"
#include "genus.h"
#include "hash.h"
#include "io.h"
#include "isometry.h"
#include "mass.h"
#include "matrix_tools.h"
#include "nbr_data.h"
#include "neighbor.h"
#include "square_matrix.h"
#include "typedefs.h"

// a look up table of bit counts for numbers up to 8 bits, for quick popcount

int bitcounts[256] = {
  0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2,
  3, 3, 4, 3, 4, 4, 5, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3,
  3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3,
  4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4,
  3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5,
  6, 6, 7, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4,
  4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5,
  6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 2, 3, 3, 4, 3, 4, 4, 5,
  3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 3,
  4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6,
  6, 7, 6, 7, 7, 8
};

/* computes the popcount of a 64-bit word (number of lit bits) */
int popcnt(W64 x)
{
  if (x <= 0xff) return bitcounts[x];
  
  int count = 0;
  while (x)
    {
      count += bitcounts[x&0xff];
      x >>= 8;
    }
  return count;
}

/* compute the discriminant of a lattice from its gram matrix  */
void disc_init_set(fmpz_t disc, const square_matrix_t q)
{
  fmpz_mat_t q_fmpz;
  
  fmpz_init(disc);
  fmpz_mat_init_set_square_matrix(q_fmpz, q);
  fmpz_mat_det(disc, q_fmpz);
  fmpz_divexact_si(disc, disc, 2);
  
  fmpz_mat_clear(q_fmpz);
  
  return;
}

/* Estimate the number of genus representatives based on its mass.
   This is only a rough estimate, in order to allocate memory initially */
size_t genus_size_estimate(const fmpq_t mass)
{
  fmpz_t genus_size_fmpz;
  size_t genus_size;
  
  fmpz_init(genus_size_fmpz);
  fmpz_set(genus_size_fmpz, fmpq_denref(mass));
  genus_size = fmpz_get_si(genus_size_fmpz);
  fmpz_clear(genus_size_fmpz);
  genus_size = (4 * genus_size) / 3; // load factor

  return genus_size;
}

/* Initialize the possible conductors for the spinor norm (all square-free divisors of the discriminant).
   Allocates memory and initializes values for the relevant parts of the genus structure */
void conductors_init(genus_t genus)
{
  slong c, value;
  size_t bits, mask;
  
  genus->num_conductors = 1LL << genus->spinor->num_primes;
  
  genus->conductors = (slong*)malloc((genus->num_conductors) * sizeof(slong));
  genus->conductors[0] = 1;
  
  bits = 0;
  mask = 1;
  for (c = 1; c < genus->num_conductors; c++) {
    if (c == 2*mask) {
      ++bits;
      mask = 1LL << bits;
    }
    value = genus->spinor->primes[bits].n;
    if (value == 4)
      value = 2;
    value *= genus->conductors[c ^ mask];
    genus->conductors[c] = value;
  }  

  return;
}

// Determine which subspaces this representative contributes.
// sets the ignore flags and returns the order of the automorphism group
slong ignore_set(bool* ignore, const genus_t genus, slong genus_idx)
{
  slong c, gen_idx, order;
  isometry_t s;
  aut_grp_t aut_grp;
  W64 vals;

  aut_grp_init_square_matrix(aut_grp, genus->genus_reps->keys[genus_idx]);
  
  for (c = 0; c < genus->num_conductors; c++)
    ignore[c] = false;

  for (gen_idx = 0; gen_idx < aut_grp->num_gens; gen_idx++) {
#ifdef DEBUG_LEVEL_FULL
    isometry_init_set_square_matrix(s, aut_grp->gens[gen_idx], 1);
    assert(isometry_is_isom(s, genus->genus_reps->keys[genus_idx], genus->genus_reps->keys[genus_idx]));
    isometry_clear(s);
#endif // DEBUG_LEVEL_FULL

    isometry_init(s);
    isometry_inv(s, genus->isoms[genus_idx]);
    isometry_mul_mateq_left(s, aut_grp->gens[gen_idx], 1);
    isometry_muleq_left(s, genus->isoms[genus_idx]);

#ifdef DEBUG_LEVEL_FULL
    assert(isometry_is_isom(s, genus->genus_reps->keys[0], genus->genus_reps->keys[0]));
#endif // DEBUG_LEVEL_FULL
	
    vals = spinor_norm_isom(genus->spinor, s);
    // !! TODO - we can break the loop after we find one, right?
    for (c = 0; c < genus->num_conductors; c++)
      if (!ignore[c] && (popcnt(vals & c) & 1))
	ignore[c] = true;
    isometry_clear(s);	
  }

  order = aut_grp->order;

  aut_grp_clear(aut_grp);
  
  return order;
}

/* Compute the dimensions of each space of orthogonal modular forms based on this genus.
   Also initializes the look up tables lut_positions and the number of automorphisms num_auts accordingly. */
void dimensions_init(genus_t genus)
{
  slong c, genus_idx, order;
  bool* ignore;
  
  // Initialize the dimensions to zero, we will compute these values below.
  genus->dims = (slong*)malloc((genus->num_conductors)*sizeof(slong));
  for (c = 0; c < genus->num_conductors; c++)
    genus->dims[c] = 0;

  // Create the lookup table values for each genus rep at each conductor.
  genus->lut_positions = (slong**)malloc((genus->num_conductors)*sizeof(slong*));
  for (c = 0; c < genus->num_conductors; c++) {
    genus->lut_positions[c] = (slong*)malloc((genus->genus_reps->num_stored) * sizeof(slong));
    for (genus_idx = 0; genus_idx < genus->genus_reps->num_stored; genus_idx++)
      genus->lut_positions[c][genus_idx] = -1;
  }

  // allocate the ignore flags and the sizes of automorphism groups
  ignore = (bool*)malloc((genus->num_conductors)*sizeof(bool));
  genus->num_auts = (slong**)malloc((genus->num_conductors)*sizeof(slong*));
  for (c = 0; c < genus->num_conductors; c++)
    genus->num_auts[c] = (slong*)malloc((genus->genus_reps->num_stored)*sizeof(slong));

  for (genus_idx = 0; genus_idx < genus->genus_reps->num_stored; genus_idx++) {
    order = ignore_set(ignore, genus, genus_idx);

    for (c = 0; c < genus->num_conductors; c++) {
      if (!ignore[c]) {
	genus->lut_positions[c][genus_idx] = genus->dims[c];
	genus->num_auts[c][genus->dims[c]] = order;
      }
      genus->dims[c] += (ignore[c] ? 0 : 1);
    }

  }

  free(ignore);
  
  return;
}

/* compute the genus of a quadratic form */
void genus_init_square_matrix(genus_t genus, const square_matrix_t q, int h)
{
  aut_grp_t aut_grp;
  square_matrix_t nbr, genus_rep;
  isometry_t s, isom;
  fmpq_t mass, acc_mass, mass_form;
  fmpz_t prime;
  int p, current, key_num;
  size_t genus_size;
  hash_table_t slow_genus;
  bool found, is_isom, genus_full;
  
#ifndef NBR_DATA
  neighbor_manager_t nbr_man;
  int i;
#else
  nbr_data_t nbr_man;
  fmpz_mat_t nbr_isom, nbr_fmpz;
#endif // NBR_DATA

#ifdef NBR_DATA
  fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
  fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
#endif // NBR_DATA
  
#ifdef NBR_DATA
  fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
  fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
#endif // NBR_DATA

  disc_init_set(genus->disc, q);

  fmpq_init(mass);
  get_mass(mass, q);
  
#ifdef DEBUG_LEVEL_FULL
  printf("mass = ");
  fmpq_print(mass);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  spinor_init_square_matrix(genus->spinor, q);

  genus_size = genus_size_estimate(mass);
  
  hash_table_init(slow_genus, genus_size);
  hash_table_add(slow_genus, q);

  genus->isoms = (isometry_t*)malloc((slow_genus->capacity) * sizeof(isometry_t));

  isometry_init(genus->isoms[0]);
  
  // initializing the conductors
  
  conductors_init(genus);

  aut_grp_init_square_matrix(aut_grp, q);
  
  fmpq_init(mass_form);
  fmpq_init(acc_mass);
  fmpq_set_si(acc_mass, 1, aut_grp->order);

  aut_grp_clear(aut_grp);

#ifdef DEBUG_LEVEL_FULL
  printf("acc_mass = ");
  fmpq_print(acc_mass);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  fmpz_init(prime);
  fmpz_set_ui(prime, 1);

  if (h == -1)
    genus_full = !fmpq_cmp(acc_mass, mass);
  else
    genus_full = (h == slow_genus->num_stored);
  
  while (!genus_full) {
    current = 0;
#ifdef DEBUG
    if (h == -1)
      assert(fmpq_cmp(acc_mass, mass) <= 0);
#endif // DEBUG
    /* !! TODO - we don't really need to restrict to good primes here, */
    /* but let's check these first */
    do {
      fmpz_nextprime(prime, prime, true);
      p = fmpz_get_ui(prime);
    }
    while (square_matrix_is_bad_prime(q,p));

    while ((current < slow_genus->num_stored) && (!genus_full)){
#ifdef DEBUG_LEVEL_FULL
      printf("current = %d\n", current);
#endif // DEBUG_LEVEL_FULL

#ifndef NBR_DATA
      i = 0;
      
      /* Right now all the isotropic vectors are paritioned to p */
      /* sets, by index as Gonzalo did */
      
      while ((i < p) && (!genus_full)) {
	nbr_process_init(nbr_man, slow_genus->keys[current], p, i, genus->isoms[current]);
	while ((!(nbr_process_has_ended(nbr_man))) && (!genus_full)) {
	  nbr_process_build_nb_and_isom(nbr, s, nbr_man);
#else
	nbr_data_init(nbr_man, slow_genus->keys[current], p, 1);
	while ((!(nbr_data_has_ended(nbr_man))) && (!genus_full)) {
	  nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
	  square_matrix_set_fmpz_mat(nbr, nbr_fmpz);
#endif // NBR_DATA
	  
#ifdef DEBUG_LEVEL_FULL
	  printf("nbr = \n");
	  square_matrix_print(nbr);
#endif // DEBUG_LEVEL_FULL
	  
	  key_num = -1;

	  found = true;
	  is_isom = false;
      
	  while ((found) && (!is_isom)) {
#ifdef DEBUG_LEVEL_FULL
	    printf("checking if it is already in the genus...\n");
#endif // DEBUG_LEVEL_FULL
	    found = hash_table_get_key(genus_rep, slow_genus, nbr, &key_num);
	    if (found) {
#ifdef DEBUG_LEVEL_FULL
	      printf("Found candidate :\n");
	      square_matrix_print(genus_rep);
	      printf("Checking for isometry...\n");
#endif // DEBUG_LEVEL_FULL
	      is_isom = is_isometric(isom, genus_rep, nbr);
	    }
	  }
	  
	  if (!is_isom) {
#ifdef DEBUG_LEVEL_FULL
	    printf("no Isometry found, adding neighbor...\n");
#endif // DEBUG_LEVEL_FULL
#ifdef NBR_DATA
	    isometry_init_set_fmpz_mat(s, nbr_isom, p); 
#endif // NBR_DATA
#ifdef DEBUG_LEVEL_FULL
	    assert(isometry_is_isom(s, slow_genus->keys[current], nbr));
#endif // DEBUG_LEVEL_FULL
	    
	    greedy(nbr, s, QF_RANK);
#ifdef DEBUG_LEVEL_FULL
	    assert(isometry_is_isom(s, slow_genus->keys[current], nbr));
#endif // DEBUG_LEVEL_FULL

	    // rellocating to have enough room
	    if (slow_genus->num_stored == slow_genus->capacity)
	      genus->isoms = (isometry_t*)realloc(genus->isoms, (slow_genus->capacity << 1) * sizeof(isometry_t));
	    
	    hash_table_add(slow_genus, nbr);
	    
	    // The genus rep isometries were initialized only to contain the
	    // isometry between the parent and its child, we now want to update
	    // these isometries so that they are rational isometries between the
	    // "mother" quadratic form and the genus rep.
	    isometry_mul(genus->isoms[slow_genus->num_stored-1], genus->isoms[current], s);

#ifdef DEBUG_LEVEL_FULL
	    assert(isometry_is_isom(genus->isoms[slow_genus->num_stored-1], q, nbr));
#endif // DEBUG_LEVEL_FULL
	    
	    isometry_clear(s);
	    
	    aut_grp_init_square_matrix(aut_grp, nbr);
	    fmpq_set_si(mass_form, 1, aut_grp->order);
	    aut_grp_clear(aut_grp);
	    fmpq_add(acc_mass, acc_mass, mass_form);
	    if (h == -1)
	      genus_full = !fmpq_cmp(acc_mass, mass);
	    else
	      genus_full = (h == slow_genus->num_stored);
	    
#ifdef DEBUG_LEVEL_FULL
	    printf("acc_mass = ");
	    fmpq_print(acc_mass);
	    printf("\n");
#endif // DEBUG_LEVEL_FULL
	  }
	  
#ifndef NBR_DATA
	  nbr_process_advance(nbr_man);
#else
	  nbr_data_get_next_neighbor(nbr_man);
#endif // NBR_DATA
	}
	
#ifndef NBR_DATA
	i++;
	nbr_process_clear(nbr_man);
	}
#else
	nbr_data_clear(nbr_man);
#endif // NBR_DATA
	current++;
      }
      
    }

    hash_table_recalibrate(genus->genus_reps, slow_genus);
    
    hash_table_clear(slow_genus);
    
    assert(genus->genus_reps->num_stored > 0);
    
    dimensions_init(genus);
    
    fmpq_clear(mass);
    fmpq_clear(acc_mass);
    fmpq_clear(mass_form);
    fmpz_clear(prime);
    
#ifdef NBR_DATA
    fmpz_mat_clear(nbr_fmpz);
    fmpz_mat_clear(nbr_isom);
#endif // NBR_DATA
  
    return;
}

/* clear the memory allocated for the genus */
void genus_clear(genus_t genus)
{
  slong c;

  fmpz_clear(genus->disc);
  
  for (c = 0; c < genus->num_conductors; c++)
    free(genus->num_auts[c]);
  for (c = 0; c < genus->genus_reps->num_stored; c++)
    isometry_clear(genus->isoms[c]);
  
  if (genus->genus_reps->num_stored > 0)
    free(genus->isoms);
  
  for (c = 0; c < genus->num_conductors; c++)
    free(genus->lut_positions[c]);

  if (genus->num_conductors > 0) {
    free(genus->num_auts);
    free(genus->lut_positions);
    free(genus->dims);
    free(genus->conductors);
    spinor_clear(genus->spinor);
  }
  
  hash_table_clear(genus->genus_reps);
  return;
}

/* Check whether two lattices, given by their gram matrices, are isometric over Q.
   If so, construct an isometry between them. */
// !! TODO - right now we are doing it very naively, by constructing neighbors
// We should replace it by proper Gram-Schmidt
bool square_matrix_is_Q_isometric(isometry_t isom, const square_matrix_t A, const square_matrix_t B)
{
  isometry_t isom_A, isom_B, isom_B_inv;
  neighbor_manager_t nbr_man;
  square_matrix_t nbr;
  int i, p;
  fmpz_t prime;
  bool found = false;
  
  fmpz_init(prime);
  isometry_init(isom_A);
  isometry_init(isom_B);
  isometry_init(isom_B_inv);
  square_matrix_init(nbr);

  p = 2;
  i = 0;
  while (!found) {
    nbr_process_init(nbr_man, A, p, i, isom_A);
    while (!nbr_process_has_ended(nbr_man)) {
      nbr_process_build_nb_and_isom(nbr, isom_A, nbr_man);
      found = is_isometric(isom_B, B, nbr);
#ifdef DEBUG_LEVEL_FULL
      assert(isometry_is_isom(isom_A, A, nbr));
#endif // DEBUG_LEVEL_FULL
      if (found) {
	isometry_inv(isom_B_inv, isom_B);
#ifdef DEBUG_LEVEL_FULL
	assert(isometry_is_isom(isom_B, B, nbr));
	assert(isometry_is_isom(isom_B_inv, nbr, B));
#endif // DEBUG_LEVEL_FULL
	isometry_mul(isom, isom_A, isom_B_inv);
#ifdef DEBUG_LEVEL_FULL
	assert(isometry_is_isom(isom, A, B));
#endif // DEBUG_LEVEL_FULL
	break;
      }
      nbr_process_advance(nbr_man);
    }
    nbr_process_clear(nbr_man);
    i++;
    if (i == p) {
      fmpz_nextprime(prime, prime, true);
      p = fmpz_get_si(prime);
      i = 0;
    }
  }

  square_matrix_clear(nbr);
  isometry_clear(isom_A);
  isometry_clear(isom_B);
  isometry_clear(isom_B_inv);
  fmpz_clear(prime);
  
  return found;
}
 
/* Find isometries from the lattice whose gram matrix is A to each of the h lattices whose gram matrices are the
   elements of the list B. Return the isometries in isoms, for which memory is assumed to be allocated. */
// for now, this is faster (do one round and collect all)
// assumes isoms are initialized and that B and isoms are the length of the class number
void square_matrix_find_isometries(isometry_t* isoms, const square_matrix_t A, const square_matrix_t* B, int h)
{
  genus_t genus;
  isometry_t isom_B, isom_B_inv;
  slong i, rep_idx;
  double theta_time, isom_time;
  int num_isom = 0;

  theta_time = isom_time = 0;
  
  genus_init_square_matrix(genus, A, h);
  assert(genus->dims[0] == h);
  isometry_init(isom_B);

  for (i = 0; i < genus->dims[0]; i++) {
    rep_idx = hash_table_index_and_isom(genus->genus_reps, B[i], isom_B, &theta_time, &isom_time, &num_isom);
#ifdef DEBUG_LEVEL_FULL
    assert(isometry_is_isom(isom_B, B[i], genus->genus_reps->keys[rep_idx]));
    assert(isometry_is_isom(genus->isoms[rep_idx], A, genus->genus_reps->keys[rep_idx]));
#endif // DEBUG_LEVEL_FULL
    isometry_inv(isom_B_inv, isom_B);
    isometry_mul(isoms[i], genus->isoms[rep_idx], isom_B_inv);
#ifdef DEBUG_LEVEL_FULL
    assert(isometry_is_isom(isoms[i], A, B[i]));
#endif // DEBUG_LEVEL_FULL
  }
  
  isometry_clear(isom_B);
  genus_clear(genus);
  
  return;
}

/* Initialize a genus from a list of gram matrices for the genus representatives */
void genus_init_set_square_matrix_vec(genus_t genus, const square_matrix_t* reps, size_t h)
{
#ifdef DEBUG
  fmpq_t mass, acc_mass, mass_form;
#endif // DEBUG
#ifdef DEBUG_LEVEL_FULL
  aut_grp_t aut_grp;
#endif // DEBUG_LEVEL_FULL
  size_t genus_size;
  hash_table_t slow_genus;
  
  size_t i;

  assert(h > 0);
  disc_init_set(genus->disc, reps[0]);

  // Here we don't really need it, just here for verification
#ifdef DEBUG
  fmpq_init(mass);
  get_mass(mass, reps[0]);
#ifdef DEBUG_LEVEL_FULL
  printf("mass = ");
  fmpq_print(mass);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  fmpq_init(mass_form);
  fmpq_init(acc_mass);
  fmpq_zero(acc_mass);
#endif // DEBUG
    
  genus_size = (4 * h) / 3; // load factor

  hash_table_init(slow_genus, genus_size);

  spinor_init_square_matrix(genus->spinor, reps[0]);
  
  for (i = 0; i < h; i++) {
#ifdef DEBUG
    assert(fmpq_cmp(acc_mass, mass) <= 0);
#endif // DEBUG
    hash_table_add(slow_genus, reps[i]);
    
#ifdef DEBUG_LEVEL_FULL
    aut_grp_init_square_matrix(aut_grp, reps[i]);
    fmpq_set_si(mass_form, 1, aut_grp->order);
    aut_grp_clear(aut_grp);
    fmpq_add(acc_mass, acc_mass, mass_form);
    printf("acc_mass = ");
    fmpq_print(acc_mass);
    printf("\n");
#endif // DEBUG_LEVEL_FULL
  }

  hash_table_recalibrate(genus->genus_reps, slow_genus);

  hash_table_clear(slow_genus);

  // constructing isometries between the forms
  genus->isoms = (isometry_t*)malloc(h * sizeof(isometry_t));

  square_matrix_find_isometries(genus->isoms, reps[0], reps, h);
  
  // initializing the conductors
  
  conductors_init(genus);

  assert(genus->genus_reps->num_stored > 0);

  dimensions_init(genus);

#ifdef DEBUG  
  fmpq_clear(mass);
  fmpq_clear(acc_mass);
  fmpq_clear(mass_form);
#endif // DEBUG
      
  return;
}

/* Initialize a genus from a list of gram matrices for the genus representatives, and a list of specific isometries
   from the firt representative to all the others. */
void genus_init_set_square_matrix_vec_and_isoms(genus_t genus, const square_matrix_t* reps,
						const isometry_t* isoms, size_t h)
{
  size_t genus_size;
  hash_table_t slow_genus;
  
  size_t i;

  assert(h > 0);
  disc_init_set(genus->disc, reps[0]);

  genus_size = (4 * h) / 3; // load factor

  hash_table_init(slow_genus, genus_size);

  spinor_init_square_matrix(genus->spinor, reps[0]);
  
  for (i = 0; i < h; i++)
    hash_table_add(slow_genus, reps[i]);

  hash_table_recalibrate(genus->genus_reps, slow_genus);

  hash_table_clear(slow_genus);

  // constructing isometries between the forms
  genus->isoms = (isometry_t*)malloc(h * sizeof(isometry_t));
  for (i = 0; i < h; i++) {
    isometry_init_set(genus->isoms[i], isoms[i]);
#ifdef DEBUG_LEVEL_FULL
    assert(isometry_is_isom(genus->isoms[i], genus->genus_reps->keys[0], genus->genus_reps->keys[i]));
#endif // DEBUG_LEVEL_FULL
  }
  
  // initializing the conductors
  
  conductors_init(genus);

  assert(genus->genus_reps->num_stored > 0);

  dimensions_init(genus);
      
  return;
}

/* Initialize the genus from a file */
void genus_init_file(genus_t genus, const char* genus_fname, size_t disc, bool with_isom)
{
  square_matrix_t* genus_reps = NULL;
  isometry_t* isoms = NULL;
  size_t genus_size;

  if (with_isom)
    genus_size = read_genus_and_isom(&genus_reps, &isoms, genus_fname, disc);
  else
    genus_size = read_genus(&genus_reps, genus_fname, disc);
  
  if ((genus_size > 0) && (with_isom))
    genus_init_set_square_matrix_vec_and_isoms(genus, genus_reps, isoms, genus_size);
  else if (genus_size > 0)
    genus_init_set_square_matrix_vec(genus, genus_reps, genus_size);
  else
    genus_init_empty(genus, disc);
  
  free(genus_reps);
  free(isoms);
  return;
}

/* Initialize an empty genus, only from the discriminant */
void genus_init_empty(genus_t genus, size_t disc)
{
  fmpz_init_set_ui(genus->disc,disc);

  hash_table_init(genus->genus_reps, 0);

  spinor_init_fmpz(genus->spinor, genus->disc);

  conductors_init(genus);

  dimensions_init(genus);
  
  return;
}

/* check whether the genus is trivial */
bool genus_is_trivial_cond(const genus_t genus)
{
  bool only_trivial;
  slong c;
  
  only_trivial = true;
  for (c = 1; c < genus->num_conductors; c++) {
    if ((c > 0) && (genus->dims[c] > 0))
      only_trivial = false;
  }

  return only_trivial;
}
