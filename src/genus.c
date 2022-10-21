#include <assert.h>
#include <carat/matrix.h>

#include <flint/fmpq.h>
#include <flint/fmpz.h>

#include "genus.h"
#include "hash.h"
#include "mass.h"
#include "matrix_tools.h"
#include "nbr_data.h"
#include "neighbor.h"
#include "typedefs.h"

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

/* compute the genus of a quadratic form */
void genus_init(genus_t genus, matrix_TYP* Q)
{
  bravais_TYP *aut_grp;
  matrix_TYP *nbr, *isom, *genus_rep;
  isometry_t s;
  fmpq_t mass, acc_mass, mass_form;
  fmpz_t prime;
  int p, current, key_num;
  size_t genus_size, genus_idx, gen_idx;
  fmpz_t genus_size_fmpz;
  hash_table_t slow_genus;
  fmpz_mat_t q_fmpz;
  size_t c, mask, bits;
  slong value;
  int n;
  bool* ignore;
  W64 vals;
  
#ifndef NBR_DATA
  neighbor_manager nbr_man;
  int i;
#else
  nbr_data_t nbr_man;
  fmpz_mat_t nbr_isom, nbr_fmpz;
#endif // NBR_DATA

  n = Q->rows;
  assert(n == Q->cols);

#ifdef NBR_DATA
  fmpz_mat_init(nbr_isom, n, n);
  fmpz_mat_init(nbr_fmpz, n, n);
#endif // NBR_DATA
  
#ifdef NBR_DATA
  fmpz_mat_init(nbr_isom, n, n);
  fmpz_mat_init(nbr_fmpz, n, n);
#endif // NBR_DATA

  fmpz_init(genus->disc);

  fmpq_init(mass);
  get_mass(mass, Q);
  
#ifdef DEBUG
  printf("mass = ");
  fmpq_print(mass);
  printf("\n");
#endif // DEBUG

  /* this is ceiling */

  fmpz_init(genus_size_fmpz);
  // fmpz_cdiv_q(genus_size_fmpz, fmpq_numref(mass), fmpq_denref(mass));
  fmpz_set(genus_size_fmpz, fmpq_denref(mass));

  genus_size = fmpz_get_si(genus_size_fmpz);

  fmpz_clear(genus_size_fmpz);
    
  genus_size = (4 * genus_size) / 3; // load factor

  hash_table_init(slow_genus, genus_size);
  hash_table_add(slow_genus, Q);

  fmpz_mat_init_set_matrix_TYP(q_fmpz, Q);
  fmpz_mat_det(genus->disc, q_fmpz);
  fmpz_divexact_si(genus->disc, genus->disc, 2);
  spinor_init(genus->spinor, q_fmpz);
  fmpz_mat_clear(q_fmpz);

  genus->isoms = (isometry_t*)malloc((slow_genus->capacity) * sizeof(isometry_t));

  isometry_init(genus->isoms[0], init_mat(n, n, "1"), 1);
  
  // initializing the conductors
  
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
    value = fmpz_get_si(fq_nmod_ctx_prime(genus->spinor->fields[bits])) * genus->conductors[c ^ mask];
    genus->conductors[c] = value;
  }  
  
  aut_grp = automorphism_group(Q);

  fmpq_init(mass_form);
  fmpq_init(acc_mass);
  fmpq_set_si(acc_mass, 1, aut_grp->order);

#ifdef DEBUG_LEVEL_FULL
  printf("acc_mass = ");
  fmpq_print(acc_mass);
  printf("\n");
#endif // DEBUG_LEVEL_FULL
  
  fmpz_init(prime);
  fmpz_set_ui(prime, 1);
  
  while (fmpq_cmp(acc_mass, mass)) {
    current = 0;
#ifdef DEBUG
    if (fmpq_cmp(acc_mass, mass) > 0) {
      printf("Error! Accumulated too much mass!\n");
      return;
    }
#endif // DEBUG
    /* !! TODO - we don't really need to restrict to good primes here, */
    /* but let's check these first */
    do {
      fmpz_nextprime(prime, prime, true);
      p = fmpz_get_ui(prime);
/* #ifdef DEBUG */
/*       printf("p = %d, Q = ", p); */
/*       print_mat(Q); */
/* #endif //DEBUG */
    }
    while (!p_mat_det(Q, p));

    while ((current < slow_genus->num_stored) && fmpq_cmp(acc_mass, mass)){
#ifdef DEBUG_LEVEL_FULL
      printf("current = %d\n", current);
#endif // DEBUG_LEVEL_FULL

#ifndef NBR_DATA
      i = 0;
      
      /* Right now all the isotropic vectors are paritioned to p */
      /* sets, by index as Gonzalo did */
      
      while ((i < p) && fmpq_cmp(acc_mass, mass)) {
	/* printf("i = %d\n", i); */
	init_nbr_process(&nbr_man, slow_genus->keys[current], p, i);
	while ((!(has_ended(&nbr_man))) && fmpq_cmp(acc_mass, mass)) {
	  nbr = build_nb(&nbr_man);
#else
	nbr_data_init(nbr_man, slow_genus->keys[current], p, 1);
	while ((!(nbr_data_has_ended(nbr_man))) && fmpq_cmp(acc_mass, mass)) {
	  nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
	  matrix_TYP_init_set_fmpz_mat(&nbr, nbr_fmpz);
#endif // NBR_DATA
	  
#ifdef DEBUG_LEVEL_FULL
	  printf("nbr = \n");
	  print_mat(nbr);
#endif // DEBUG_LEVEL_FULL
	  
	  key_num = -1;
	  isom = NULL;
	  // just something that is not NULL
	  genus_rep = slow_genus->keys[current];
      
	  while ((genus_rep != NULL) && (isom == NULL)) {
#ifdef DEBUG_LEVEL_FULL
	    printf("checking if it is already in the genus...\n");
#endif // DEBUG_LEVEL_FULL
	    genus_rep = hash_table_get_key(slow_genus, nbr, &key_num);
	    if (genus_rep != NULL) {
#ifdef DEBUG_LEVEL_FULL
	      printf("Found candidate :\n");
	      print_mat(genus_rep);
	      printf("Checking for isometry...\n");
#endif // DEBUG_LEVEL_FULL
	      isom = is_isometric(genus_rep, nbr);
	    }
	  }
	  
	  if (isom == NULL) {
#ifdef DEBUG_LEVEL_FULL
	    printf("no Isometry found, adding neighbor...\n");
#endif // DEBUG_LEVEL_FULL
#ifdef NBR_DATA
	    isometry_init_fmpz_mat(s, nbr_isom, p);
#else
	    // !! TODO - should complete here with the isometry for the neighbor
	    isometry_init(s, init_mat(n,n, "1"), 1);
	    
#endif // NBR_DATA
	    assert(isometry_is_isom(s, slow_genus->keys[current], nbr));
	    
	    greedy(nbr, s->s, n, n);
	    s->s_inv = mat_inv(s->s);
	    assert(isometry_is_isom(s, slow_genus->keys[current], nbr));
	    //assert(is_isometry(s, slow_genus->keys[current], nbr, p));
	    // The genus rep isometries were initialized only to contain the
	    // isometry between the parent and its child, we now want to update
	    // these isometries so that they are rational isometries between the
	    // "mother" quadratic form and the genus rep.
	    isometry_init(genus->isoms[slow_genus->num_stored], s->s, 1);
	    isometry_mul(genus->isoms[slow_genus->num_stored], genus->isoms[current],
			 genus->isoms[slow_genus->num_stored]);
	  
	    assert(isometry_is_isom(genus->isoms[slow_genus->num_stored], Q, nbr));
	    
	    isometry_clear(s);
	    
	    hash_table_add(slow_genus, nbr);
	    aut_grp = automorphism_group(nbr);
	    fmpq_set_si(mass_form, 1, aut_grp->order);
	    fmpq_add(acc_mass, acc_mass, mass_form);
#ifdef DEBUG_LEVEL_FULL
	    printf("acc_mass = ");
	    fmpq_print(acc_mass);
	    printf("\n");
#endif // DEBUG_LEVEL_FULL
	  }
	  
#ifndef NBR_DATA
	  advance_nbr_process(&nbr_man);
#else
	  nbr_data_get_next_neighbor(nbr_man);
#endif // NBR_DATA
	}
	
#ifndef NBR_DATA
	i++;
	free_nbr_process(&nbr_man);
	}
#else
	nbr_data_clear(nbr_man);
#endif // NBR_DATA
	current++;
      }
     
    }

    hash_table_recalibrate(genus->genus_reps, slow_genus);

    hash_table_clear(slow_genus);
    
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
    genus->num_auts = (slong**)malloc((genus->num_conductors)*sizeof(slong*));

    assert(genus->genus_reps->num_stored > 0);

    ignore = (bool*)malloc((genus->num_conductors)*sizeof(bool));
    for (c = 0; c < genus->num_conductors; c++)
      genus->num_auts[c] = (slong*)malloc((genus->genus_reps->num_stored)*sizeof(slong));
    
    for (genus_idx = 0; genus_idx < genus->genus_reps->num_stored; genus_idx++) {
      // Determine which subspaces this representative contributes.
      aut_grp = automorphism_group(genus->genus_reps->keys[genus_idx]);
      
      for (c = 0; c < genus->num_conductors; c++)
	ignore[c] = false;

      for (gen_idx = 0; gen_idx < aut_grp->gen_no; gen_idx++) {
	assert(is_isometry(aut_grp->gen[gen_idx],
			   genus->genus_reps->keys[genus_idx],
			   genus->genus_reps->keys[genus_idx], 1));
	isometry_inv(s, genus->isoms[genus_idx]);
	isometry_mul_mat_left(s, aut_grp->gen[gen_idx], s);
	isometry_mul(s, genus->isoms[genus_idx], s);
	
	assert(isometry_is_isom(s, Q, Q));
	
	vals = spinor_norm_isom(genus->spinor, s);
	// !! TODO - we can break the loop after we find one, right?
	for (c = 0; c < genus->num_conductors; c++)
	  if (!ignore[c] && (popcnt(vals & c) & 1))
	    ignore[c] = true;
	isometry_clear(s);
	
      }

      for (c = 0; c < genus->num_conductors; c++) {
	if (!ignore[c]) {
	  genus->lut_positions[c][genus_idx] = genus->dims[c];
	  genus->num_auts[c][genus->dims[c]] = aut_grp->order;
	}
	genus->dims[c] += (ignore[c] ? 0 : 1);
      }
    }

#ifdef DEBUG
    printf("The possible conductors are: \n");
    for (c = 0; c < genus->num_conductors; c++)
      printf("%ld ", genus->conductors[c]);
    printf("\n");
    printf("The corresponding dimensions are: ");
    for (c = 0; c < genus->num_conductors; c++)
      printf("%ld ", genus->dims[c]);
    printf("\n");
#endif // DEBUG
    
    // should also deallocate aut_grp somehow - not clear how to do that
    free(ignore);
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

void genus_clear(genus_t genus)
{
  slong c;

  fmpz_clear(genus->disc);
  
  for (c = 0; c < genus->num_conductors; c++)
    free(genus->num_auts[c]);
  for (c = 0; c < genus->genus_reps->num_stored; c++)
    isometry_clear(genus->isoms[c]);
  free(genus->isoms);
 
  free(genus->num_auts);
  for (c = 0; c < genus->num_conductors; c++)
    free(genus->lut_positions[c]);
  free(genus->lut_positions);
  free(genus->dims);
  free(genus->conductors);
  spinor_clear(genus->spinor);
  hash_table_clear(genus->genus_reps);
  return;
}
