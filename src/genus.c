#include <carat/matrix.h>

#include <flint/fmpq.h>
#include <flint/fmpz.h>

#include "genus.h"
#include "hash.h"
#include "mass.h"
#include "matrix_tools.h"
#ifndef NBR_DATA
#include "neighbor.h"
#else
#include "nbr_data.h"
#endif // NBR_DATA
#include "typedefs.h"

/* compute the genus of a quadratic form */
void genus_init(genus_t genus, matrix_TYP* Q)
{
  bravais_TYP *aut_grp;
  matrix_TYP *nbr, *isom, *genus_rep, *s;
  fmpq_t mass, acc_mass, mass_form;
  fmpz_t prime;
  int p, current, key_num;
  size_t genus_size;
  fmpz_t genus_size_fmpz;
  hash_table_t slow_genus;

#ifndef NBR_DATA
  neighbor_manager nbr_man;
  int i;
#else
  nbr_data_t nbr_man;
  fmpz_mat_t nbr_isom, nbr_fmpz;

  fmpz_mat_init(nbr_isom, N, N);
  fmpz_mat_init(nbr_fmpz, N, N);
#endif // NBR_DATA
  
  /* until we implement the mass formula, have it fixed */

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
      fmpz_nextprime(prime, prime, TRUE);
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
	    s = init_mat(5,5,"1");
	    greedy(nbr, s, 5, 5);
	    free_mat(s);
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
    
  fmpq_clear(mass);
  fmpq_clear(acc_mass);
  fmpq_clear(mass_form);
  fmpz_clear(prime);

#ifdef NBR_DATA
  fmpz_mat_clear(nbr_fmpz);
  fmpz_mat_clear(nbr_isom);
#endif // NBR_DATA

  hash_table_recalibrate(genus->genus_reps, slow_genus);
  
  return;
}

void genus_clear(genus_t genus)
{
  hash_table_clear(genus->genus_reps);
  return;
}
