#include <gmp.h>

#include "arith.h"
#include "genus.h"
#include "hash.h"
#include "matrix_tools.h"
#include "neighbor.h"

typedef mpz_t Z;

/* compute the genus of a quadratic form */
hash_table* get_genus_reps(matrix_TYP* Q)
{
  bravais_TYP *aut_grp;
  matrix_TYP *nbr, *isom, *genus_rep, *s;
  rational mass, acc_mass, mass_form;
  Z prime;
  int i, p, current, next_idx, key_num;
  //  int options[6] = {0};
  size_t genus_size;
  hash_table* genus;
  neighbor_manager nbr_man;

  /* until we implement the mass formula, have it fixed */
  
  mass.z = 31;
  mass.n = 96;

#ifdef DEBUG_LEVEL_FULL
  printf("mass = %d / %d \n", mass.z, mass.n);
#endif // DEBUG_LEVEL_FULL

  /* this is ceiling */
  genus_size = (mass.z + mass.n - 1)/ mass.n;

  genus_size = (4 * genus_size) / 3; // load factor

  genus = create_hash(genus_size);
  
  add(genus, Q);

  current = 0;
  next_idx = 1;
  
  aut_grp = automorphism_group(Q);
  
  acc_mass.z = 1;
  acc_mass.n = aut_grp->order;

#ifdef DEBUG_LEVEL_FULL
  printf("acc_mass = %d / %d \n", acc_mass.z, acc_mass.n);
#endif // DEBUG_LEVEL_FULL

  mpz_init(prime);
  mpz_set_ui(prime, 1);
  
  while (rational_lt(acc_mass, mass)) {
    /* !! TODO - we don't really need to restrict to good primes here, */
    /* but let's check these first */
    do {
      mpz_nextprime(prime, prime);
      p = mpz_get_ui(prime);
    }
    while (!p_mat_det(Q, p));

    while ((current < genus->num_stored) && rational_lt(acc_mass, mass)){
#ifdef DEBUG_LEVEL_FULL
      printf("current = %d\n", current);
#endif // DEBUG_LEVEL_FULL
      i = 0;
      
      /* Right now all the isotropic vectors are paritioned to p */
      /* sets, by index as Gonzalo did */
      
      while ((i < p) && rational_lt(acc_mass, mass)) {
	/* printf("i = %d\n", i); */
	init_nbr_process(&nbr_man, genus->keys[current], p, i);
	while ((!(has_ended(&nbr_man))) && rational_lt(acc_mass, mass)) {
	  nbr = q61_nb(&nbr_man);
#ifdef DEBUG_LEVEL_FULL
	  printf("nbr = \n");
	  print_mat(nbr);
#endif // DEBUG_LEVEL_FULL
	  
	  key_num = -1;
	  isom = NULL;
	  // just something that is not NULL
	  genus_rep = genus->keys[current];
      
	  while ((genus_rep != NULL) && (isom == NULL)) {
#ifdef DEBUG_LEVEL_FULL
	    printf("checking if it is already in the genus...\n");
#endif // DEBUG_LEVEL_FULL
	    genus_rep = get_key(genus, nbr, &key_num);
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
	    add(genus, nbr);
	    aut_grp = automorphism_group(nbr);
	    mass_form.z = 1;
	    mass_form.n = aut_grp->order;
	    acc_mass = rational_sum(acc_mass, mass_form);
#ifdef DEBUG_LEVEL_FULL
	    printf("acc_mass = %d / %d \n", acc_mass.z, acc_mass.n);
#endif // DEBUG_LEVEL_FULL
	  }
	  
	  advance_nbr_process(&nbr_man);
	}
	  
	i++;
	free_nbr_process(&nbr_man);
      }

      current++;
    }
     
  }

  mpz_clear(prime);
  
  return genus;
}
