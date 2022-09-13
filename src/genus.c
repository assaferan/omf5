#include <carat/matrix.h>

#include "genus.h"
#include "hash.h"
#include "mass.h"
#include "matrix_tools.h"
#include "neighbor.h"

/* compute the genus of a quadratic form */
hash_table* get_genus_reps(matrix_TYP* Q)
{
  bravais_TYP *aut_grp;
  matrix_TYP *nbr, *isom, *genus_rep, *s;
  fmpq_t mass, acc_mass, mass_form;
  fmpz_t prime;
  int i, p, current, next_idx, key_num;
  size_t genus_size;
  fmpz_t genus_size_fmpz;
  hash_table* genus;
  neighbor_manager nbr_man;

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

  genus = create_hash(genus_size);
  
  add(genus, Q);

  current = 0;
  next_idx = 1;
  
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
#ifdef DEBUG
    if (fmpq_cmp(acc_mass, mass) > 0) {
      printf("Error! Accumulated too much mass!\n");
      return genus;
    }
#endif // DEBUG
    /* !! TODO - we don't really need to restrict to good primes here, */
    /* but let's check these first */
    do {
      fmpz_nextprime(prime, prime, TRUE);
      p = fmpz_get_ui(prime);
#ifdef DEBUG
      printf("p = %d, Q = ", p);
      print_mat(Q);
#endif //DEBUG
    }
    while (!p_mat_det(Q, p));

    while ((current < genus->num_stored) && fmpq_cmp(acc_mass, mass)){
#ifdef DEBUG_LEVEL_FULL
      printf("current = %d\n", current);
#endif // DEBUG_LEVEL_FULL
      i = 0;
      
      /* Right now all the isotropic vectors are paritioned to p */
      /* sets, by index as Gonzalo did */
      
      while ((i < p) && fmpq_cmp(acc_mass, mass)) {
	/* printf("i = %d\n", i); */
	init_nbr_process(&nbr_man, genus->keys[current], p, i);
	while ((!(has_ended(&nbr_man))) && fmpq_cmp(acc_mass, mass)) {
	  nbr = build_nb(&nbr_man);
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
	    fmpq_set_si(mass_form, 1, aut_grp->order);
	    fmpq_add(acc_mass, acc_mass, mass_form);
#ifdef DEBUG_LEVEL_FULL
	    printf("acc_mass = ");
	    fmpq_print(acc_mass);
	    printf("\n");
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

  fmpq_clear(mass);
  fmpq_clear(acc_mass);
  fmpq_clear(mass_form);
  fmpz_clear(prime);
  
  return genus;
}

