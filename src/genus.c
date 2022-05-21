#include <gmp.h>

#include "carat/autgrp.h"
#include "carat/symm.h"

#include "arith.h"
#include "genus.h"
#include "hash.h"
#include "matrix_tools.h"
#include "neighbor.h"

typedef mpz_t Z;

/* initialize the neighbors data */
nbrs_data* q61_init(int p, int k)
{
  nbrs_data* dtm;
  
  int s, i;

  matrix_TYP* genus[8];

  /* Set of representatives for the quadratic forms */
  
  int genus_coeffs[8][15] =
    {
     {2,1,2,1,0,2,0,1,0,16,0,0,0,0,2},
     {2,1,2,0,0,2,1,1,1,6,-1,-1,0,-3,6},
     {2,1,2,0,0,2,0,0,0,4,1,0,0,-1,6},
     {2,-1,2,-1,1,2,1,-1,-1,6,-1,1,0,2,8},
     {2,1,2,1,0,4,-1,-1,0,4,0,0,1,0,4},
     {2,1,2,-1,0,2,0,0,0,6,-1,0,1,-1,6},
     {2,-1,2,0,0,2,1,0,-1,4,1,0,0,0,8},
     {2,0,2,0,0,2,0,0,1,4,0,1,1,2,6}
    };

  /* printf("in q61_init. Initializing genus representatives to:\n"); */
  
  for (i = 0; i < 8; i++) {
    genus[i] = init_sym_matrix(genus_coeffs[i]);
    /* print_mat(genus[i]); */
  }
  
  dtm = (nbrs_data*) malloc(sizeof(nbrs_data));

  /* Hardcode detect class */
  // !! This depends on the size of our hash
  // In general replace by some hash
  dtm->th61 = (int*) malloc(20*sizeof(int));

  /* printf("Computing theta series...\n"); */
  for (i = 0; i < 8; i++) {
    s = hash_form(genus[i]);

    /* printf("s = %d\n", s); */
    dtm->th61[s] = i;
  }

  dtm->Q = init_sym_matrix(genus_coeffs[k]);

  dtm->v = get_isotropic_vector(dtm->Q, p);

  return dtm;
}

/*
int hash_form(matrix_TYP* Q)
{
  int norm, x;
  int num_short[3];
  
  for (norm = 2; norm <= 6; norm += 2) {
    short_vectors(Q, norm, norm, 0, 1, &(num_short[norm/2 - 1]));
  }

  x = num_short[0] - num_short[1] + num_short[2] ;

  return x;
}
*/

/* identify the genus representative of Q (returns the index) */
int q61_id(matrix_TYP* Q, int* th61)
{
  return th61[hash_form(Q)];
}

/* compute the genus of a quadratic form */
hash_table* get_genus_reps(matrix_TYP* Q)
{
  bravais_TYP *aut_grp;
  matrix_TYP *iso_vec, *nbr, *isom, *genus_rep;
  rational mass, acc_mass, mass_form;
  Z prime;
  int i, p, current, next_idx, key_num;
  int options[6] = {0};
  size_t genus_size;
  hash_table* genus;

  /* until we implement the mass formula, have it fixed */
  
  mass.z = 31;
  mass.n = 96;

  printf("mass = %d / %d \n", mass.z, mass.n);

  /* this is ceiling */
  genus_size = (mass.z + mass.n - 1)/ mass.n;

  genus_size = (4 * genus_size) / 3; // load factor

  genus = create_hash(genus_size);
  
  add(genus, Q);

  current = 0;
  next_idx = 1;
  
  aut_grp = autgrp(&Q, 1, NULL, NULL, 0, options);
  
  acc_mass.z = 1;
  acc_mass.n = aut_grp->order;

  printf("acc_mass = %d / %d \n", acc_mass.z, acc_mass.n);

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

    while (current < genus->num_stored) {
      iso_vec = get_isotropic_vector(genus->keys[current], p);

      i = 0;
      nbr = q61_nb(genus->keys[current], p, iso_vec);

      key_num = -1;
      isom = NULL;
      
      while ((genus_rep != NULL) && (isom == NULL)) {
	genus_rep = get_key(genus, nbr, &key_num);
	if (genus_rep != NULL)
	  isom = isometry(&genus_rep, &nbr, 1,
			  NULL, NULL, NULL, 0, options);
      }

      if (isom == NULL) {
	add(genus, nbr);
	aut_grp = autgrp(&nbr, 1, NULL, NULL, 0, options);
	mass_form.z = 1;
	mass_form.n = aut_grp->order;
	acc_mass = rational_sum(acc_mass, mass_form);
      }
      
      /* Right now all the isotropic vectors are paritioned to p */
      /* sets, by index as Gonzalo did */

      current++;
    }
    
    
  }

  mpz_clear(prime);
  
  return genus;
}
