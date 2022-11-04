#include <assert.h>

#include "genus.h"
#include "hash.h"
#include "isometry.h"
#include "nbr_data.h"
#include "neighbor.h"
#include "square_matrix.h"
#include "weight.h"

// Just update the matrix acting at each index
int process_isotropic_vector_isometries(neighbor_manager_t nbr_man,
					square_matrix_t** large_hecke,
					const genus_t genus,
					double* theta_time, double* isom_time,
					double* total_time, int* num_isom, int gen_idx)
{
  int i, j;
  clock_t cputime;
  square_matrix_t nbr, aut_tr, orbit_sum;
  isometry_t s_nbr, s_inv;
  isometry_t hash_isom;
  bool non_singular;
  int vec_cmp;
  slong orb_size, stab_size;
  isometry_t* orbit_isom;
  vector_t* orbit;
  bool found;
  vector_t g_vec;
  
  assert(genus->genus_reps->num_stored > 0);

  orbit = (vector_t*)malloc(nbr_man->num_auts * sizeof(vector_t));
  orbit_isom = (isometry_t*)malloc(nbr_man->num_auts * sizeof(isometry_t));
  isometry_init(orbit_isom[0]);

  // right now we grow the orbit.
  // Could also compute on the fly and in the end divide by the stabilizer's size
  // weird that it isn't already reduced?
  vector_mod_p(nbr_man->iso_vec, nbr_man->p);
  // this might be unnecessary, but we want to be on the safe side for now
  normalize_mod_p(nbr_man->iso_vec, nbr_man->p);
  vector_set(orbit[0], nbr_man->iso_vec);
  // Applying the automorphism group to the isotropic vector
  stab_size = 0;
  orb_size = 1;
  for (i = 0; i < nbr_man->num_auts; i++) {
    square_matrix_mul_vec_left(g_vec, nbr_man->iso_vec, nbr_man->auts[i]);
    vector_mod_p(g_vec, nbr_man->p);
    normalize_mod_p(g_vec, nbr_man->p);
    vec_cmp = vector_cmp(g_vec, nbr_man->iso_vec);
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
	if (vector_cmp(g_vec, orbit[j]) == 0) {
	  found = true;
	  break;
	}
      if (!found) {
	vector_set(orbit[orb_size], g_vec);
	square_matrix_transpose(aut_tr, nbr_man->auts[i]);
	isometry_init_set_square_matrix(orbit_isom[orb_size], aut_tr, 1);
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
  
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
  
  cputime = clock();
  i = hash_table_index_and_isom(genus->genus_reps, nbr, hash_isom, theta_time, isom_time, num_isom);
  (*total_time) += clock() - cputime;

  assert( (i >= 0) && (i <= genus->dims[0]) );

  // TODO - compute only the upper half, and complete using hermitian property
  // Determine which subspaces this representative contributes.

  // TODO - reduce the number of actual computations required here

  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[gen_idx], nbr));
  for (j = 0; j < orb_size; j++)
    assert(isometry_is_isom(orbit_isom[j], genus->genus_reps->keys[gen_idx], genus->genus_reps->keys[gen_idx]));
  assert(isometry_is_isom(genus->isoms[gen_idx], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[gen_idx]));
  
  isometry_muleq_left(s_nbr, genus->isoms[gen_idx]);
  for (j = 0; j < orb_size; j++)
     isometry_muleq_left(orbit_isom[j], genus->isoms[gen_idx]);
  
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], nbr));
  for (j = 0; j < orb_size; j++)
     assert(isometry_is_isom(orbit_isom[j], genus->genus_reps->keys[0], genus->genus_reps->keys[gen_idx]));

  isometry_inv(s_inv, genus->isoms[gen_idx]);
  for (j = 0; j < orb_size; j++)
     isometry_muleq_right(orbit_isom[j], s_inv);
  
  for (j = 0; j < orb_size; j++)
    assert(isometry_is_isom(orbit_isom[j], genus->genus_reps->keys[0],
			    genus->genus_reps->keys[0]));
  
  assert(isometry_is_isom(hash_isom, nbr, genus->genus_reps->keys[i]));
  
  isometry_muleq_right(s_nbr, hash_isom);
  
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0], genus->genus_reps->keys[i]));
  assert(isometry_is_isom(genus->isoms[i], genus->genus_reps->keys[0],
			  genus->genus_reps->keys[i]));
  
  isometry_inv(s_inv, genus->isoms[i]);
  
  assert(isometry_is_isom(s_inv, genus->genus_reps->keys[i],
			  genus->genus_reps->keys[0]));
  
  isometry_muleq_right(s_nbr, s_inv);
  
  assert(isometry_is_isom(s_nbr, genus->genus_reps->keys[0],
			  genus->genus_reps->keys[0]));

  // !! TODO - handle denominators !!!
  square_matrix_init(orbit_sum);
  square_matrix_zero(orbit_sum);
  for (j = 0; j < orb_size; j++) {
    square_matrix_add(orbit_sum, orbit_sum, orbit_isom[j]->s);
  }
  square_matrix_muleq_right(orbit_sum, s_nbr->s);
  square_matrix_add(large_hecke[gen_idx][i],large_hecke[gen_idx][i], orbit_sum); 

  isometry_clear(s_nbr);
  isometry_clear(s_inv);
  isometry_clear(hash_isom);

  square_matrix_clear(nbr);

  free(orbit);
  free(orbit_isom);
  
  return orb_size;
}


slong hecke_col_isometries(square_matrix_t** large_hecke,
			   int p, int gen_idx, const genus_t genus)
{
  square_matrix_t Q;
  neighbor_manager_t nbr_man;
  int lc, i;
  int num_isom;
  double theta_time, isom_time, total_time;

  num_isom = lc = 0;
  theta_time = isom_time = total_time = 0;

  lc = 0;
  
  square_matrix_set(Q,genus->genus_reps->keys[gen_idx]);

#ifdef DEBUG_LEVEL_FULL
  printf("initialized Q: \n");
  square_matrix_print(Q);
#endif // DEBUG_LEVEL_FULL

  for (i = 0; i < p; i++) {
    nbr_process_init(nbr_man, Q, p, i);

    while (!(nbr_process_has_ended(nbr_man))) {
      lc += process_isotropic_vector_isometries(nbr_man, large_hecke, genus,
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


square_matrix_t** hecke_matrices_isometries(const genus_t genus, int p)
{
  square_matrix_t** hecke;
 
  int gen_idx;
  slong i, j;

  hecke = (square_matrix_t**)malloc(genus->dims[0] * sizeof(square_matrix_t*));

  for (i = 0; i < genus->dims[0]; i++) {
    hecke[i] = (square_matrix_t*)malloc(genus->dims[0] * sizeof(square_matrix_t));
    for (j = 0; j < genus->dims[0]; j++)
      square_matrix_zero(hecke[i][j]);
  }
  
  for (gen_idx = 0; gen_idx < genus->dims[0]; gen_idx++)
    hecke_col_isometries(hecke, p, gen_idx, genus);
  
  return hecke;
}
