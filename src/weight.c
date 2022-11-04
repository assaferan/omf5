#include <assert.h>

#include <flint/fmpq_mat.h>

#include "genus.h"
#include "hash.h"
#include "isometry.h"
#include "matrix_tools.h"
#include "nbr_data.h"
#include "neighbor.h"
#include "square_matrix.h"
#include "weight.h"

// Just update the matrix acting at each index
int process_isotropic_vector_isometries(neighbor_manager_t nbr_man,
					isometry_t** large_hecke,
					const genus_t genus,
					double* theta_time, double* isom_time,
					double* total_time, int* num_isom, int gen_idx)
{
  int i, j;
  clock_t cputime;
  square_matrix_t nbr, aut_tr;
  isometry_t s_nbr, s_inv;
  isometry_t hash_isom;
  bool non_singular;
  int vec_cmp;
  slong orb_size, stab_size;
  isometry_t* orbit_isom;
  isometry_t orbit_sum;
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

  isometry_zero(orbit_sum);
  for (j = 0; j < orb_size; j++) {
    isometry_add(orbit_sum, orbit_sum, orbit_isom[j]);
  }
  isometry_muleq_right(orbit_sum, s_nbr);
  isometry_add(large_hecke[gen_idx][i],large_hecke[gen_idx][i], orbit_sum);

  isometry_clear(s_nbr);
  isometry_clear(s_inv);
  isometry_clear(hash_isom);

  square_matrix_clear(nbr);

  free(orbit);
  free(orbit_isom);
  
  return orb_size;
}


slong hecke_col_isometries(isometry_t** large_hecke,
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


isometry_t** hecke_matrices_isometries(const genus_t genus, int p)
{
  isometry_t** hecke;
 
  int gen_idx;
  slong i, j;

  hecke = (isometry_t**)malloc(genus->dims[0] * sizeof(isometry_t*));

  for (i = 0; i < genus->dims[0]; i++) {
    hecke[i] = (isometry_t*)malloc(genus->dims[0] * sizeof(isometry_t));
    for (j = 0; j < genus->dims[0]; j++)
      isometry_zero(hecke[i][j]);
  }
  
  for (gen_idx = 0; gen_idx < genus->dims[0]; gen_idx++)
    hecke_col_isometries(hecke, p, gen_idx, genus);
  
  return hecke;
}

// right now just for the standard representation
void invariant_subspace(fmpq_mat_t sub, const aut_grp_t grp)
{
  slong i, j;
  fmpq_mat_t new_sub, gen1_Q, gen2_Q, one, ker;
  fmpz_mat_t gen1_Z, gen2_Z;

  fmpq_mat_init(one, QF_RANK, QF_RANK);
  fmpq_mat_one(one);
  fmpq_mat_init(sub, QF_RANK, QF_RANK);
  fmpq_mat_one(sub);

  fmpq_mat_init(gen1_Q, QF_RANK, QF_RANK);
  fmpq_mat_init(gen2_Q, QF_RANK, QF_RANK);
  
  for (i = 0; i < grp->num_gens; i++) {
    fmpz_mat_init_set_square_matrix(gen1_Z, grp->gens[i]);
    fmpq_mat_set_fmpz_mat(gen1_Q, gen1_Z);
    if (grp->det_gens[i] == 1) {
      fmpq_mat_sub(gen1_Q, gen1_Q, one);
      fmpq_mat_kernel(ker, gen1_Q);
      fmpq_mat_meet(new_sub, sub, ker);
      fmpq_mat_clear(sub);
      fmpq_mat_clear(ker);
      fmpq_mat_init_set(sub, new_sub);
      fmpq_mat_clear(new_sub);
    }
    else {
       for (j = 0; j < grp->num_gens; j++) {
	 fmpz_mat_init_set_square_matrix(gen2_Z, grp->gens[j]);
	 fmpq_mat_set_fmpz_mat(gen2_Q, gen2_Z);
	 fmpq_mat_mul(gen2_Q, gen1_Q, gen2_Q);
	 if (grp->det_gens[j] == -1) {
	   fmpq_mat_sub(gen2_Q, gen2_Q, one);
	   fmpq_mat_kernel(ker, gen2_Q);
	   fmpq_mat_meet(new_sub, sub, ker);
	   fmpq_mat_clear(sub);
	   fmpq_mat_clear(ker);
	   fmpq_mat_init_set(sub, new_sub);
	   fmpq_mat_clear(new_sub);
	 }
	 fmpz_mat_clear(gen2_Z);
	 
       }
    }
    fmpz_mat_clear(gen1_Z);
    
  }

  fmpq_mat_clear(gen1_Q);
  fmpq_mat_clear(gen2_Q);
  fmpq_mat_clear(one);

  return;
}

matrix_TYP* hecke_matrices_std_rep(const genus_t genus, int p)
{
  slong i,j,ii,jj;
  isometry_t** large_hecke;
  matrix_TYP* hecke;
  fmpq_mat_t* inv_subs;
  aut_grp_t auts;
  fmpq_mat_t hecke_tr, hecke_res, hecke_block_Q;
  fmpz_mat_t hecke_block_Z;
  slong dim, dim1, dim2;

  large_hecke = hecke_matrices_isometries(genus, p);

  printf("Invariant subspaces are: \n");
  inv_subs = (fmpq_mat_t*)malloc((genus->dims[0]) * sizeof(fmpq_mat_t));
  for (i = 0; i < genus->dims[0]; i++) {
    aut_grp_init_square_matrix(auts, genus->genus_reps->keys[i]);
    invariant_subspace(inv_subs[i], auts);
    aut_grp_clear(auts);
    fmpq_mat_print(inv_subs[i]);
    printf(",");
  }
  printf("Base change isometries are: \n");
  // print the base change isometries
  for (i = 0; i < genus->dims[0]; i++) {
    isometry_print(genus->isoms[i]);
    // transfer the invariant subspace to the original space
    fmpq_mat_muleq(inv_subs[i], inv_subs[i], genus->isoms[i]);
  }
  printf("Invariant subspaces are: \n");
  for (i = 0; i < genus->dims[0]; i++) {
    fmpq_mat_print(inv_subs[i]);
    printf(",");
  }

  dim = 0;
  for (i = 0; i < genus->dims[0]; i++) {
    dim += fmpq_mat_nrows(inv_subs[i]);
  }
  printf("dim = %ld\n", dim);

  hecke = init_mat(dim, dim, "");

  dim1 = 0;
  for (i = 0; i < genus->dims[0]; i++) {
    dim2 = 0;
    for (j = 0; j < genus->dims[0]; j++) {
      fmpq_mat_init(hecke_tr, QF_RANK, QF_RANK);
      fmpz_mat_init_set_square_matrix(hecke_block_Z, large_hecke[i][j]->s);
      fmpq_mat_init_set_fmpz_mat(hecke_block_Q, hecke_block_Z);
      fmpq_mat_transpose(hecke_tr, hecke_block_Q);
      fmpq_mat_scalar_div(hecke_tr, hecke_tr, large_hecke[i][j]->denom);
      fmpq_mat_init(res, fmpq_mat_nrows(inv_subs[i]), fmpq_mat_nrows(inv_subs[i]));
      restrict_mat(res, hecke_tr, inv_subs[i]);
      for (ii = 0; ii < fmpq_mat_nrows(res); ii++)
	for (jj = 0; jj < fmpq_mat_ncols(res); jj++)
	  hecke[dim1 + ii][dim2 + jj] = fmpq_get_si(fmpq_mat_entry(res, ii, jj));
      fmpq_mat_clear(res);
      fmpq_mat_clear(hecke_tr);
      fmpq_mat_clear(hecke_block_Q);
      fmpz_mat_clear(hecke_block_Z);
      dim2 += fmpq_mat_nrows(inv_subs[i]);
    }
    dim1 += fmpq_mat_nrows(inv_subs[i]);
  }
  
  for (i = 0; i < genus->dims[0]; i++) {
    fmpq_mat_clear(inv_subs[i]);
  }
  free(inv_subs);

  return hecke;
}
