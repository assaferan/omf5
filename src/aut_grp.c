#include <carat/datei.h>

#include "aut_grp.h"
#include "matrix_tools.h"
#include "square_matrix.h"

void aut_grp_init_set_bravais_TYP(aut_grp_t grp, const bravais_TYP* brav)
{
  slong i;

  grp->det_gens = (short *) malloc(brav->gen_no * sizeof(short));
  
  grp->gens = (square_matrix_t*)malloc(brav->gen_no * sizeof(square_matrix_t));

  for (i = 0; i < brav->gen_no; i++) {
    square_matrix_set_matrix_TYP(grp->gens[i], brav->gen[i]);
    // !! TODO - check, maybe need to take first prime not dividing the denominator ?
    grp->det_gens[i] = ((p_mat_det(brav->gen[i], 3) == 1) ? 1 : -1);
  }

  grp->num_gens = brav->gen_no;
  grp->order = brav->order;
  
  return;
}

void aut_grp_clear(aut_grp_t grp)
{
  slong i;
  
  for (i = 0; i < grp->num_gens; i++)
    square_matrix_clear(grp->gens[i]);

  free(grp->gens);
  free(grp->det_gens);

  return;
}

void aut_grp_init_square_matrix(aut_grp_t grp, const square_matrix_t mat)
{
  bravais_TYP* brav;
  matrix_TYP* Q;

  Q = matrix_TYP_init_set_square_matrix(mat);
  brav = automorphism_group(Q);

  aut_grp_init_set_bravais_TYP(grp, brav);

  free_bravais(brav);
  free_mat(Q);

  return;
}

bool square_matrix_in_list(const square_matrix_t mat, const square_matrix_t* mat_list, slong len)
{
  slong i;

  for (i = 0; i < len; i++)
    if (square_matrix_is_equal(mat, mat_list[i]))
      return true;

  return false;
}

void aut_grp_get_elements_and_dets(square_matrix_t* elts, short* dets, const aut_grp_t grp)
{
  slong next_idx, elt_idx, gen_idx;
  square_matrix_t cand;
  
  for (elt_idx = 0; elt_idx < grp->order; elt_idx++)
    square_matrix_init(elts[elt_idx]);
  square_matrix_init(cand);

  // starting with 1
  square_matrix_one(elts[0]);
  dets[0] = 1;
 
  next_idx = 1;
  for (elt_idx = 0; (elt_idx < next_idx) && (next_idx < grp->order); elt_idx++) {
    // Will it be faster with taking powers?
    for (gen_idx = 0; gen_idx < grp->num_gens; gen_idx++) {
      square_matrix_mul(cand, elts[elt_idx], grp->gens[gen_idx]);
      if (!square_matrix_in_list(cand, elts, next_idx)) {
	square_matrix_set(elts[next_idx], cand);
	dets[next_idx] = grp->det_gens[gen_idx] * dets[elt_idx];
	next_idx++;
      }
    }
  }

  for (elt_idx = 0; elt_idx < grp->order; elt_idx++) {
    square_matrix_set(cand, elts[elt_idx]);
    square_matrix_transpose(elts[elt_idx], cand);
  }
  
  square_matrix_clear(cand);
  
  return; 
}

void aut_grp_get_elements(square_matrix_t* elts, const aut_grp_t grp)
{
  short* dets;

  dets = (short*)malloc(grp->order * sizeof(short));

  aut_grp_get_elements_and_dets(elts, dets, grp);

  free(dets);
  
  return;
}

// This is not very efficient, but at the moment will suffice
void aut_grp_get_special_elements(square_matrix_t* elts, const aut_grp_t grp)
{
  square_matrix_t* all_elts;
  short* dets;
  slong src_idx, dest_idx;

  all_elts = (square_matrix_t*)malloc(grp->order * sizeof(square_matrix_t));
  dets = (short*)malloc(grp->order * sizeof(short));

  aut_grp_get_elements_and_dets(all_elts, dets, grp);

  dest_idx = 0;
  for (src_idx = 0; src_idx < grp->order; src_idx++) {
    if (dets[src_idx] == 1)
      square_matrix_set(elts[dest_idx++], all_elts[src_idx]);
  }
  
  free(dets);
  free(all_elts);
  
  return;
}
