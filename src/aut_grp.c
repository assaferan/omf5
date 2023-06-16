/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : aut_grp.c
 *
 * Description: functions handling the automorphism group
 *              of a lattice
 *
 **********************************************************
 */

// Required packages dependencies

#include <carat/datei.h>

// Self dependencies

#include "aut_grp.h"
#include "matrix_tools.h"
#include "square_matrix.h"

// implementations

// initialize and set an automorphism group from the equivalent type in CARAT
// copies the generators, their number and the order of the group

void aut_grp_init_set_bravais_TYP(aut_grp_t grp, const bravais_TYP* brav)
{
  slong i;
  
  grp->gens = (square_matrix_t*)malloc(brav->gen_no * sizeof(square_matrix_t));

  for (i = 0; i < brav->gen_no; i++)
    square_matrix_set_matrix_TYP(grp->gens[i], brav->gen[i]);

  grp->num_gens = brav->gen_no;
  grp->order = brav->order;
  
  return;
}

// clearing the memory allocated for the group

void aut_grp_clear(aut_grp_t grp)
{
  slong i;
  
  for (i = 0; i < grp->num_gens; i++)
    square_matrix_clear(grp->gens[i]);

  free(grp->gens);

  return;
}

// initialize the automorphism group from a gram matrix

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

// check whether a square matrix is inside a list
// !! TODO - move it to the files for square_matrix, and make it public

bool square_matrix_in_list(const square_matrix_t mat, const square_matrix_t* mat_list, slong len)
{
  slong i;

  for (i = 0; i < len; i++)
    if (square_matrix_is_equal(mat, mat_list[i]))
      return true;

  return false;
}

// get all elements of an automorphism group,
// and write them into an array of square matrices.

void aut_grp_get_elements(square_matrix_t* elts, const aut_grp_t grp)
{
  slong next_idx, elt_idx, gen_idx;
  square_matrix_t cand;

  for (elt_idx = 0; elt_idx < grp->order; elt_idx++)
    square_matrix_init(elts[elt_idx]);
  square_matrix_init(cand);

  // starting with 1
  square_matrix_one(elts[0]);
 
  next_idx = 1;
  for (elt_idx = 0; (elt_idx < next_idx) && (next_idx < grp->order); elt_idx++) {
    // Will it be faster with taking powers?
    for (gen_idx = 0; gen_idx < grp->num_gens; gen_idx++) {
      square_matrix_mul(cand, elts[elt_idx], grp->gens[gen_idx]);
      if (!square_matrix_in_list(cand, elts, next_idx)) {
	square_matrix_set(elts[next_idx], cand);
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
