#include <carat/datei.h>

#include "aut_grp.h"
#include "matrix_tools.h"
#include "square_matrix.h"

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

void aut_grp_clear(aut_grp_t grp)
{
  slong i;
  
  for (i = 0; i < grp->num_gens; i++)
    square_matrix_clear(grp->gens[i]);

  free(grp->gens);

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
