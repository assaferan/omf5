#ifndef __AUT_GRP_H__
#define __AUT_GRP_H__

#include <carat/typedef.h>

#include "square_matrix.h"

typedef struct {
  square_matrix_t* gens;
  slong num_gens;
  slong order;
} aut_grp_struct;

typedef aut_grp_struct aut_grp_t[1];

void aut_grp_init_set_bravais_TYP(aut_grp_t grp, const bravais_TYP* brav);

void aut_grp_init_square_matrix(aut_grp_t grp, const square_matrix_t mat);

void aut_grp_get_elements(square_matrix_t* elts, const aut_grp_t grp);

void aut_grp_clear(aut_grp_t grp);

#endif // __AUT_GRP_H__
