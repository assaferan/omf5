/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : aut_grp.h
 *
 * Description: automorphism groups of lattices
 *
 **********************************************************
 */

#ifndef __AUT_GRP_H__
#define __AUT_GRP_H__

// Required packages dependencies

#include <carat/typedef.h>

// Self dependencies

#include "square_matrix.h"

/*************************************************************
 *
 * Type: aut_grp_t
 *
 * Description: automorphism group of a lattice.
 *
 * Fields:
 *     + gens (square_matrix_t*) - a list of generators
 *     + num_gens (slong)        - number of generators
 *     + order (slong)           - order of group
 *
 * Used to hold the data of the automorphism group.
 ************************************************************
 */

typedef struct {
  square_matrix_t* gens;
  slong num_gens;
  slong order;
} aut_grp_struct;

typedef aut_grp_struct aut_grp_t[1];

/***************************************************************
 *
 * Function: aut_grp_init_set_bravais_TYP
 *
 * Description: Conversion from bravais_TYP (from CARAT).
 *
 * Arguments:
 *     + brav (const bravais_TYP*) - CARAT's data structure
 *
 * Returns:
 *     + grp (aut_grp_t) - the automorphism group
 *
 **************************************************************
 */

void aut_grp_init_set_bravais_TYP(aut_grp_t grp, const bravais_TYP* brav);

/********************************************************************
 *
 * Function: aut_grp_init_square_matrix
 *
 * Description: Constructs the automorphism group of the lattice
 *              with given gram matrix.
 *
 * Arguments:
 *     + gram (const square_matrix_t) - gram matrix of the lattice
 *
 * Returns:
 *     + grp (aut_grp_t) - the automorphism group
 *
 *******************************************************************
 */

void aut_grp_init_square_matrix(aut_grp_t grp, const square_matrix_t gram);

/********************************************************************
 *
 * Function: aut_grp_get_elements
 *
 * Description: Construct all elements of the group.
 *
 * Arguments:
 *     + grp (const aut_grp_t) - the automorphism group
 *
 * Returns:
 *     + elts (square_matrix_t*) - a list of matrices, which form
 *                                 all the elements of the group
 *
 *******************************************************************
 */

void aut_grp_get_elements(square_matrix_t* elts, const aut_grp_t grp);

/***********************************************************************
 *
 * Function: aut_grp_clear
 *
 * Description: Clear all memory allocated for the automorphism group.
 *
 * Arguments:
 *     +  grp (aut_grp_t) - the automorphism group
 *
 **********************************************************************
 */

void aut_grp_clear(aut_grp_t grp);

#endif // __AUT_GRP_H__
