/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : genus.h
 *
 * Description: Genus of a lattice.
 *              Used also to record relevant data about
 *              the space of orthogonal modular forms
 *              supported on this genus with different
 *              conductors.
 *
 **********************************************************
 */

#ifndef __GENUS_H__
#define __GENUS_H__

// Required packages dependencies

#include <carat/matrix.h>

// Self dependencies

#include "hash.h"
#include "isometry.h"
#include "spinor.h"
#include "square_matrix.h"

/***********************************************************************
 *
 * Type: genus_t
 *
 * Description: this struct holds data for the spaces of
 *              orthogonal modular forms associated to a
 *              certain genus, with different conductors.
 *
 * Fields:
 *     + genus_reps (hash_table_t) - a hash table for the
 *                                   representatives of the genus,
 *                                   these are represented by
 *                                   gram matrices of the lattices,
 *                                   and the hash table is for
 *                                   faster isometry testing
 *     + spinor (spinor_t) - a data structure used to record
 *                           all data relevant to compute the
 *                           spinor norm for isometries of this
 *                           quadratic space (see spinor.h)
 *     + disc (fmpz_t) - the discriminant of lattices in the genus
 *     + dims (slong*) - dimensions of each of the spaces of
 *                       orthogonal modular forms, by conductor
 *     + conductors (slong*) - list of relevant conductors
 *                             (square free divisors of disc)
 *     + num_auts (slong**) - the number of automorphisms of each
 *                            lattice, spread in a redundant way,
 *                            by conductor
 *     + lut_positions (slong**) - look up tables for every conductor,
 *                                 identifying the lattices in the genus
 *                                 on which the space is supported.
 *     + num_conductors (slong) - number of possible conductors.
 *     + isoms (isometry_t*) - isometries (over QQ) between the
 *                             different genus representatives
 *
 ***********************************************************************
 */

typedef struct
{
  hash_table_t genus_reps;
  spinor_t spinor;
  fmpz_t disc;
  
  slong* dims;
  slong* conductors;
  slong** num_auts;
  slong** lut_positions;
  slong num_conductors;

  isometry_t* isoms; // isometries corresponding to the genus representatives
  
} genus_struct;

typedef genus_struct genus_t[1];

/* compute the genus of a quadratic form */
void genus_init_square_matrix(genus_t genus, const square_matrix_t Q, int h);

/* set the genus from a list */
void genus_init_set_square_matrix_vec(genus_t genus, const square_matrix_t* reps, size_t genus_size);
void genus_init_set_square_matrix_vec_and_isoms(genus_t genus, const square_matrix_t* reps,
						const isometry_t* isoms, size_t genus_size);

void genus_init_empty(genus_t genus, size_t disc);

void genus_init_file(genus_t genus, const char* genus_fname, size_t disc, bool has_isoms);

void genus_clear(genus_t genus);

bool genus_is_trivial_cond(const genus_t genus);

#endif // __GENUS_H__
