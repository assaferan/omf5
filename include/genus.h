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

/**********************************************************************
 *
 * Function: genus_init_square_matrix
 *
 * Description: Initializes a genus from the gram matrix of a lattice.
 *              It uses Kneser's neighbors method to produce other
 *              lattices in the genus, and the mass formula for
 *              completeness.
 *
 * Arguments:
 *     + Q (const square_matrix_t) - the gram matrix of the lattice
 *     + h (int) - the number of genus representatives
 *                 If h is unknown in advance, set h = -1
 *
 * Returns:
 *     + genus (genus_t) - the resulting genus structure
 *
 **********************************************************************
 */


/* compute the genus of a quadratic form */
void genus_init_square_matrix(genus_t genus, const square_matrix_t Q, int h);

/**************************************************************************
 *
 * Function: genus_init_square_matrix_vec
 *
 * Description: Initializes a genus from a list of
 *              gram matrices of genus representatives.
 *              
 * Arguments:
 *     + reps (const square_matrix_t*) - the gram matrices of the lattices
 *     + genus_size (size_t) - the number of genus representatives
 *
 * Returns:
 *     + genus (genus_t) - the resulting genus structure
 *
 **************************************************************************
 */

/* set the genus from a list */
void genus_init_set_square_matrix_vec(genus_t genus, const square_matrix_t* reps, size_t genus_size);


/**************************************************************************
 *
 * Function: genus_init_square_matrix_vec_and_isoms
 *
 * Description: Initializes a genus from a list of
 *              gram matrices of genus representatives,
 *              together with specified isometries (over Q) between them.
 *              Useful for consistency checks.
 *              
 * Arguments:
 *     + reps (const square_matrix_t*) - the gram matrices of the lattices
 *     + isoms (const isometry_t*) - the isometries between the lattices
 *     + genus_size (size_t) - the number of genus representatives
 *
 * Returns:
 *     + genus (genus_t) - the resulting genus structure
 *
 **************************************************************************
 */

void genus_init_set_square_matrix_vec_and_isoms(genus_t genus, const square_matrix_t* reps,
						const isometry_t* isoms, size_t genus_size);

/**************************************************************************
 *
 * Function: genus_init_empty
 *
 * Description: Initializes an empty genus from the discriminant.
 *              Suffices for figuring out the possible conductors
 *              and spinor norm data, without computing the genus
 *              representatives.
 *              
 * Arguments:
 *     + disc (size_t) - the discriminant of the lattice
 *
 * Returns:
 *     + genus (genus_t) - the resulting genus structure
 *
 **************************************************************************
 */

void genus_init_empty(genus_t genus, size_t disc);

/********************************************************************************
 *
 * Function: genus_init_file
 *
 * Description: Initializes a genus from a file.
 *              
 * Arguments:
 *     + genus_fname (const char*) - the name of the file containing the genus.
 *                                   If has_isoms is true, the file is assumed
 *                                   to be a list of pairs of lists of matrices,
 *                                   the first corresponding to gram matrices
 *                                   of the genus representatives, and the second
 *                                   to the isometries, indexed by the level.
 *                                   If has_isoms is false, the file is assumed
 *                                   to consist of a list of lists of numbers,
 *                                   each list corresponding to a symmetric matrix
 *                                   by listing its lower triangular part,
 *                                   so that each list describes the genus
 *                                   representatives.
 *                                   The list is again indexed by level.
 *     + disc (size_t) - the discriminant of the lattice
 *     + has_isoms (bool) - whether the genus file contains the isometries
 *
 * Returns:
 *     + genus (genus_t) - the resulting genus structure
 *
 ******************************************************************************
 */

void genus_init_file(genus_t genus, const char* genus_fname, size_t disc, bool has_isoms);

/**************************************************************************
 *
 * Function: genus_clear
 *
 * Description: Clears memory allocated to initialize genus.
 *              
 * Arguments:
 *     + genus (genus_t) - the genus structure
 *
 **************************************************************************
 */

void genus_clear(genus_t genus);

/**************************************************************************
 *
 * Function: genus_is_trivial_cond
 *
 * Description: Checks whether the spaces of orthogonal modular forms
 *              for this genus with non-trivial spinor norm are all
 *              zero-dimensional. In such a case, we do not have to
 *              record the isometries.
 *              
 * Arguments:
 *     + genus (genus_t) - the genus structure
 *
 * Returns:
 *     + (bool) - True if spaces for all non-trivial conductors
 *                are zero-dimensional.
 *
 **************************************************************************
 */

bool genus_is_trivial_cond(const genus_t genus);

#endif // __GENUS_H__
