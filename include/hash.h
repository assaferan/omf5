/*****************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : hash.h
 *
 * Description: Hash table for quickly testing
 *              isometry between lattices.
 *
 *****************************************************************
 */

#ifndef __HASH_H__
#define __HASH_H__

// System dependencies

#include <stdint.h>

// Required packacges dependencies

#include <carat/matrix.h>
#include <flint/fmpq.h>

// Self dependencies

#include "isometry.h"
#include "square_matrix.h"
#include "typedefs.h"

// hash_t is the type used for the hash function.
// We do not anticipate more than 2^32 forms in the genus...
// Maybe changing to W16 could improve that, have to check

typedef W32 hash_t;

/***********************************************************************
 *
 * Type: hash_table_t
 *
 * Description: A (closed) hash table for isometry classes of lattices.
 *              The keys are gram matrices of lattices, representing
 *              the isometry classes.
 *
 * Fields:
 *     + keys (square_matrix_t*) - the list of keys stored in the table  
 *     + vals (hash_t*) - the list of corresponding hash values (full word)
 *                        Note that the actual hashing index is obtained by
 *                        considering only some of its bits
 *     + key_ptr (hash_t*) - for any hash value v, key_ptr[v] is the index
 *                           in the vals array in which v appears
 *     + counts (W32*) - the number of times every hash value appears,
 *                       for resolving collisions
 *     + probs (fmpq_t*) - probabilities of hitting values in the table
 *     + num_stored (hash_t) - number of keys stored in the table
 *     + mask (hash_t) - mask used to determine which bits of the full
 *                       hash value to use
 *     + capacity (hash_t) - capacity of the table
 *                           (how many entries can the allocated memory
 *                            currently support)
 *     + theta_prec (W32) - Our hash function uses the theta series
 *                          associated to the lattice. This determines
 *                          how many terms from the theta series
 *                          to consider.
 *     + red_on_isom (bool) - a flag which determines, when checking
 *                            for an isometry, whether to first reduce
 *                            the lattice using a naive reduction
 *                            algorithm.
 *
 ***********************************************************************
 */

typedef struct {
  square_matrix_t* keys;
  hash_t* vals;
  hash_t* key_ptr;
  W32* counts;
  fmpq_t* probs;

  hash_t num_stored;
  hash_t mask; 
  hash_t capacity;
  W32 theta_prec;

  bool red_on_isom;
} hash_table;

typedef hash_table hash_table_t[1];

/**********************************************************************
 *
 * Function: hash_form
 *
 * Description: Hash function - hash the lattice into a hash_value.
 *
 * Arguments:
 *     + Q (const square_matrix_t) - the gram matrix of the lattice
 *     + theta_prec (W32) - the number of terms in the theta series
 *                          f the lattice to use.
 *
 * Returns:
 *     + (hash_t) - the hash value
 *
 **********************************************************************
 */

/* hash the form Q into an index between 0 and hash_size */
hash_t hash_form(const square_matrix_t Q, W32 theta_prec);

/**********************************************************************
 *
 * Function: hash_table_init
 *
 * Description: Initializes a hash table of given size.
 *
 * Arguments:
 *     + hash_size (hash_t) - size of the hash table
 *
 * Returns:
 *     + table (hash_table_t) - a hash table of that size
 *
 **********************************************************************
 */

void hash_table_init(hash_table_t table, hash_t hash_size);

/**********************************************************************
 *
 * Function: hash_table_recalibrate
 *
 * Description: Change the number of terms we are taking from the
 *              theta series of the lattices, to optimize performance.
 *
 * Arguments:
 *     + old_table (const hash_table_t) - current hash_table
 *
 * Returns:
 *     + new_table (hash_table_t) - a hash_table with the same
 *                                  keys, where the hash function
 *                                  uses the new number of terms
 *                                  from the theta series
 *
 **********************************************************************
 */

void hash_table_recalibrate(hash_table_t new_table, const hash_table_t old_table);

/**********************************************************************
 *
 * Function: hash_table_clear
 *
 * Description: Clears memory allocated for hash table.
 *
 * Arguments:
 *     + table (hash_table_t) - the hash table
 *
 **********************************************************************
 */

void hash_table_clear(hash_table_t table);

/**********************************************************************
 *
 * Function: hash_table_add
 *
 * Description: Adds an entry to the hash table.
 *
 * Arguments:
 *     + table (hash_table_t) - the hash table
 *     + key (const square_matrix_t) - the gram matrix of a lattice
 *                                     in the isometry class.
 *
 * Returns:
 *     + (int) - 0 if there is already a lattice from this isometry
 *               class in the table, 1 otherwise.
 *
 **********************************************************************
 */

int hash_table_add(hash_table_t table, const square_matrix_t key);

/**********************************************************************
 *
 * Function: hash_table_get_key
 *
 * Description: Given a gram matrix of a lattice, obtain a gram
 *              matrix of a lattice stored in the table, with the same
 *              hash value. (not necessarily isometric)
 *
 * Arguments:
 *     + table (const hash_table_t) - the hash table
 *     + key (const square_matrix_t) - the gram matrix of the lattice
 *     + index (int*) - a pointer for resolving collisions,
 *                      which points to the first index in the table
 *                      from which to search (for the
 *                      first occurrence, set *index = -1)
 *
 * Returns:
 *     + new_key (square_matrix_t) - a gram matrix of a lattice with
 *                                   the same hash value.
 *     + (bool) - True if such a lattice was found.
 *
 **********************************************************************
 */

bool hash_table_get_key(square_matrix_t new_key, const hash_table_t table, const square_matrix_t key, int* index);

/**********************************************************************
 *
 * Function: hash_table_exists
 *
 * Description: Given a gram matrix of a lattice, check whether
 *              there is a lattice stored in the table with the same
 *              hash value. If check_isom is set, check whether
 *              there exists such an isometric lattice.
 *
 * Arguments:
 *     + table (const hash_table_t) - the hash table
 *     + key (const square_matrix_t) - the gram matrix of the lattice
 *     + check_isom (int) - 1 for checking an isometry as well
 *
 * Returns:
 *     + (int) - 1 if such a lattice was found, 0 otherwise.
 *
 **********************************************************************
 */

int hash_table_exists(const hash_table_t table, const square_matrix_t key, int check_isom);

/**********************************************************************
 *
 * Function: hash_table_indexof
 *
 * Description: Given a gram matrix of a lattice, find the index of
 *              a lattice stored in the table with the same
 *              hash value. If check_isom is set, find the index of
 *              an isometric lattice stored in the table.
 *
 * Arguments:
 *     + table (const hash_table_t) - the hash table
 *     + key (const square_matrix_t) - the gram matrix of the lattice
 *     + check_isom (int) - 1 for checking an isometry as well
 *
 * Returns:
 *     + (int) - the index in the table in which such a lattice occurs.
 *               If not found, returns -1.
 *     + theta_time (double*) - for recalibration, the amount of time
 *                              spent computing the theta series
 *     + isom_time (double*) - for recalibration, the amount of time
 *                             spent checking isometries
 *     + num_isom (int *) - for recalibration, the number of isometry
 *                          tests performed
 *
 **********************************************************************
 */

int hash_table_indexof(const hash_table_t table, const square_matrix_t key, int check_isom, double* theta_time, double* isom_time, int* num_isom);

/**********************************************************************
 *
 * Function: hash_table_index_and_isom
 *
 * Description: Given a gram matrix of a lattice, find the index of
 *              a lattice stored in the table with the same
 *              hash value. If check_isom is set, find the index of
 *              an isometric lattice stored in the table.
 *              In addition, returns an isometry between the
 *              input lattice and the stored lattice.
 *
 * Arguments:
 *     + table (const hash_table_t) - the hash table
 *     + key (const square_matrix_t) - the gram matrix of the lattice
 *     + check_isom (int) - 1 for checking an isometry as well
 *
 * Returns:
 *     + (int) - the index in the table in which such a lattice occurs.
 *               If not found, returns -1.
 *     + isom (isometry_t) - an isometry between key and the stored
 *                           lattice
 *     + theta_time (double*) - for recalibration, the amount of time
 *                              spent computing the theta series
 *     + isom_time (double*) - for recalibration, the amount of time
 *                             spent checking isometries
 *     + num_isom (int *) - for recalibration, the number of isometry
 *                          tests performed
 *
 **********************************************************************
 */

int hash_table_index_and_isom(const hash_table_t table, const square_matrix_t key, isometry_t isom,
			      double* theta_time, double* isom_time, int* num_isom);

/**********************************************************************
 *
 * Function: hash_table_expand
 *
 * Description: Expand the hash table (increase its capacity).
 *
 * Arguments:
 *     + table (hash_table_t) - the hash table
 *
 **********************************************************************
 */

void hash_table_expand(hash_table_t table);

/**********************************************************************
 *
 * Function: hash_table_insert
 *
 * Description: Inserts an isometry class to the hash table, at
 *              a specific index.
 *
 * Arguments:
 *     + table (hash_table_t) - the hash table
 *     + key (const square_matrix_t) - the gram matrix of the lattice
 *     + val (hash_t) - the value of the hash function
 *     + index (int) - the index in which to insert the class
 *     + do_push_back (int) - 1 to store the lattice, 0 otherwise.
 *
 **********************************************************************
 */

// !! TODO - this seems to be an internal function, why do we need it in the include file?
int hash_table_insert(hash_table_t table, const square_matrix_t key, hash_t val,
		      int index, int do_push_back);

#endif // __HASH_H__
