/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : tq_table.h
 *
 * Description: The table of T_q: the class of the q-neighbor of
 *              each class of the genus at each isotropic line mod q,
 *              with an isometry onto it (see issue #63).
 *
 **********************************************************
 */

#ifndef __TQ_TABLE_H__
#define __TQ_TABLE_H__

#include "genus.h"
#include "isometry.h"
#include "square_matrix.h"
#include "typedefs.h"

/***********************************************************************
 *
 * Type: tq_table_t
 *
 * Description: for a prime q which does not divide the discriminant, and
 *              each class j of the genus and isotropic line w of L_j mod q,
 *              the class j' of the q-neighbor N of L_j at w and an isometry
 *              from L_j onto L_{j'} which takes L_j to N first:
 *              isom = s_nbr * hash_isom, where s_nbr : L_j -> N is the
 *              isometry of nbr_data_build_neighbor and hash_isom : N -> L_{j'}
 *              that of hash_table_index_and_isom. In the convention of
 *              isometry_is_isom, it is an isometry from keys[j] to keys[j'].
 *              The lines are indexed by tq_table_line_index.
 *
 * Fields:
 *     + q (int) - the prime
 *     + num_classes (slong) - the number of classes of the genus
 *     + size (slong) - q^5, the number of indices of lines
 *     + cls (int**) - cls[j][idx] is j', or -1 if idx is not an
 *                     isotropic line of L_j outside its radical
 *     + isom (isometry_t**) - isom[j][idx], when cls[j][idx] >= 0
 *     + spin (W64**) - spin[j][idx], the entry of spin_vals which
 *                      hecke_col_nbr_data_all_conductors gives for this
 *                      neighbor: (j' << num_primes) | spinor bits
 *     + num_lines (slong*) - the number of lines of each class
 *
 ***********************************************************************
 */

typedef struct
{
  int q;
  slong num_classes;
  slong size;
  int** cls;
  isometry_t** isom;
  W64** spin;
  slong* num_lines;
} tq_table_struct;

typedef tq_table_struct tq_table_t[1];

// The index of the line of the vector w mod q (w not 0 mod q): the normalized vector, with its first
// nonzero entry 1, read as an integer in base q.
slong tq_table_line_index(const vector_t w, int q);

// Builds the table of T_q for the genus, with one isometry test for each class and line. q must be a
// prime which does not divide the discriminant.
void tq_table_init(tq_table_t T, const genus_t genus, int q);

void tq_table_clear(tq_table_t T);

#endif // __TQ_TABLE_H__
