/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : decomposition.h
 *
 * Description: data of Hecke module decomposition
 *
 **********************************************************
 */

#ifndef __DECOMPOSITION_H__
#define __DECOMPOSITION_H__

// Required packages dependencies

#include <flint/fmpq_mat.h>

// Self dependencies

#include "genus.h"

/****************************************************************
 *
 * Type: decomposition_t
 *
 * Description: this struct holds data for the decomposition
 *              of the spaces as Hecke modules.
 *
 * Fields:
 *     + bases (fmpq_mat_t**) - a list of matrices
 *                              for each conductor,
 *                              the rows of each matrix
 *                              form the basis for a submodule
 *     + hecke (fmpq_mat_t**) - a list of Hecke operators
 *                              for each prime
 *     + num (slong*) - number of submodules in the decomposition
 *                      for each conductor
                        (num[c] is the size of bases[c]) 
 *     + num_primes (slong) - number of primes for which we
 *                            computed Hecke operators.
 *                            (size of hecke)
 *     + num_conductors (slong) - number of conductors
 *                               (size of bases)
 *
 * Used to hold the data of a decomposition to Hecke submodules.
 ****************************************************************
 */

typedef struct {
  fmpq_mat_t** bases; // bases for the vector spaces, by conductor
  fmpq_mat_t** hecke;
  slong* num; // number of elements in the decomposition, by conductor
  slong num_primes;
  slong num_conductors;
} decomposition_struct;

typedef decomposition_struct decomposition_t[1];

/******************************************************************
 *
 * Function: decomposition_init
 *
 * Description: Allocates memory for a decomposition struct.
 *
 * Arguments:
 *     + num (slong) - number of different conductors (spaces)
 *
 * Returns:
 *     + decomp (decomposition_t) - the empty decomposition data
 *
 ******************************************************************
 */

void decomposition_init(decomposition_t decomp, slong num);

/******************************************************************
 *
 * Function: decomposition_clear
 *
 * Description: Deallocates memory for a decomposition struct.
 *
 * Arguments:
 *     + decomp (decomposition_t) - the decomposition data
 *
 ******************************************************************
 */

void decomposition_clear(decomposition_t decomp);

/******************************************************************
 *
 * Function: decompose
 *
 * Description: Decomposes the space of orthogonal modular forms
 *              of trivial weight for a specific genus,
 *              and conductor c, specifying the appropriate
 *              spinor norm twist.
 *
 * Arguments:
 *     + genus (const genus_t) - the genus on which the space
 *                               of modular forms is supported
 *     + c (slong) - the conductor
 *
 * Returns:
 *     + decomp (decomposition_t) - the full decomposition data
 *                                  including bases for all
 *                                  irreducible Hecke submodules
 *
 ******************************************************************
 */

void decompose(decomposition_t decomp, const genus_t genus, slong c);

#endif // __DECOMPOSITION_H__
