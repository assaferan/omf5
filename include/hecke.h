/*****************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : hecke.h
 *
 * Description: Functions for computing Hecke operators.
 *
 *****************************************************************
 */

#ifndef __HECKE_H__
#define __HECKE_H__

// Required packages dependencies

#include <carat/matrix.h>

#include <antic/nf_elem.h>

// Self dependencies

#include "decomposition.h"
#include "eigenvalues.h"
#include "genus.h"
#include "matrix_tools.h"
#include "nbr_data.h"
#include "neighbor.h"

// !! TODO - merge these when we finish debugging

/**********************************************************************
 *
 * Function: process_isotropic_vector_nbr_data
 *
 * Description: Update the values of the Hecke operator,
 *              given an isotropic vector, corresponding to a neighbor.
 *              In this version the neighbor manager is of type nbr_data.
 *
 * Arguments:
 *     + nbr_man (nbr_data_t) - a neighbor manager
 *     + T (int*) - a column of the matrix representing the Hecke operator
 *                  corresponding to the lattice whose neighbors we
 *                  are enumerating over
 *     + genus (const genus_t) - the associated genus
 *
 * Returns:
 *     + (int) - always returns 0
 *     + theta_time (double*) - time spent computing theta series
 *     + isom_time (double*) - time spent checking isometries
 *     + total_time (double*) - total time spent
 *     + num_isom (int*) - number of isometry tests done
 *
 **********************************************************************
 */

int process_isotropic_vector_nbr_data(nbr_data_t nbr_man, int* T, const genus_t genus, 
				      double* theta_time, double* isom_time, double* total_time, int* num_isom);

/**********************************************************************
 *
 * Function: process_isotropic_vector
 *
 * Description: Update the values of the Hecke operator,
 *              given an isotropic vector, corresopnding to a neighbor.
 *              In this version the neighbor manager is of type
 *              neighbor_manager.
 *
 * Arguments:
 *     + nbr_man (neighbor_manager_t) - a neighbor manager
 *     + T (int*) - a column of the matrix representing the Hecke operator
 *                  corresponding to the lattice whose neighbors we
 *                  are enumerating over
 *     + genus (const genus_t) - the associated genus
 *
 * Returns:
 *     + (int) - returns the size of the orbit of the isotropic vector
 *               under the automorphism group (number of neighbors covered)
 *     + theta_time (double*) - time spent computing theta series
 *     + isom_time (double*) - time spent checking isometries
 *     + total_time (double*) - total time spent
 *     + num_isom (int*) - number of isometry tests done
 *
 **********************************************************************
 */

int process_isotropic_vector(neighbor_manager_t nbr_man, int* T, const genus_t genus,
			     double* theta_time, double* isom_time, double* total_time, int* num_isom);

/**********************************************************************
 *
 * Function: process_neighbour_chunk
 *
 * Description: Update the values of the Hecke operator,
 *              on a chunk of neighbours, referred to by an index.
 *              In this version the neighbor manager is of type
 *              neighbor_manager.
 *
 * Arguments:
 *     + T (int*) - a column of the matrix representing the Hecke operator
 *                  corresponding to the lattice whose neighbors we
 *                  are enumerating over
 *     + p (int) - the prime (computing T_p)
 *     + i (int) - index of the chunk, specifying isotropic vectors
 *     + gen_idx (int) - the index in the genus of the lattice for
 *                       which we are computing the Hecke column
 *     + genus (const genus_t) - the associated genus
 *
 * Returns:
 *     + (int) - returns the number of neighbors processed
 *     + theta_time (double*) - time spent computing theta series
 *     + isom_time (double*) - time spent checking isometries
 *     + total_time (double*) - total time spent
 *     + num_isom (int*) - number of isometry tests done
 *
 **********************************************************************
 */

int process_neighbour_chunk(int* T, int p, int i, int gen_idx, const genus_t genus,
			    double* theta_time, double* isom_time, double* total_time, int* num_isom);

/**********************************************************************
 *
 * Function: process_neighbour_chunk_nbr_data
 *
 * Description: Update the values of the Hecke operator,
 *              on all neighbours.
 *              In this version the neighbor manager is of type
 *              nbr_data.
 *
 * Arguments:
 *     + T (int*) - a column of the matrix representing the Hecke operator
 *                  corresponding to the lattice whose neighbors we
 *                  are enumerating over
 *     + p (int) - the prime
 *     + k (int) - the index of the Hecke operator (computing T_{p,k})
 *     + gen_idx (int) - the index in the genus of the lattice for
 *                       which we are computing the Hecke column
 *     + genus (const genus_t) - the associated genus
 *
 * Returns:
 *     + (int) - returns the number of neighbors processed
 *     + theta_time (double*) - time spent computing theta series
 *     + isom_time (double*) - time spent checking isometries
 *     + total_time (double*) - total time spent
 *     + num_isom (int*) - number of isometry tests done
 *
 **********************************************************************
 */

int process_neighbour_chunk_nbr_data(int* T, int p, int k, int gen_idx, const genus_t genus, 
				     double* theta_time, double* isom_time, double* total_time, int* num_isom);

/**********************************************************************
 *
 * Function: hecke_col_nbr_data
 *
 * Description: Compute a column of a Hecke operator.
 *              In this version, the neighbor manager used is of type
 *              nbr_data.
 *
 * Arguments:
 *     + p (int) - the prime
 *     + k (int) - the index of hte Hecke operator (computing T_{p,k})
 *     + gen_idx (int) - the index in the genus of the lattice for
 *                       which we are computing the Hecke column
 *     + genus (const genus_t) - the associated genus
 *
 * Returns:
 *     + T (int*) - the column of the matrix representing the
 *                  Hecke operator T_{p,k}
 *                  corresponding to the lattice whose neighbors we
 *                  are enumerating over (genus[gen_idx])
 *
 **********************************************************************
 */

void hecke_col_nbr_data(int* T, int p, int k, int gen_idx, const genus_t genus);

/**********************************************************************
 *
 * Function: hecke_col_nbr_data_all_conductors
 *
 * Description: Compute a column of a Hecke operator for every conductor.
 *              In this version, the neighbor manager used is of type
 *              nbr_data.
 *
 * Arguments:
 *     + p (int) - the prime
 *     + k (int) - the index of hte Hecke operator (computing T_{p,k})
 *     + gen_idx (int) - the index in the genus of the lattice for
 *                       which we are computing the Hecke column
 *     + genus (const genus_t) - the associated genus
 *
 * Returns:
 *     + (slong) - number of neighbors processed
 *     + spin_vals (W64*) - a vector, where for every neighbor, we
 *                          store all the values of the spinor norm
 *                          as a 64-bit word - lit bits are for primes
 *                          where the spinor norm is -1. 
 *
 **********************************************************************
 */

slong hecke_col_nbr_data_all_conductors(W64* spin_vals, int p, int k, int gen_idx, const genus_t genus);

void hecke_col(int* T, int p, int gen_idx, const genus_t genus);

slong hecke_col_all_conductors(W64* spin_vals, int p, int gen_idx, const genus_t genus);

int** hecke_col_all_conds_sparse(int p, int col_idx, const genus_t genus);

matrix_TYP** hecke_matrices_all_conductors(const genus_t genus, int p);

void get_hecke_ev_nbr_data(nf_elem_t e, const genus_t genus, const eigenvalues_t evs, int p, int k, int ev_idx);

void get_hecke_ev_nbr_data_all_conductors(nf_elem_t e, const genus_t genus, const eigenvalues_t evs, int p, int k, int ev_idx, slong c);

void get_hecke_ev(nf_elem_t e, const genus_t genus, const eigenvalues_t evs, int p, int ev_idx);
void get_hecke_ev_all_conductors(nf_elem_t e, const genus_t genus,
				 const eigenvalues_t evs,
				 int p, int ev_idx, slong ev_cond);

matrix_TYP* hecke_matrix(const genus_t genus, int p);

void get_hecke_fmpq_mat_all_conductors(fmpq_mat_t* hecke_fmpq_mat, const genus_t genus, int p, int k);

void hecke_eigenforms(eigenvalues_t evs, const decomposition_t D, const genus_t genus, slong c);

eigenvalues_t* hecke_eigenforms_all_conductors(const genus_t genus);

#endif // __HECKE_H__
