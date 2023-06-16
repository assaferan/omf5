/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : eigenvalues.h
 *
 * Description: data of eigenforms and their eigenvalues
 *              for a space of orthogonal modular forms.
 *
 **********************************************************
 */

#ifndef __EIGENVALUES_H__
#define __EIGENVALUES_H__

// Required packages dependencies

#include <antic/nf.h>
#include <antic/nf_elem.h>

#include <carat/typedef.h>

#include <flint/fmpz.h>

// Self dependencies

#include "genus.h"
#include "typedefs.h"

/****************************************************************
 *
 * Type: aut_type
 *
 * Description: automorphic type of an eigenform.
 *
 * An eigenform f gives rise to a automorphic representation pi,
 * whose type is one of F (Eisenstein), P (Saito-Kurokawa),
 * Y (Yoshida), or G (in the orthogonal complement of the above).
 * We added the type O to indicate oldforms.
 *
 ****************************************************************
 */

typedef enum {
  G,
  Y, // I thought we were not supposed to see Yoshida lifts, because paramodular ?
  // we see them in the spaces of orthogonal modular forms, but they do not correspond to paramodular forms
  P,
  F,
  O // oldform
} aut_type;

/****************************************************************
 *
 * Type: eigenvalues_t
 *
 * Description: this struct holds data for the eigenforms
 *              and their eigenvalues at a specific prime. 
 *
 * Fields:
 *     + nfs (nf_t*) - a list of the number fields K(f), over
 *                     which the Hecke eigenvalues of f are
 *                     defined
 *     + eigenvals (nf_elem_t*) - a list of Hecke eigenvalues,
 *                                storing the eigenvalue of a
 *                                single Hecke operator
 *     + eigenvecs (nf_elem_t**) - a list of vectors,
 *                                 representing the eigenforms
 *                                 in the basis of indicator
 *                                 functions for the genus
 *     + num (slong) - number of eigenforms
 *     + dim (slong) - dimension of the space
 *     + lift_type (aut_type*) - automorphic type of the form
 *
 * Used to hold the data of eigenforms and their eigenvalues,
 ****************************************************************
 */

typedef struct {
  nf_t* nfs;
  nf_elem_t* eigenvals; // right now only stores the eigenvalue of the Hecke operator used to construct it. rethink that
  nf_elem_t** eigenvecs; 
  slong num;
  slong dim;
  // bool* is_lift;
  aut_type* lift_type;
} eigenvalues;

typedef eigenvalues eigenvalues_t[1];

/****************************************************************
 *
 * Type: all_eigen
 *
 * Description: Data for all eigenvalues.
 *
 * Fields:
 *     + evs (eigenvalues_t*) - a list of eigenvalue structures,
 *                              one for each prime
 *     + num (slong) - number of eigenvalues
 *
 ****************************************************************
 */

typedef struct {
  eigenvalues_t* evs;
  slong num;
} all_eigen;

/******************************************************************
 *
 * Function: eigenvalues_init
 *
 * Description: Allocates memory for an eigenvalue struct.
 *
 * Arguments:
 *     + num (slong) - number of eigenvectors
 *     + dim (slong) - dimension of the space (length of vectors)
 *
 * Returns:
 *     + evs (eigenvalues_t) - the empty eigenvalues data
 *
 ******************************************************************
 */

void eigenvalues_init(eigenvalues_t evs, slong num, slong dim);

/*********************************************************************
 *
 * Function: eigenvalues_init_set_mat
 *
 * Description: Get eigenvectors and eigenvalues of a matrix.
 *
 * Arguments:
 *     + mat (matrix_TYP*) - the matrix (CARAT type)
 *
 * Returns:
 *     + evs (eigenvalues_t) - the eigenvectors and eigenvalues data
 *
 *********************************************************************
 */

void eigenvalues_init_set_mat(eigenvalues_t evs, matrix_TYP* mat);

/*********************************************************************
 *
 * Function: eigenvalues_set_lifts
 *
 * Description: Find which eigenforms are lifts,
 *              and set the automorphic types accordingly.
 *              Done by considering factorization of L-polynomials.
 *
 * Arguments:
 *     + evs (eigenvalues_t) - the eigenvectors
 *     + prec (slong) - number of primes at which to test the L-poly
 *     + c (slong) - the conductor
 *     + genus (const genus_t) - the genus 
 *
 * Returns:
 *     + evs (eigenvalues_t) - the eigenvectors with the aut_type set
 *
 *********************************************************************
 */

void eigenvalues_set_lifts(eigenvalues_t evs, slong prec, slong c, const genus_t genus);

/*********************************************************************
 *
 * Function: ev_get_pivot
 *
 * Description: Find a nonzero entry in the vector
 *              representing the eigenform, such that the
 *              corresponding lattice in the genus has a
 *              maximal number of automorphisms.
 *
 * Arguments:
 *     + evec (const nf_elem_t*) - the eigenvector representing f
 *     + K (nf_t) - the number field K(f)
 *     + genus (const genus_t) - the genus 
 *     + cond (slong) - the conductor
 *
 * Returns:
 *     + (slong) - index of the entry.
 *
 *********************************************************************
 */

slong ev_get_pivot(const nf_elem_t* evec, const nf_t K, const genus_t genus, slong cond);

/******************************************************************
 *
 * Function: eigenvalues_clear
 *
 * Description: Deallocates memory for an eigenvalues struct.
 *
 * Arguments:
 *     + evs (eigenvalues_t) - the eigenvalues data
 *
 ******************************************************************
 */

void eigenvalues_clear(eigenvalues_t evs);

#endif // __EIGENVALUES_H__
