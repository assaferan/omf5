#ifndef __WEIGHT_H__
#define __WEIGHT_H__

#include "genus.h"
#include "square_matrix.h"

square_matrix_t** hecke_matrices_isometries(const genus_t genus, int p);

#include <flint/fmpq_mat.h>
#include <flint/fmpq_mpoly.h>

/****************************************************************************
 *
 * Type: weight_rep_t
 *
 * Description: The representation of weight (a, 0) of the orthogonal group of
 *              the first lattice of a genus, with gram matrix q0: the
 *              polynomials of degree a in QF_RANK variables which are
 *              harmonic for the laplacian of q0^(-1), on which g acts by
 *              f -> det(g)^a f(g^(-1) x). The factor det(g)^a makes -1 act
 *              trivially, so the forms are those of the special orthogonal
 *              group, as in Magma.
 *
 *     + a (int) - the weight
 *     + nmons (slong) - the number of monomials of degree a
 *     + exps (ulong*) - their exponents, QF_RANK for each monomial
 *     + dim (slong) - the dimension of the representation
 *     + basis (fmpq_mat_t) - its basis, dim x nmons, in reduced row echelon
 *                            form, so that the coordinates of a vector of it
 *                            are its entries at the pivots
 *     + pivots (slong*) - the pivot of each row of basis
 *     + ctx (fmpq_mpoly_ctx_t) - the polynomial ring
 *     + polys (fmpq_mpoly_struct*) - the polynomials of the basis
 *
 ***************************************************************************
 */

typedef struct
{
  int a;
  slong nmons;
  ulong* exps;
  slong dim;
  fmpq_mat_t basis;
  slong* pivots;
  fmpq_mpoly_ctx_t ctx;
  fmpq_mpoly_struct* polys;
} weight_rep_struct;

typedef weight_rep_struct weight_rep_t[1];

void weight_rep_init(weight_rep_t W, int a, const square_matrix_t q0);

void weight_rep_clear(weight_rep_t W);

// R (W->dim x W->dim) is the matrix of the action of the isometry s of q0, in the basis of W:
// its column k is the coordinates of the image of the k-th basis vector.
void weight_rep_action(fmpq_mat_t R, const weight_rep_t W, const isometry_t s);

/****************************************************************************
 *
 * Type: weight_space_t
 *
 * Description: The spaces of algebraic modular forms of weight W on a genus,
 *              one for each conductor c: at the lattice L_i, the vectors w of
 *              W with chi_c(g) rho(g) w = w for every automorphism g of L_i,
 *              taken to the first lattice by genus->isoms[i].
 *
 *     + num_conductors (slong) - as in the genus
 *     + num_classes (slong) - the number of lattices of the genus
 *     + H (fmpq_mat_t**) - H[c][i], a basis of the space at L_i, as rows in the
 *                          coordinates of W, in reduced row echelon form
 *     + pivots (slong***) - pivots[c][i], the pivots of H[c][i]
 *     + offsets (slong**) - offsets[c][i], the index of the first basis vector
 *                           of L_i among those of the conductor c
 *     + dims (slong*) - dims[c], the dimension of the space of conductor c
 *
 ***************************************************************************
 */

typedef struct
{
  slong num_conductors;
  slong num_classes;
  fmpq_mat_struct** H;
  slong*** pivots;
  slong** offsets;
  slong* dims;
} weight_space_struct;

typedef weight_space_struct weight_space_t[1];

void weight_space_init(weight_space_t S, const genus_t genus, const weight_rep_t W);

void weight_space_clear(weight_space_t S);

// The matrices of T_p, p not dividing the discriminant, on the spaces of each conductor, with the
// normalization p^a, which makes the eigenvalues algebraic integers. Returns an array of
// genus->num_conductors matrices, to be cleared and freed by the caller.
// The neighbors are taken by orbits of the automorphisms of each lattice, so only one neighbor of each
// orbit is tested for isometry.
fmpq_mat_struct* hecke_matrices_weight(const genus_t genus, const weight_rep_t W, const weight_space_t S, int p);

// The same, with every neighbor tested for isometry (the neighbors of nbr_data). For tests.
fmpq_mat_struct* hecke_matrices_weight_by_neighbor(const genus_t genus, const weight_rep_t W, const weight_space_t S, int p);

// Prints the matrices of hecke_matrices_weight as omf5 -hecke does, {c : [[...], ...], ...}
int print_hecke_matrices_weight(const genus_t genus, int a, int p);

#endif // __WEIGHT_H__
