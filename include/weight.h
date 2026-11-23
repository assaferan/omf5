#ifndef __WEIGHT_H__
#define __WEIGHT_H__

#include "genus.h"
#include "square_matrix.h"

square_matrix_t** hecke_matrices_isometries(const genus_t genus, int p);

#include <flint/fmpq_mat.h>
#include <flint/fmpz_mpoly.h>

// the variables of the polynomials of a weight: x_1, ..., x_QF_RANK and y_1, ..., y_QF_RANK
#define WEIGHT_NVARS (2*QF_RANK)

/****************************************************************************
 *
 * Type: weight_rep_t
 *
 * Description: The representation of weight (a, b), a >= b >= 0, of the
 *              orthogonal group of the first lattice of a genus, with gram
 *              matrix q0: the polynomials f(x, y) in two vectors of QF_RANK
 *              variables, of degree a in x and b in y, which are harmonic
 *              for the laplacians of q0^(-1) in x, in y and in both, and are
 *              killed by sum_i x_i d/dy_i, on which g acts by
 *              f -> det(g)^(a + b) f(g^(-1) x, g^(-1) y). These form the
 *              irreducible representation of highest weight (a, b), of
 *              dimension (2a + 3)(2b + 1)(a - b + 1)(a + b + 2) / 6; with
 *              b = 0, the harmonic polynomials of degree a in x. The factor
 *              det(g)^(a + b) makes -1 act trivially, so the forms are those
 *              of the special orthogonal group, as in Magma.
 *
 *     + a, b (int) - the weight
 *     + nmons (slong) - the number of monomials of bidegree (a, b)
 *     + exps (ulong*) - their exponents, WEIGHT_NVARS for each monomial, x first
 *     + dim (slong) - the dimension of the representation
 *     + pivots (slong*) - the pivots of its basis, dim x nmons, in reduced row
 *                         echelon form, so that the coordinates of a vector
 *                         of it are its entries at the pivots
 *     + nx, ny (slong) - the numbers of monomials of degree a and b in
 *                        QF_RANK variables; the monomial k of bidegree
 *                        (a, b) is the product of the monomial k / ny of
 *                        degree a in x and k % ny of degree b in y
 *     + ex, ey (ulong*) - their exponents, QF_RANK for each monomial
 *     + basis_z (fmpz_mat_t) - the rows of the basis times their denominators,
 *                              nx x (dim ny): the entry (kx, k * ny + ky)
 *                              is the coefficient of the basis vector k at
 *                              the monomial kx * ny + ky
 *     + basis_den (fmpz*) - the denominator of each row of the basis
 *
 ***************************************************************************
 */

typedef struct
{
  int a;
  int b;
  slong nmons;
  ulong* exps;
  slong dim;
  slong* pivots;
  slong nx, ny;
  ulong *ex, *ey;
  fmpz_mat_t basis_z;
  fmpz* basis_den;
} weight_rep_struct;

typedef weight_rep_struct weight_rep_t[1];

void weight_rep_init(weight_rep_t W, int a, int b, const square_matrix_t q0);

// The number of entries of the largest integer matrices which weight_rep_init builds for the weight
// (a, b): the basis of the representation, dim x nmons, or a block of the conditions on the
// polynomials of bidegree (a, b). The memory is about 200 bytes for each.
double weight_rep_matrix_size(int a, int b);

// The largest size for which omf5 computes the forms of a weight (a, b) with b > 0: 2.5 10^7 entries,
// about 5 GB, which allows the weight (6, 6) and not (7, 7). The weights (a, 0) are not limited, as
// before the weights (a, b).
#define WEIGHT_MAX_MATRIX_SIZE 2.5e7

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
int print_hecke_matrices_weight(const genus_t genus, int a, int b, int p);

// Prints the eigenforms of the weight (a, b) of each conductor, with their eigenvalues at the primes ps
// (which have to be prime to the discriminant), as omf5 -prec does: {c : [{ 'aut_rep_type' : ...,
// 'field_poly' : ..., 'lambda_p' : [...], 'lambda_p_square' : [], }, ...], ...}. The types are not
// determined: U, or O for an old form. An eigenvalue of an old form is 'NULL' when the operator is not
// an element of the field of the form on its space (see ev_is_scalar_on_old_space). The matrices which
// the decomposition computed are used again, and the others are computed one prime at a time.
int print_eigenvalues_weight(const genus_t genus, int a, int b, const int* ps, int num_ps);

#endif // __WEIGHT_H__
