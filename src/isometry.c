/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : isometry.c
 *
 * Description: Data structure for isometries.
 *
 **********************************************************
 */

// System dependencies

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

// Required packages dependencies

#include <carat/autgrp.h>
#include <carat/symm.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

// Self dependencies

#include "isometry.h"
#include "matrix_tools.h"

/* Initialize an isometry - sets it to the identity. */
void isometry_init(isometry_t isom)
{
  // by default we initialize isometries to one
  square_matrix_one(isom->s);
  square_matrix_one(isom->s_inv);
  isom->denom = 1;
  isom->inv_denom = 1;

  return;
}

/* clone an isometry */
void isometry_init_set(isometry_t dest, const isometry_t src)
{
  square_matrix_set(dest->s, src->s);
  square_matrix_set(dest->s_inv, src->s_inv);
  dest->denom = src->denom;
  dest->inv_denom = src->inv_denom;

  return;
}

/* initialize an isometry from a rational matrix (using square_matrix_t) */
void isometry_init_set_square_matrix(isometry_t isom, const square_matrix_t s, Z64 denom)
{
  // maybe we want to save time copying - check later
  square_matrix_set(isom->s,s);
  isom->denom = denom;
  isom->inv_denom = square_matrix_inv(isom->s_inv, s, isom->denom);

  return;
}

/* initialize an isometry from a rational matrix (using fmpz_mat_t) */
void isometry_init_set_fmpz_mat(isometry_t isom, const fmpz_mat_t s, Z64 denom)
{

  square_matrix_set_fmpz_mat(isom->s, s);
  isom->denom = denom;
  isom->inv_denom = square_matrix_inv(isom->s_inv, isom->s, isom->denom);

  return;
}

/* Sets prod / prod_denom to the product of matL / denomL and matR / denomR.
   The numerators and the denominators are 64 bit integers. When the products
   of the numerators and of the denominators fit in them, these are the
   result, and the fraction is not reduced. Otherwise the product is found
   with big integers and the fraction is reduced, and it is an error when it
   does not fit after that. The products used to be found with no check, and
   the denominators were of type int, so that composing isometries with large
   denominators gave matrices which were not isometries.
   prod can be one of matL and matR, and prod_denom can point to one of the
   denominators. */
static void rational_matrix_mul(square_matrix_t prod, Z64* prod_denom,
				const square_matrix_t matL, Z64 denomL,
				const square_matrix_t matR, Z64 denomR)
{
  square_matrix_t res;
  Z64 res_denom, term;
  fmpz_mat_t L, R, P;
  fmpz_t d, g;
  int i, j, k;
  bool fits;

  fits = !__builtin_mul_overflow(denomL, denomR, &res_denom);
  for (i = 0; fits && (i < QF_RANK); i++)
    for (j = 0; fits && (j < QF_RANK); j++) {
      res[i][j] = 0;
      for (k = 0; fits && (k < QF_RANK); k++)
	fits = !__builtin_mul_overflow(matL[i][k], matR[k][j], &term) &&
	  !__builtin_add_overflow(res[i][j], term, &(res[i][j]));
    }

  if (!fits) {
    fmpz_mat_init_set_square_matrix(L, matL);
    fmpz_mat_init_set_square_matrix(R, matR);
    fmpz_mat_init(P, QF_RANK, QF_RANK);
    fmpz_init(d);
    fmpz_init(g);

    fmpz_mat_mul(P, L, R);
    fmpz_set_si(d, denomL);
    fmpz_mul_si(d, d, denomR);

    // the denominators are not zero, so neither is g
    fmpz_mat_content(g, P);
    fmpz_gcd(g, g, d);
    fmpz_mat_scalar_divexact_fmpz(P, P, g);
    fmpz_divexact(d, d, g);

    fits = fmpz_fits_si(d);
    for (i = 0; fits && (i < QF_RANK); i++)
      for (j = 0; fits && (j < QF_RANK); j++)
	fits = fmpz_fits_si(fmpz_mat_entry(P, i, j));
    if (!fits) {
      fprintf(stderr, "Error: a product of isometries does not fit in 64 bits\n");
      exit(EXIT_FAILURE);
    }
    square_matrix_set_fmpz_mat(res, P);
    res_denom = fmpz_get_si(d);

    fmpz_mat_clear(L);
    fmpz_mat_clear(R);
    fmpz_mat_clear(P);
    fmpz_clear(d);
    fmpz_clear(g);
  }

  square_matrix_set(prod, res);
  *prod_denom = res_denom;

  return;
}

/* multiply (compose) two isometries.
   we do not reduce the denominator when doing so, unless the product does
   not fit otherwise.
   prod can be one of s1 and s2. */
void isometry_mul(isometry_t prod, const isometry_t s1, const isometry_t s2)
{
  rational_matrix_mul(prod->s, &(prod->denom), s1->s, s1->denom, s2->s, s2->denom);
  rational_matrix_mul(prod->s_inv, &(prod->inv_denom), s2->s_inv, s2->inv_denom, s1->s_inv, s1->inv_denom);

  return;
}

/* multiply in place, from the right */
void isometry_muleq_right(isometry_t sL, const isometry_t sR)
{
  isometry_mul(sL, sL, sR);

  return;
}

/* multiply in place, from the left */
void isometry_muleq_left(isometry_t sR, const isometry_t sL)
{
  isometry_mul(sR, sL, sR);

  return;
}

/* multiply from the left by a rational matrix */
void isometry_mul_mat_left(isometry_t prod, const square_matrix_t s1, Z64 s1_denom, const isometry_t s2)
{
  rational_matrix_mul(prod->s, &(prod->denom), s1, s1_denom, s2->s, s2->denom);
  prod->inv_denom = square_matrix_inv(prod->s_inv, prod->s, prod->denom);

  return;
}

/* multiply from the right by a rational matrix */
void isometry_mul_mat_right(isometry_t prod, const isometry_t s1, const square_matrix_t s2, Z64 s2_denom)
{
  rational_matrix_mul(prod->s, &(prod->denom), s1->s, s1->denom, s2, s2_denom);
  prod->inv_denom = square_matrix_inv(prod->s_inv, prod->s, prod->denom);

  return;
}

/* multiply from the left, in place, by a rational matrix */
void isometry_mul_mateq_left(isometry_t sR, const square_matrix_t sL, Z64 sL_denom)
{
  isometry_mul_mat_left(sR, sL, sL_denom, sR);

  return;
}

/* multiply from the right, in place, by a rational matrix */
void isometry_mul_mateq_right(isometry_t sL, const square_matrix_t sR, Z64 sR_denom)
{
  isometry_mul_mat_right(sL, sL, sR, sR_denom);

  return;
}

/* invert */
void isometry_inv(isometry_t inv, const isometry_t isom)
{
  
  square_matrix_set(inv->s, isom->s_inv);
  square_matrix_set(inv->s_inv, isom->s);
  inv->denom = isom->inv_denom;
  inv->inv_denom = isom->denom;

  return;
}

/* clear memory allocated */
void isometry_clear(isometry_t isom)
{
  // did not allocate any memory
  return;
}

/* compute the new gram matrix after applying the isometry (g^t * Q * g) */
// !! TODO - we can make this operation faster by doing it all at once
void isometry_transform_gram(square_matrix_t gtQg, const isometry_t g, const square_matrix_t Q)
{
  int i,j,k,l;
  square_matrix_t /* gt, Qg, */ gtQg2;

  /*
  square_matrix_transpose(gt, g->s);
  square_matrix_mul(Qg, Q, g->s);
  */
  
  for (i = 0; i < QF_RANK; i++)
    for (j = i; j < QF_RANK; j++) {
      gtQg2[i][j] = 0;
      for (k = 0; k < QF_RANK; k++)
	for (l = 0; l < QF_RANK; l++)
	  gtQg2[i][j] += Q[k][l] * (g->s)[k][i] * (g->s)[l][j];
    }

  // square_matrix_mul(gtQg, gt, Qg);
  
  for (i = 0; i < QF_RANK; i++) {
    for (j = 0; j < i; j++)
      gtQg[i][j] = gtQg2[j][i];
    for (j = i; j < QF_RANK; j++)
      gtQg[i][j] = gtQg2[i][j];
  }

  //  assert(square_matrix_is_equal(gtQg, gtQg2));
  
  // the square of the denominator might overflow
  // square_matrix_div_scalar(gtQg, gtQg, (g->denom)*(g->denom));
  square_matrix_div_scalar(gtQg, gtQg, g->denom);
  square_matrix_div_scalar(gtQg, gtQg, g->denom);

  return;
}

/* Checks whether s is an isometry from q1 to q2, i.e. if s^t * q1 * s = q2. */
// !! TODO - s should be const here
bool isometry_is_isom(isometry_t s, const square_matrix_t q1, const square_matrix_t q2)
{
  fmpz_mat_t S, St, Q1, Q2;
  bool is_isom;

#ifdef DEBUG_LEVEL_FULL
  printf("s = \n");
  square_matrix_print(s->s);
  printf("\n");
  printf("q1 = \n");
  square_matrix_print(q1);
  printf("\n");
  printf("q2 = \n");
  square_matrix_print(q2);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  // The matrix of s is S / d for S = s->s and d = s->denom, so s is an isometry when
  // S^t * q1 * S = d^2 * q2. This is not checked by isometry_transform_gram, since it rounds
  // when it divides by d, and its entries might overflow. Then what is not an isometry could
  // be found to be one.
  if (s->denom == 0)
    return false;

  fmpz_mat_init_set_square_matrix(S, s->s);
  fmpz_mat_init_set_square_matrix(Q1, q1);
  fmpz_mat_init_set_square_matrix(Q2, q2);
  fmpz_mat_init(St, QF_RANK, QF_RANK);

  fmpz_mat_transpose(St, S);
  fmpz_mat_mul(Q1, St, Q1);
  fmpz_mat_mul(Q1, Q1, S);
  fmpz_mat_scalar_mul_si(Q2, Q2, s->denom);
  fmpz_mat_scalar_mul_si(Q2, Q2, s->denom);

  is_isom = fmpz_mat_equal(Q1, Q2);

  fmpz_mat_clear(S);
  fmpz_mat_clear(St);
  fmpz_mat_clear(Q1);
  fmpz_mat_clear(Q2);

  return is_isom;
}

/* Checks whether q1 and q2 are isometric, and if so computes an isometry betwenn them.
   Uses isometry testing from the CARAT package. */
/* The search of CARAT in the bases in which the lattices are given. It lists the vectors of both
   lattices up to the largest diagonal entry of q1, so it is slow and needs much memory when that
   entry is large, as it is in a basis which is far from reduced. */
static bool is_isometric_in_bases(isometry_t s, const square_matrix_t q1, const square_matrix_t q2)
{
  int i, Qmax;
  matrix_TYP *SV1, *SV2, *isom, *Q1, *Q2;
  square_matrix_t s_sm, s_tr;
  int options[6] = {0};

  assert(square_matrix_is_positive_definite(q1));
  assert(square_matrix_is_positive_definite(q2));
  
  Q1 = matrix_TYP_init_set_square_matrix(q1);
  Q2 = matrix_TYP_init_set_square_matrix(q2);
  
  Qmax = Q1->array.SZ[0][0];
  for(i=1;i < Q1->cols;i++) {
    if(Q1->array.SZ[i][i] > Qmax)
      Qmax = Q1->array.SZ[i][i];
  }
  SV1 = short_vectors(Q1, Qmax, 0, 0, 0, &i);
  SV2 = short_vectors(Q2, Qmax, 0, 0, 0, &i);

  isom = isometry(&Q1, &Q2, 1, SV1, SV2, NULL, 0, options);
  if (isom != NULL) {
    square_matrix_set_matrix_TYP(s_sm, isom);
    square_matrix_transpose(s_tr, s_sm);
    isometry_init_set_square_matrix(s, s_tr, 1);
    free_mat(isom);
  }

  free_mat(SV1);
  free_mat(SV2);
  free_mat(Q1);
  free_mat(Q2);

  return (isom != NULL);
}

/* CARAT lists the vectors of both lattices up to the largest diagonal entry of q1, which is large
   when the basis of q1 is far from reduced, and it stops when the coordinates of a vector of q2 are
   above 10000, as they can be when the basis of q2 is far from reduced. So a lattice whose basis is
   far from reduced is replaced by a reduction of it for the search, and the isometry is moved back
   to the bases given. Whether the lattices are isometric does not depend on the bases. When the
   isometry or its inverse does not fit in 64 bits in the bases given, the search is done in them. */
bool is_isometric(isometry_t s, const square_matrix_t q1, const square_matrix_t q2)
{
  square_matrix_t r1, r2, s_given, s_inv;
  fmpz_mat_t t1, t2;
  isometry_t s_red;
  Z64 inv_denom;
  bool is_reduced, is_isom;

  square_matrix_set(r1, q1);
  square_matrix_set(r2, q2);
  fmpz_mat_init(t1, QF_RANK, QF_RANK);
  fmpz_mat_init(t2, QF_RANK, QF_RANK);
  fmpz_mat_one(t1);
  fmpz_mat_one(t2);

  is_reduced = false;
  if (square_matrix_is_far_from_reduced(q1) && square_matrix_reduce_basis(r1, t1))
    is_reduced = true;
  if (square_matrix_is_far_from_reduced(q2) && square_matrix_reduce_basis(r2, t2))
    is_reduced = true;

  if (!is_reduced)
    is_isom = is_isometric_in_bases(s, q1, q2);
  else {
    isometry_init(s_red);
    is_isom = is_isometric_in_bases(s_red, r1, r2);
    if (is_isom) {
      if (square_matrix_change_basis_fits(s_given, t1, s_red->s, t2)
	  && square_matrix_inv_fits(s_inv, &inv_denom, s_given, 1)) {
	square_matrix_set(s->s, s_given);
	square_matrix_set(s->s_inv, s_inv);
	s->denom = 1;
	s->inv_denom = inv_denom;
      }
      else
	is_isom = is_isometric_in_bases(s, q1, q2);
    }
    isometry_clear(s_red);
  }

  assert(!is_isom || isometry_is_isom(s, q1, q2));

  fmpz_mat_clear(t1);
  fmpz_mat_clear(t2);

  return is_isom;
}

// Matrix operations

/* swap two basis vectors */
void isometry_swap_vecs(isometry_t isom, int col1, int col2)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  Z64 denom;
#endif // DEBUG
  int row;
  
  for (row = 0; row < QF_RANK; row++)
    square_matrix_swap_elts(isom->s, row, col1, row, col2);
  for (row = 0; row < QF_RANK; row++)
    square_matrix_swap_elts(isom->s_inv, col1, row, col2, row);

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG
  
  return;
}

/* add a mupltiple of one basis vector to another */
void isometry_add_vec(isometry_t isom, int target_col, Z64 c, int src_col)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  Z64 denom;
#endif // DEBUG
  int row;

  for (row = 0; row < QF_RANK; row++)
    isom->s[row][target_col] += c * isom->s[row][src_col];

  for (row = 0; row < QF_RANK; row++)
    isom->s_inv[src_col][row] -= c * isom->s_inv[target_col][row];

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG
  
  return;
}

/* replace a basis vector by a given vector */
void isometry_replace_vec(isometry_t isom, int col, const vector_t x)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  Z64 denom;
#endif // DEBUG
  int row;
  
  for (row = 0; row < QF_RANK; row++)
    isom->s[row][col] = isom->s[row][0];
  for (row = 0; row < QF_RANK; row++)
    isom->s[row][0] = x[row];
  isom->inv_denom = isom->s[col][0];

  // and updating s_inv
  isom->s_inv[0][0] = 0;
  isom->s_inv[col][0] = 1;
  isom->s_inv[0][col] = 1;
  for (row = 1; row < QF_RANK; row++)
    isom->s_inv[row][col] = -isom->s[row][0];
  if (col != 0)
    isom->s_inv[col][col] = -isom->s[0][0];

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG
  
  return;
}

/* divide a basis vector by a scalar */
void isometry_vec_scalar_div(isometry_t isom, int col, Z64 scalar)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  Z64 denom;
#endif // DEBUG
  
  int i,j;

  isom->denom *= scalar;
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      if (j != col)
	isom->s[i][j] *= scalar;

  for (i = 0; i < QF_RANK; i++)
    isom->s_inv[col][i] *= scalar;

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG

  return;
}

/* multiply a basis vector by a scalar */
void isometry_vec_scalar_mul(isometry_t isom, int col, Z64 scalar)
{
#ifdef DEBUG
  square_matrix_t s_inv;
  Z64 denom;
#endif // DEBUG
  
  int i,j;

  isom->inv_denom *= scalar;
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      if (j != col)
	isom->s_inv[j][i] *= scalar;

  for (i = 0; i < QF_RANK; i++)
    isom->s[i][col] *= scalar;

#ifdef DEBUG
  denom = square_matrix_inv(s_inv, isom->s, isom->denom);
  assert(square_matrix_is_equal(s_inv, isom->s_inv));
  assert(denom == isom->inv_denom);
#endif // DEBUG

  return;
}

/* Check if two isometries are equal */
bool isometry_is_equal(const isometry_t isom1, const isometry_t isom2)
{
  square_matrix_t scaled_isom1, scaled_isom2;

  square_matrix_mul_scalar(scaled_isom1, isom1->s, isom2->denom);
  square_matrix_mul_scalar(scaled_isom2, isom2->s, isom1->denom);

  return square_matrix_is_equal(scaled_isom1, scaled_isom2);
}
