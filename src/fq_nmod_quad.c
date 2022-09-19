#include "flint/fq_nmod.h"
#include "flint/fq_nmod_mat.h"
#include "flint/fq_nmod_poly.h"

#include "fq_nmod_quad.h"
#include "fq_nmod_mat.h"

#include "typedefs.h"

// This file implements utilities for quadratic forms over F_p

// return q(vec) 
void fq_nmod_quad_evaluate_p2(fq_nmod_t value, const fq_nmod_mat_t q, const fq_nmod_mat_t vec, const fq_nmod_ctx_t F)
{
  slong i, j;
  fq_nmod_t prod;

  fq_nmod_init(prod, F);
  fq_nmod_zero(value, F);
  
  for (i = 0; i < fq_nmod_mat_nrows(q); i++) {
    fq_nmod_mul(prod, fq_nmod_mat_entry(q, i, i), fq_nmod_mat_entry(vec, 0, i), F);
    fq_nmod_add(value, value, prod, F);
    for (j = i + 1; j < fq_nmod_mat_nrows(q); j++) {
      fq_nmod_mul(prod, fq_nmod_mat_entry(vec, 0, i), fq_nmod_mat_entry(vec, 0, j), F);
      fq_nmod_mul(prod, prod, fq_nmod_mat_entry(q, i, j), F);
      fq_nmod_add(value, value, prod, F);
    }
  }

  fq_nmod_clear(prod, F);
  return;
}

// return q(vec) 
void fq_nmod_quad_evaluate(fq_nmod_t value, const fq_nmod_mat_t q, const fq_nmod_mat_t vec, const fq_nmod_ctx_t F)
{
  fq_nmod_mat_t vec_t, vec_q, vec_q_vec_t;
  fq_nmod_t two;

  fq_nmod_init(two, F);
  fq_nmod_set_si(two, 2, F);

  if (fq_nmod_is_zero(two, F))
    return fq_nmod_quad_evaluate_p2(value, q, vec, F);
  
  fq_nmod_mat_init(vec_t, fq_nmod_mat_ncols(vec), fq_nmod_mat_nrows(vec), F);
  fq_nmod_mat_transpose(vec_t, vec, F);

  fq_nmod_mat_init(vec_q, fq_nmod_mat_nrows(vec), fq_nmod_mat_ncols(q), F);
  fq_nmod_mat_mul(vec_q, vec, q, F);

  fq_nmod_mat_init(vec_q_vec_t, fq_nmod_mat_nrows(vec), fq_nmod_mat_nrows(vec));
  fq_nmod_mat_mul(vec_q_vec_t, vec_q, vec_t, F);

  fq_nmod_set(value, fq_nmod_mat_entry(vec_q_vec_t, 0, 0), F);

  fq_nmod_div(value, value, two, F);
  
  fq_nmod_mat_clear(vec_q_vec_t, F);
  fq_nmod_mat_clear(vec_q, F);
  fq_nmod_mat_clear(vec_t, F);
  fq_nmod_clear(two, F);
  return;
}

void fq_nmod_quad_transform(fq_nmod_mat_t new_gram, const fq_nmod_mat_t gram, const fq_nmod_mat_t isom, const fq_nmod_ctx_t F)
{
  fq_nmod_mat_t isom_gram, isom_t;

  fq_nmod_mat_init(isom_gram, N, N, F);
  fq_nmod_mat_init(isom, N, N, F);
  
  fq_nmod_mat_transpose(isom_t, isom, F);
  fq_nmod_mat_mul(isom_gram, isom, gram, F);
  fq_nmod_mat_mul(new_gram, isom_gram, isom_t, F);

  fq_nmod_mat_clear(isom_gram, F);
  fq_nmod_mat_clear(isom_t, F);
  
  return;
}

bool fq_nmod_quad_isotropic_vector_p2(fq_nmod_mat_t vec, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, slong start)
{
  slong i,j,dim;
  fq_nmod_t a,b,c,d,e,f,g;

  for (j = start; j < n; j++) {
    if (fq_nmod_is_zero(fq_nmod_mat_entry(q,j,j),F)) {
      fq_nmod_one(fq_nmod_mat_entry(vec,0,j),F);
      return true;
    }
  }

  if (start + 2 == N) {
    fq_nmod_init(a,F);
    fq_nmod_init(b,F);
    fq_nmod_init(c,F);

    fq_nmod_set(a, fq_nmod_mat_entry(q,start,start), F);
    fq_nmod_set(b, fq_nmod_mat_entry(q,start,start+1),F);
    fq_nmod_set(c, fq_nmod_mat_entry(q,start+1,start+1),F);

    if (fq_nmod_is_zero(b,F)) {
      fq_nmod_div(b,a,c,F);
      assert(fq_nmod_is_square(b,F));
      fq_nmod_sqrt(b,b,F);
      fq_nmod_one(fq_nmod_mat_entry(vec,0,start),F);
      fq_nmod_set(fq_nmod_mat_entry(vec,0,start+1),b,F);
      fq_nmod_clear(c,F);
      fq_nmod_clear(b,F);
      fq_nmod_clear(a,F);
      return true;
    }

    // if we are here we should have a = b = c = 1, which is anisotropic
    assert(fq_nmod_is_one(a,F));
    assert(fq_nmod_is_one(b,F));
    assert(fq_nmod_is_one(c,F));
    
    fq_nmod_clear(c,F);
    fq_nmod_clear(b,F);
    fq_nmod_clear(a,F);
    return false;
  }

  // If we can find a pair of orthogonal basis vectors,
  //  we can easily construct an isotropic vector.
  for (i = start; i < N-1; i++)
    for (j = i+1; j < N; j++) {
      if (fq_nmod_is_zero(fq_nmod_mat_entry(q,i,j),F)) {
	fq_nmod_init(b,F);
	fq_nmod_div(b,fq_nmod_mat_entry(q,j,j),fq_nmod_mat_entry(q,i,i),F);
	assert(fq_nmod_is_square(b,F));
	fq_nmod_sqrt(b,b,F);
	fq_nmod_set(fq_nmod_mat_entry(vec,0,i),b,F);
	fq_nmod_one(fq_nmod_mat_entry(vec,0,j),F);
	fq_nmod_clear(b,F);
	return true;
      }
    }

  dim = N - start;
  if (dim == 1) {
    if (fq_nmod_is_zero(fq_nmod_mat_entry(q,start,start),F)) {
      fq_nmod_one(fq_nmod_mat_entry(vec,0,start),F);
      return true;
    }
    return false;
  }

  if (dim == 2) {
    for (i = 0; i < 2; i++)
      if (fq_nmod_is_zero(fq_nmod_mat_entry(q,start+i,start+i))) {
	fq_nmod_one(fq_nmod_mat_entry(vec,0,start+i),F);
	return true;
      }
    fq_nmod_init(a,F);
    fq_nmod_init(b,F);
    fq_nmod_init(c,F);

    fq_nmod_set(a,fq_nmod_mat_entry(q,start,start),F);
    fq_nmod_set(b,fq_nmod_mat_entry(q,start,start+1),F);
    fq_nmod_set(c,fq_nmod_mat_entry(q,start+1,start+1),F);
    if (fq_nmod_is_zero(b,F)) {
      fq_nmod_one(fq_nmod_mat_entry(vec,0,start),F);
      fq_nmod_one(fq_nmod_mat_entry(vec,0,start+1),F);
      fq_nmod_clear(c,F);
      fq_nmod_clear(b,F);
      fq_nmod_clear(a,F);
      return true;
    }
    // In this case a = b = c = 1, so the form is anisotropic
    
    fq_nmod_clear(c,F);
    fq_nmod_clear(b,F);
    fq_nmod_clear(a,F);
    return false;
  }
  assert(dim >= 3);
  // Otherwise, while the formulation is a bit more
  //  complicated, we can produce an isotropic vector
  //  by taking a linear combination of the first three
  //  basis vectors as follows:

  // Convenient references.
  fq_nmod_init(a,F);
  fq_nmod_init(b,F);
  fq_nmod_init(c,F);
  fq_nmod_init(d,F);
  fq_nmod_init(e,F);
  fq_nmod_init(f,F);
  
  fq_nmod_set(a,fq_nmod_mat_entry(q,start,start),F);
  fq_nmod_set(b,fq_nmod_mat_entry(q,start+1,start+1),F);
  fq_nmod_set(c,fq_nmod_mat_entry(q,start+2,start+2),F);
  fq_nmod_set(d,fq_nmod_mat_entry(q,start+1,start+2),F);
  fq_nmod_set(e,fq_nmod_mat_entry(q,start,start+2),F);
  fq_nmod_set(f,fq_nmod_mat_entry(q,start,start+1),F);

  // g = (b*e*e/f/f + c + e*d/f)/a in stages:
  
  // 1. g = e*d/f 
  fq_nmod_mul(g,d,e,F);
  fq_nmod_div(g,g,f,F);

  // 2. save d = e/f for later
  fq_nmod_div(d,e,f,F);

  // 3. b = b*e^2/f^2
  fq_nmod_sqr(f,f,F);
  fq_nmod_sqr(e,e,F);
  fq_nmod_mul(b,b,e,F);
  fq_nmod_div(b,b,f,F);

  // 4. g = (g+b+c)/a
  fq_nmod_add(g,g,b,F);
  fq_nmod_add(g,g,c,F);
  fq_nmod_div(g,g,a,F);

  assert(fq_nmod_is_square(g,F));

  fq_nmod_set(fq_nmod_mat_entry(vec,0,start),g,F);
  fq_nmod_set(fq_nmod_mat_entry(vec,0,start+1),d,F);
  fq_nmod_one(fq_nmod_mat_entry(vec,0,start+2),F);
  
  fq_nmod_clear(f,F);
  fq_nmod_clear(e,F);
  fq_nmod_clear(d,F);
  fq_nmod_clear(c,F);
  fq_nmod_clear(b,F);
  fq_nmod_clear(a,F);

  return true;
}

bool fq_nmod_quad_isotropic_vector(fq_nmod_mat_t vec, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, slong start, bool deterministic)
{
  fq_nmod_mat_t sub_B, rad, basis, subM, v;
  slong row, col, i, j, k, dim;
  fq_nmod_t a, b, c, d, ac, scalar;
  fmpz_t p, x, y;
  flint_rand_t state;
 
  fq_nmod_mat_init(sub_B,N-start,N-start,F);
  for (row = 0; row < N - start; row++)
    for (col = 0; col < N - start; col++)
      fq_nmod_set(fq_nmod_mat_entry(sub_B, row, col), fq_nmod_mat_entry(q, start+row,start+col),F);

  if ( (fmpz_cmp_si(fq_nmod_ctx_prime(F),2) != 0 ) && (fq_nmod_mat_rref(q,F) < N) ) {
    fq_nmod_mat_kernel(rad, sub_B, F);
    for (i = 0; i < N - start; i++)
      fq_nmod_set(fq_nmod_mat_entry(vec,0,start+i), fq_nmod_mat_entry(rad,0,i), F);
    fq_nmod_mat_clear(rad);
    fq_nmod_mat_clear(sub_B);
    return true;
  }

  dim = N - start;

  if (dim == 1) {
    if (fq_nmod_is_zero(fq_nmod_mat_entry(q,start,start),F)) {
      fq_nmod_one(fq_nmod_mat_entry(vec, 0, start), F);
      fq_nmod_mat_clear(sub_B);
      return true;
    }
    fq_nmod_mat_clear(sub_B);
    return false;
  }

  if (fmpz_cmp_si(fq_nmod_ctx_prime(F),2) == 0) {
    fq_nmod_mat_clear(sub_B);
    return fq_nmod_quad_isotropic_vector_p2(vec, q, F, start);
  }

  if (dim == 2) {
    fq_nmod_init(a, F);
    fq_nmod_init(b, F);
    fq_nmod_init(c, F);
    fq_nmod_init(d, F);
    fq_nmod_init(ac, F);

    fq_nmod_set(a, fq_nmod_mat_entry(q, start, start), F);
    fq_nmod_set(b, fq_nmod_mat_entry(q, start, start+1), F);
    fq_nmod_set(c, fq_nmod_mat_entry(q, start+1, start+1), F);
    // The form is isotropic if and only if b^2-ac is a square.
    // d = b^2 - ac
    fq_nmod_sqr(d, b, F);
    fq_nmod_mul(ac, a, c, F);
    fq_nmod_sub(d, d, ac, F);
    // If not a square, this form is anisotropic.
    if (!(is_square(d_lift, F))) {
      fq_nmod_mat_clear(sub_B);
      fq_nmod_clear(a, F);
      fq_nmod_clear(b, F);
      fq_nmod_clear(c, F);
      fq_nmod_clear(d, F);
      fq_nmod_clear(ac, F);
      return false;
    }
    // Since a ne 0 and the form is isotropic, we're done.
    fq_nmod_sqrt(d,d,F);
    // vec[start] = -((b+d)/a)
    fq_nmod_add(d,d,b,F);
    fq_nmod_div(d,d,a,F);
    fq_nmod_neg(d,d,F);
    fq_nmod_set(fq_nmod_mat_entry(vec,0,start),d,F);
    fq_nmod_one(fq_nmod_mat_entry(vec,0,start+1),F);
    
    fq_nmod_mat_clear(sub_B);
    fq_nmod_clear(a, F);
    fq_nmod_clear(b, F);
    fq_nmod_clear(c, F);
    fq_nmod_clear(d, F);
    fq_nmod_clear(ac, F);
    return true;
  }

  assert(dim >= 3);

  // Check the diagonal
  for (i = start; i < n; i++)
    if (fq_nmod_is_zero(fq_nmod_mat_entry(q,i,i),F)) {
      fq_nmod_one(fq_nmod_mat_entry(vec,0,i),F);
      fq_nmod_mat_clear(sub_B);
      return true;
    }

  // isometry on the submatrix of 3 first variables
  fq_nmod_mat_init(basis,3,3,F);
  fq_nmod_mat_init(subM,3,3,F);

  fq_nmod_mat_zero(basis,F);
  for (i = 0; i < 3; i++)
    fq_nmod_one(fq_nmod_mat_entry(basis,i,i),F);
  
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      fq_nmod_set(fq_nmod_mat_entry(subM,i,j), fq_nmod_mat_entry(q,start+i,start+j),F);

  fq_nmod_init(scalar,F);

  // clear the off-diagonal entries
  for (i = 0; i < 2; i++)
    for (j = i+1; j < 3; j++) {
      fq_nmod_div(scalar, fq_nmod_mat_entry(subM,i,j), fq_nmod_mat_entry(subM,i,i), F);
      fq_nmod_neg(scalar,scalar,F);
      fq_nmod_mat_add_col(subM,j,i,scalar,F);
      fq_nmod_mat_add_row(subM,j,i,scalar,F);
      fq_nmod_mat_add_row(basis,j,i,scalar,F);
      if (fq_nmod_is_zero(fq_nmod_mat_entry(subM,j,j),F)) {
	for (k = 0; k < 3; k++)
	  fq_nmod_set(fq_nmod_mat_entry(vec,0,start+k), fq_nmod_mat_entry(basis,j,k),F);
	fq_nmod_clear(scalar,F);
	fq_nmod_mat_clear(subM,F);
	fq_nmod_mat_clear(basis,F);
	return true;
      }
    }

  // Check if the first two variables alone are isotropic.
  fq_nmod_init(d,F);

  fq_nmod_mul(d, fq_nmod_mat_entry(subM,0,0), fq_nmod_mat_entry(subM,1,1),F);
  fq_nmod_neg(d,d,F);

  if (fq_nmod_is_square(d,F)) {
    fq_nmod_sqrt(d,d,F);
    fq_nmod_init(a,F);
    for (k = 0; k < 3; k++) {
      fq_nmod_div(a,fq_nmod_mat_entry(q,start,start),d,F);
      fq_nmod_mul(a, a, fq_nmod_mat_entry(basis,1,k), F);
      fq_nmod_add(fq_nmod_mat_entry(vec,0,start+k), a, fq_nmod_mat_entry(basis,0,k), F);
    }
    fq_nmod_clear(a,F);
    fq_nmod_clear(d,F);
    fq_nmod_clear(scalar,F);
    fq_nmod_mat_clear(subM,F);
    fq_nmod_mat_clear(basis,F);
    return true;
  }

  if (deterministic) {
    fmpz_init(p);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_set(p, fq_nmod_ctx_prime(F));
    // The quadratic form over three variables.
    fq_nmod_mat_init(v,1,3,F);

    for (fmpz_zero(x); !(fmpz_equal(x,p)); fmpz_add_si(x,x,1))
      for (fmpz_zero(y); !(fmpz_equal(y,p)); fmpz_add_si(y,y,1)) {
	fq_nmod_set_fmpz(fq_nmod_mat_entry(v,0,0),x,F);
	fq_nmod_set_fmpz(fq_nmod_mat_entry(v,0,1),y,F);
	fq_nmod_one(fq_nmod_mat_entry(v,0,2),F);
	fq_nmod_quad_evaluate(d,subM,v,F);
	if (fq_nmod_is_zero(d,F)) {
	  // Found an isotropic vector, return it.
	  for (j = 0; j < 3; j++) {
	    // vec[start+j] = v[0]*basis(0,j) + v[1]*basis(1,j) + basis(2,j);
	    fq_nmod_mul(fq_nmod_mat_entry(vec,0,start+j),fq_nmod_mat_entry(v,0,0),fq_nmod_mat_entry(basis,0,j),F);
	    fq_nmod_mul(d,fq_nmod_mat_entry(v,0,1),fq_nmod_mat_entry(basis,1,j),F);
	    fq_nmod_add(fq_nmod_mat_entry(vec,0,start+j),fq_nmod_mat_entry(vec,0,start+j),d,F);
	    fq_nmod_add(fq_nmod_mat_entry(vec,0,start+j),fq_nmod_mat_entry(vec,0,start+j),fq_nmod_mat_entry(basis,2,j),F);
	  }
	  fq_nmod_mat_clear(v,F);
	  fmpz_clear(y);
	  fmpz_clear(x);
	  fmpz_clear(p);
	  fq_nmod_clear(d,F);
	  fq_nmod_clear(scalar,F);
	  fq_nmod_mat_clear(subM,F);
	  fq_nmod_mat_clear(basis,F);
	  return true;
	}
      }
    
    fq_nmod_mat_clear(v,F);
    fmpz_clear(y);
    fmpz_clear(x);
    fmpz_clear(p);
  }

  // If we're fine with a probabilitistic means of finding
  //  isotropic vectors, we can find them much faster.

  fq_nmod_init(a,F);
  fq_nmod_init(b,F);
  fq_nmod_init(c,F);

  fq_nmod_set(a,fq_nmod_mat_entry(subM,0,0),F);
  fq_nmod_set(b,fq_nmod_mat_entry(subM,1,1),F);
  fq_nmod_set(c,fq_nmod_mat_entry(subM,2,2),F);

  fq_nmod_mat_init(v,1,2,F);
  
  do {
    do {
      do {
	for (i = 0; i < 2; i++)
	  fq_nmod_randtest(fq_nmod_mat_entry(v,0,i),state,F);
      } while(fq_nmod_is_zero(fq_nmod_mat_entry(v,0,0),F) && fq_nmod_is_zero(fq_nmod_mat_entry(v,0,1),F));
      // d = -(a*v[0]*v[0] + b*v[1]*v[1])/c;
      fq_nmod_sqr(fq_nmod_mat_entry(v,0,0),fq_nmod_mat_entry(v,0,0),F);
      fq_nmod_sqr(fq_nmod_mat_entry(v,0,1),fq_nmod_mat_entry(v,0,1),F);
      fq_nmod_mul(a,a,fq_nmod_mat_entry(v,0,0),F);
      fq_nmod_mul(b,b,fq_nmod_mat_entry(v,0,1),F);
      fq_nmod_add(d,a,b,F);
      fq_nmod_div(d,d,c,F);
      fq_nmod_neg(d,d,F);
    } while(!fq_nmod_is_square(d,F));

    fq_nmod_sqrt(d,d,F);
    nonzero = false;
    for (j = 0; j < 3; j++) {
      // vec[start+j] = v[0]*basis(0,j) + v[1]*basis(1,j) + basis(2,j);
      fq_nmod_mul(fq_nmod_mat_entry(vec,0,start+j),fq_nmod_mat_entry(v,0,0),fq_nmod_mat_entry(basis,0,j),F);
      fq_nmod_mul(d,fq_nmod_mat_entry(v,0,1),fq_nmod_mat_entry(basis,1,j),F);
      fq_nmod_add(fq_nmod_mat_entry(vec,0,start+j),fq_nmod_mat_entry(vec,0,start+j),d,F);
      fq_nmod_add(fq_nmod_mat_entry(vec,0,start+j),fq_nmod_mat_entry(vec,0,start+j),fq_nmod_mat_entry(basis,2,j),F);
    }
  } while(!nonzero);
  
  fq_nmod_mat_clear(v,F);
  
  fq_nmod_clear(c,F);
  fq_nmod_clear(b,F);
  fq_nmod_clear(a,F);
  
  fq_nmod_clear(d,F);
  fq_nmod_clear(scalar,F);
  fq_nmod_mat_clear(subM,F);
  fq_nmod_mat_clear(basis,F);

  return true;
}

void fq_nmod_quad_split_hyperbolic_plane(const fq_nmod_mat_t vec, fq_nmod_mat_t gram, fq_nmod_mat_t basis, const fq_nmod_mat_t q, const fq_nmod_ctx_t F,  slong start)
{
  fq_nmod_mat_t original_gram, vec;
  slong i, pivot, idx;
  fmpz_t p;
  bool is_in_radical;
  fq_nmod_t scalar, two;

#ifdef DEBUG
  fq_nmod_mat_t tmp, basis_t;
#endif // DEBUG
  
  // The change of basis which preserving the isometry.
  fq_nmod_mat_zero(basis,F);
  for (i = 0; i < fq_nmod_mat_nrows(basis); i++)
    fq_nmod_one(fq_nmod_mat_entry(basis,i,i),F);

  // Make a copy of the Gram matrix.
  fq_nmod_mat_set(gram,q,F);

  fmpz_init_set(p, fq_nmod_ctx_prime(F));

  // Set the diagonal entries to zero when in characteristic 2.
  // This is because we are decomposing the associated bilinear form
  if (fmpz_equal_si(p,2)) {
    for (i = start; i < n; i++)
      fq_nmod_zero(fq_nmod_mat_entry(gram,i,i),F);
  }

  fq_nmod_mat_init_set(original_gram,gram,F);
  pivot = start;
  while (fq_nmod_is_zero(fq_nmod_mat_entry(vec,0,pivot),F)) pivot++;

  assert(pivot < N);

  // Perform the necessary basis changes so that vec becomes the first
  //  basis vector.

  fq_nmod_mat_mul_row(basis,pivot,fq_nmod_mat_entry(vec,0,pivot),F);
  fq_nmod_mat_mul_row(gram,pivot,fq_nmod_mat_entry(vec,0,pivot),F);
  fq_nmod_mat_mul_col(gram,pivot,fq_nmod_mat_entry(vec,0,pivot),F);
  for (i = pivot+1; i < N; i++) {
    fq_nmod_mat_add_row(basis,pivot,i,fq_nmod_mat_entry(vec,0,i),F);
    fq_nmod_mat_add_row(gram,pivot,i,fq_nmod_mat_entry(vec,0,i),F);
    fq_nmod_mat_add_col(gram,pivot,i,fq_nmod_mat_entry(vec,0,i),F);
  }
  fq_nmod_mat_swap_rows(basis,start,pivot,F);
  fq_nmod_mat_swap_rows(gram,start,pivot,F);
  fq_nmod_mat_swap_cols(gram,start,pivot,F);

  is_in_radical = true;

  // Find a basis vector which is not orthogonal to our isotropic vector.
  idx = start;
  for (; idx < n; idx++)
    // !! TODO - replace by is_zero_row
    if (!fq_nmod_is_zero(fq_nmod_mat_entry(gram,start,idx),F)) {
      is_in_radical = false;
      break;
    }
  
  // If the first row is entirely zero, then this vector belongs to the
  //  radical of the form.

  if (is_in_radical) {
    if (fmpz_equal_si(p,2)) {
      fq_nmod_mat_init(vec,1,N,F);
      // Recover the quadratic form along the diagonal.
      for (i = start; i < N; i++) {
	for (col = 0; col < N; col++)
	  fq_nmod_set(fq_nmod_mat_entry(vec,0,col),fq_nmod_mat_entry(basis,i,col),F);
	fq_nmod_quad_evaluate(fq_nmod_mat_entry(gram,i,i),q,vec,F);
      }
      fq_nmod_mat_clear(vec,F);
    }
    fq_nmod_mat_clear(original_gram,F);
    fmpz_clear(p);
    return;
  }
  
  // Swap this basis vector with the second basis vector.
  fq_nmod_mat_swap_rows(basis,start+1,idx,F);
  fq_nmod_mat_swap_rows(gram,start+1,idx,F);
  fq_nmod_mat_swap_cols(gram,start+1,idx,F);

  // Normalize the second basis vector so the (0,1)-entry is 1.
  fq_nmod_inv(scalar,fq_nmod_mat_entry(gram,start,start+1),F);
  fq_nmod_mat_mul_row(basis,start+1,scalar,F);
  fq_nmod_mat_mul_row(gram,start+1,scalar,F);
  fq_nmod_mat_mul_col(gram,start+1,scalar,F);

  // Determine the appropriate scalar for clearing out the (1,1)-entry.
  fq_nmod_init(scalar,F);
  if (fmpz_equal_si(p,2)) {
    fq_nmod_mat_init(vec,1,N,F);
    for (col = 0; col < N; col++)
      fq_nmod_set(fq_nmod_mat_entry(vec,0,col),fq_nmod_mat_entry(basis,start+1,col),F);
    fq_nmod_quad_evaluate(scalar,q,vec,F);
    fq_nmod_mat_clear(vec,F);
  }
  else {
    fq_nmod_init(two,F);
    fq_nmod_one(two,F);
    fq_nmod_add(two,two,two,F);

    fq_nmod_div(scalar,fq_nmod_mat_entry(gram,start+1,start+1),two,F);
    fq_nmod_neg(scalar,scalar,F);

    fq_nmod_clear(two,F);
  }

  // Clear the (1,1)-entry in the Gram matrix.
  fq_nmod_mat_add_row(basis,start+1,start,scalar,F);
  fq_nmod_mat_add_row(gram,start+1,start,scalar,F);
  fq_nmod_mat_add_col(gram,start+1,start,scalar,F);

  // Clear the remaining entries in the Gram matrix.
  for (i = start+2; i < N; i++) {
    // Clear first row/column.
    fq_nmod_neg(scalar,fq_nmod_mat_entry(gram,start,i),F);
    fq_nmod_mat_add_row(basis,i,start+1,scalar,F);
    fq_nmod_mat_add_row(gram,i,start+1,scalar,F);
    fq_nmod_mat_add_col(gram,i,start+1,scalar,F);

    // Clear second row/column.
    fq_nmod_neg(scalar,fq_nmod_mat_entry(gram,start+1,i),F);
    fq_nmod_mat_add_row(basis,i,start,scalar,F);
    fq_nmod_mat_add_row(gram,i,start,scalar,F);
    fq_nmod_mat_add_col(gram,i,start,scalar,F); 
  }

#ifdef DEBUG
  fq_nmod_mat_init(basis_t, fq_nmod_mat_ncols(basis), fq_nmod_mat_nrows(basis), F);
  fq_nmod_mat_init(tmp, fq_nmod_mat_nrows(basis), fq_nmod_mat_ncols(original_gram), F);
  fq_nmod_mat_transpose(basis_t, basis, F);
  fq_nmod_mat_mul(tmp, basis, original_gram, F);
  fq_nmod_mat_mul(tmp, tmp, basis_t, F);
  assert(fq_nmod_mat_equal(tmp, gram, F));
  fq_nmod_mat_clear(tmp, F);
  fq_nmod_mat_clear(basis_t, F);
#endif // DEBUG

  // In characteristic 2, we need to recover the diagonal entries by
  //  evaluating the basis via the quadratic form.

  if (fmpz_equal_si(p,2)) {
    for (i = start; i < N; i++) {
      fq_nmod_mat_init(vec,1,N,F);
      for (col = 0; col < N; col++)
	fq_nmod_set(fq_nmod_mat_entry(vec,0,col),fq_nmod_mat_entry(basis,i,col),F);
      fq_nmod_quad_evaluate(fq_nmod_mat_entry(gram,i,i),q,vec,F);
      fq_nmod_mat_clear(vec,F);
    }
  }
  
  fq_nmod_clear(scalar,F);
  fq_nmod_mat_clear(original_gram,F);
  fmpz_clear(p);
  return;
}

void fq_nmod_quad_hyperbolize(fq_nmod_mat_t gram, fq_nmod_mat_t basis, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, bool deterministic, slong start)
{
  fq_nmod_mat_t vec, original_gram;
  fq_nmod_t scalar;
  bool found;
  slong i, dim, lower_dim;
  bool is_square[2];

  fq_nmod_mat_init(vec,1,N,F);
  found = fq_nmod_quad_isotropic_vector(vec,q,F,start,deterministic);
  dim = N - start;

  // The space is anisotropic.
  if (!found) {
    fq_nmod_mat_init_set(original_gram,gram,F);
    fq_nmod_init(scalar,F);
    if (dim == 1) {
      // Check if the (0,0)-entry is a square.
      if (fq_nmod_is_square(fq_nmod_mat_entry(gram,start,start),F)) {
	// If so, make it a 1.
	fq_nmod_sqrt(scalar, fq_nmod_mat_entry(gram,start,start),F);
	fq_nmod_inv(scalar,scalar,F);
	fq_nmod_mat_mul_row(basis,start,scalar,F);
	fq_nmod_mat_mul_col(gram,start,scalar,F);
	fq_nmod_mat_mul_row(gram,start,scalar,F);
      }
      fq_nmod_clear(scalar,F);
      fq_nmod_mat_clear(original_gram,F);
      fq_nmod_mat_clear(vec,F);
      return;
    }
    if (fmpz_equal_si(fq_nmod_ctx_prime(F),2)) {
      // Make the (0,0)-entry equal to 1.
      assert(fq_nmod_is_square(fq_nmod_mat_entry(gram,start,start),F));
      fq_nmod_sqrt(scalar, fq_nmod_mat_entry(gram,start,start),F);
      fq_nmod_inv(scalar,scalar,F);
      fq_nmod_mat_mul_row(basis,start,scalar,F);
      fq_nmod_mat_mul_col(gram,start,scalar,F);
      fq_nmod_mat_mul_row(gram,start,scalar,F);

      // Make the (0,1)-entry equal to 1.
      fq_nmod_inv(scalar,fq_nmod_mat_entry(gram,start,start+1),F);
      fq_nmod_mat_mul_row(basis,start+1,scalar,F);
      fq_nmod_mat_mul_col(gram,start+1,scalar,F);
      fq_nmod_mat_mul_row(gram,start+1,scalar,F);
      
      fq_nmod_clear(scalar,F);
      fq_nmod_mat_clear(original_gram,F);
      fq_nmod_mat_clear(vec,F);
      return;
    }

    // Clear the (0,1)-entry.
    fq_nmod_div(scalar,fq_nmod_mat_entry(gram,start,start+1),fq_nmod_mat_entry(gram,start,start),F);
    fq_nmod_neg(scalar,scalar,F);
    fq_nmod_mat_add_row(basis,start+1,start,scalar,F);
    fq_nmod_mat_add_row(gram,start+1,start,scalar,F);
    fq_nmod_mat_add_col(gram,start+1,start,scalar,F);

    // If the (1,1)-entry is a square, make it the first entry.
    if (fq_nmod_is_square(fq_nmod_mat_entry(gram,start+1,start+1),F)) {
      fq_nmod_mat_swap_rows(basis,start,start+1,F);
      fq_nmod_mat_swap_rows(gram,start,start+1,F);
      fq_nmod_mat_swap_cols(gram,start,start+1,F);
    }

    for (i = 0; i < 2; i++) {
      // Check if the (i,i)-entry is a square then clear it, if so.
      is_square[i] = fq_nmod_is_square(fq_nmod_mat_entry(gram,start+i,start+i),F);
      if (is_square[i]) {
	fq_nmod_sqrt(scalar, fq_nmod_mat_entry(gram,start+i,start+i),F);
	fq_nmod_inv(scalar,scalar,F);
	fq_nmod_mat_mul_row(basis,start+i,scalar,F);
	fq_nmod_mat_mul_col(gram,start+i,scalar,F);
	fq_nmod_mat_mul_row(gram,start+i,scalar,F);
      }
    }

    // If neither are squares, make them -1 (note that this occurs
    //  if and only if -1 is not a square).
    if ((!is_square[0]) && (!is_square[1])) {
      for (i = 0; i < 2; i++) {
	fq_nmod_neg(scalar,fq_nmod_mat_entry(gram,start+i,start+i),F);
	assert(fq_nmod_is_square(scalar,F));
	fq_nmod_sqrt(scalar,scalar,F);
	fq_nmod_inv(scalar,scalar,F);
	fq_nmod_mat_mul_row(basis,start+i,scalar,F);
	fq_nmod_mat_mul_col(gram,start+i,scalar,F);
	fq_nmod_mat_mul_row(gram,start+i,scalar,F);
      }
    }
    
    fq_nmod_clear(scalar,F);
    fq_nmod_mat_clear(original_gram,F);
    return;
  }

  fq_nmod_init(scalar,F);
  fq_nmod_quad_evaluate(scalar, q, vec, F);
  assert(fq_nmod_is_zero(scalar,F));

  // Attempt to split a hyperbolic plane from the form.
  fq_nmod_quad_split_hyperbolic_plane(vec,gram,basis,q,F,start);

  // Determine how many dimensions we need to split off.
  lower_dim = is_zero_row(gram,0) ? 1 : 2;

  if (dim > lower_dim) {
    // Split the hyperbolic plane from the form.
    fq_nmod_mat_init_set(q_split,gram,F);
    // !! TODO - check maybe we have to replace basis here
    fq_nmod_mat_init(newbasis,N,N,F);
    fq_nmod_mat_zero(newbasis,F);
    for (i = 0; i < N; i++)
      fq_nmod_one(fq_nmod_mat_entry(newbasis,i,i),F);
    fq_nmod_quad_hyperbolize(gram,newbasis,q,F,deterministic,start+lower_dim);
    fq_nmod_mat_mul(basis,newbasis,basis,F);
    fq_nmod_mat_clear(newbasis,F);
    fq_nmod_mat_clear(q_split,F);
  }

#ifdef DEBUG_LEVEL_FULL
  printf("After hyperbolize_form with start = %d.\n", start);
  printf("Resulting gram matrix is \n");
  fq_nmod_mat_print_pretty(gram,F);
  printf("Resulting basis is \n");
  fq_nmod_mat_print_pretty(basis,F);
#endif // DEBUG_LEVEL_FULL
  
  fq_nmod_clear(scalar,F);
  fq_nmod_mat_clear(vec,F);
  return;
}

void fq_nmod_quad_decompose(fq_nmod_mat_t gram, fq_nmod_mat_t basis, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, bool deterministic)
{
  size_t rad, pos, i, j, upper;
  fq_nmod_mat_t basis_t;
  
#ifdef DEBUG_LEVEL_FULL
  fq_nmod_mat_t temp1, temp2, temp3, basis_i;
  fq_nmod_t value;

  fq_nmod_init(value);
  fmpq_mat_init(basis_i, 1, N, F);
  fmpq_mat_init(temp1, N, N, F);
  fmpq_mat_init(temp2, N, N, F);
  fmpq_mat_init(temp3, N, N, F);
#endif // DEBUG_LEVEL_FULL
  
  fq_nmod_mat_one(basis, F);
  fq_nmod_quad_hyperbolize(gram, basis, q, F, deterministic, 0);
  
#ifdef DEBUG_LEVEL_FULL
  printf("After hyperbolize_form.\n");
  printf("Resulting gram matrix is\n");
  fq_nmod_mat_print_pretty(gram, F);
  printf("Resulting basis is\n");
  fq_nmod_mat_print_pretty(basis, F);
  if (fmpz_cmp_si(fq_nmod_ctx_prime(F), 2) == 0)  { // p = 2
    // Verify that the basis evaluates correctly on the form.
    for (i = 0; i < N; i++) {
      // can use a window here, but let's just clone for now
      for (j = 0; j < N; j++)
	fq_nmod_set(fq_nmod_mat_entry(basis_i,0,j), fq_nmod_mat_entry(basis,i,j) ,F);
      fq_nmod_quad_evaluate(value, q, basis_i, F);
      assert(fq_nmod_equal(value, fq_nmod_mat_entry(gram,i,i), F));
    }
    // Zero out the diagonal to verify the bilinear form is correct.
    // Recall that in characteristic 2, the associated bilinear form
    // is b_q(x,y) = q(x+y) - q(x) - q(y), with zeros on the diagonal
    fq_nmod_mat_set(temp1, gram, F);
    fq_nmod_mat_set(temp2, q, F);
    for (i = 0; i < N; i++) {
      fq_nmod_zero(fq_nmod_mat_entry(temp1, i, i), F);
      fq_nmod_zero(fq_nmod_mat_entry(temp2, i, i), F);
    }
  }
  else {
    fq_nmod_mat_set(temp1, gram, F);
    fq_nmod_mat_set(temp2, q, F);
  }
  // Verify that the bilinear forms are isometric.
  fq_nmod_quad_transform(temp3,temp2,basis,F);
  assert(fq_nmod_mat_equal(temp3, temp1,F));
#endif // DEBUG_LEVEL_FULL

  // Let's bubble the basis vectors which belong to the radical to the
  //  end of the basis list.
  rad = 0;
  pos = N;
  while (pos >= 1) {
    if (is_zero_row(gram, pos-1)) {
      rad++;
      for (i = pos; i < N; i++) {
	fq_nmod_mat_swap_rows(basis,i-1,i,F);
	fq_nmod_mat_swap_rows(gram,i-1,i,F);
	fq_nmod_mat_swap_cols(gram,i-1,i,F);
      }
    }
    pos--;
  }

  // Let's put the hyperbolic planes in our standard antidiagonal form.

  // The upper index of the hyperbolic space.
  upper = n + 1 - rad;
  do {
    upper--;
  } while ((upper >= 1) && (!fq_nmod_is_zero(fq_nmod_mat_entry(gram,upper-1,upper-1),F)));

  // Indices of the basis vectors we'll be swapping.
  i = 1;
  j = upper;

  // Keep swapping basis vectors until j is less than or equal to i. Note
  //  that if there are no hyperbolic planes, this does nothing.
  while (i < j) {
    fq_nmod_mat_swap_rows(basis,i-1,j-2,F);
    fq_nmod_mat_swap_cols(gram,i-1,j-2,F);
    fq_nmod_mat_swap_rows(gram,i-1,j-2,F);
    i += 2;
    j -= 2;
  }
  // Since we did everything with row vectors, we need to transpose the
  //  basis, so that the subsequent code that utilizes it doesn't break.
  fq_nmod_mat_init(basis_t,N,N,F);
  fq_nmod_mat_transpose(basis_t, basis, F);
  fq_nmod_mat_set(basis, basis_t, F);
  
#ifdef DEBUG_LEVEL_FULL
  fq_nmod_clear(value);
  fq_nmod_mat_clear(temp1, F);
  fq_nmod_mat_clear(temp2, F);
  fq_nmod_mat_clear(temp3, F);
  fq_nmod_mat_clear(basis_i, F);
#endif // DEBUG_LEVEL_FULL
  return;
}

void fq_nmod_poly_set_fq_nmod_quad(fq_nmod_poly_t poly, const fq_nmod_mat_t q, const fq_nmod_ctx_t F, const fq_nmod_mpoly_ctx_t R)
{
  slong i,j,gen_idx;
  fq_nmod_mpoly_t mon;
  fq_nmod_t two;

  fq_nmod_mpoly_init(mon,R);

  fq_nmod_mpoly_zero(poly,R);
  gen_idx = 0;
  for (i = 0; i < N; i++)
    for (j = i; j < N; j++) {
      fq_nmod_mpoly_gen(mon,gen_idx,R);
      fq_nmod_mpoly_scalar_mul_fq_nmod(mon,mon,fq_mat_entry(q,i,j),R);
      fq_nmod_poly_add(poly,poly,mon,R);
      gen_idx++; 
    }

  if (!fmpz_equal_si(fq_ctx_prime(F),2)) {
    fq_nmod_init(two,F);
    fq_nmod_one(two,F);
    fq_nmod_add(two,two,two,F);
    for (i = 0; i < N; i++) {
      // the index of x_{i,i}
      gen_idx = i*N - i*(i+1)/2;
      fq_nmod_mpoly_gen(mon,gen_idx,R);
      fq_nmod_mpoly_scalar_mul_fq_nmod(mon,mon,fq_mat_entry(q,i,i),R);
      fq_nmod_mpoly_scalar_div_fq_nmod(mon,mon,two,R);
      fq_nmod_mpoly_sub(poly,poly,mon,R);
    }
    fq_nmod_clear(two,F);
  }

  fq_nmod_mpoly_clear(mon,R);
  
  return;
}
