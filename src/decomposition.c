#include <assert.h>
#include <flint/fmpq_mat.h>

#include "arith.h"
#include "decomposition.h"
#include "fmpq_poly.h"
#include "genus.h"
#include "hecke.h"

void decomposition_init(decomposition_t decomp, slong num)
{
  decomp->num = (slong*)malloc(num*sizeof(slong));
  decomp->bases = (fmpq_mat_t**)malloc(num*sizeof(fmpq_mat_t*));
  decomp->hecke = (fmpq_mat_t**)malloc(num*sizeof(fmpq_mat_t*));
  decomp->num_primes = 0;
  decomp->num_conductors = num;
  
  return;
}

// TODO - get rid of the unnecessary recursion
bool decomposition_finite_subspace(decomposition_t decomp, const genus_t genus, const fmpq_mat_t basis_V,
				   const int* ps, slong idx, slong num_ps, slong c)
{
  slong i, dim_V, next_idx;
  fmpq_mat_t T, fT, W;
  fmpq_poly_t f;
  fmpq_poly_factor_t fac;
  bool is_complete, is_complete_W;

#ifdef DEBUG
  printf("decomposing subspace with basis ");
  fmpq_mat_print(basis_V);
  printf("idx = %ld, ps[idx] = %d, num_ps = %ld\n", idx, ps[idx], num_ps);
#endif // DEBUG
  
  dim_V = fmpq_mat_nrows(basis_V);
  if (dim_V == 0) {
    return true;
  }

  if (idx >= num_ps) {
    (decomp->num[c])++;
    decomp->bases[c] = (fmpq_mat_t*)realloc(decomp->bases[c], (decomp->num[c])*sizeof(fmpq_mat_t));
    fmpq_mat_init(decomp->bases[c][decomp->num[c]-1], dim_V, dim_V);
    fmpq_mat_one(decomp->bases[c][decomp->num[c]-1]);
    return false;
  }

  fmpq_mat_init(T, dim_V, dim_V);
  fmpq_mat_init(fT, dim_V, dim_V);
  fmpq_poly_init(f);
#ifdef DEBUG
  printf("idx = %ld\n",idx);
  printf("decomp->hecke[idx] = %lx\n", (unsigned long)decomp->hecke[idx]);
#endif // DEBUG
  if (idx >= decomp->num_primes) {
#ifdef DEBUG
    printf("computing more hecke operators\n");
#endif // DEBUG
    (decomp->num_primes)++;
    decomp->hecke[idx] = (fmpq_mat_t*)malloc((genus->num_conductors)*sizeof(fmpq_mat_t));
    get_hecke_fmpq_mat_all_conductors(decomp->hecke[idx], genus, ps[idx], 1);
  }
#ifdef DEBUG
  printf("idx = %ld, hecke operator is", idx);
  fmpq_mat_print(decomp->hecke[idx][c]);
  printf("restricting hecke operator to space\n");
#endif // DEBUG
  restrict_mat(T, decomp->hecke[idx][c], basis_V);
#ifdef DEBUG
  printf("restricted matrix is");
  fmpq_mat_print(T);
#endif // DEBUG
  
  fmpq_mat_charpoly(f, T);
  fmpq_poly_factor(fac, f);

  is_complete = true;
  for (i = 0; i < fac->num; i++) {
    fmpq_poly_evaluate_fmpq_mat(fT, &(fac->p[i]), T);
    kernel_on(W, fT, basis_V);
    assert(W->rows != 0);
    // optimally, we will already compute the eigenvectors at this point.
    if (fac->exp[i] == 1) { // test for irreduciblity
      (decomp->num[c])++;
      decomp->bases[c] = (fmpq_mat_t*)realloc(decomp->bases[c], (decomp->num[c])*sizeof(fmpq_mat_t));
      fmpq_mat_init_set(decomp->bases[c][decomp->num[c]-1], W);
      is_complete_W = true;
    }
    else {
      next_idx = (fmpq_mat_nrows(W) == dim_V) ? idx+1 : 0;
      is_complete_W = decomposition_finite_subspace(decomp, genus, W, ps, next_idx, num_ps, c);
    }
    is_complete = (is_complete) && (is_complete_W);
    fmpq_mat_clear(W);
  }
  
  fmpq_mat_clear(fT);
  fmpq_mat_clear(T);
  fmpq_poly_factor_clear(fac);
  fmpq_poly_clear(f);
  
  return is_complete;
}

bool decomposition_finite(decomposition_t decomp, const genus_t genus, const int* ps,
			  slong num_ps, slong c)
{
  fmpq_mat_t basis_M;
  slong dim;
  bool is_complete;
  
  decomp->num[c] = 0;
  decomp->bases[c] = NULL;

  if (num_ps == 0)
    dim = 0;
  else
    dim = genus->dims[c];

  fmpq_mat_init(basis_M, dim, dim);
  fmpq_mat_one(basis_M);

  is_complete = decomposition_finite_subspace(decomp, genus, basis_M, ps, 0, num_ps, c);
  
  fmpq_mat_clear(basis_M);
  
  return is_complete;
}

void decompose(decomposition_t decomp, const genus_t genus, slong c)
{
  slong bound, num_ps;
  int* ps;
  bool is_complete;

  is_complete = false;
  bound = 10;
  
  while (!(is_complete)) {
    num_ps = primes_up_to(&ps, bound);
#ifdef DEBUG
    printf("trying to decompose with bound = %ld\n", bound);
#endif // DEBUG
    is_complete = decomposition_finite(decomp, genus, ps, num_ps, c);
    free(ps);
    bound *= 2;
  }

  return;
}

void decomposition_clear(decomposition_t decomp)
{
  slong c, i;
  
  for (c = 0; c < decomp->num_conductors; c++) {
    for (i = 0; i < decomp->num[c]; i++)
      fmpq_mat_clear(decomp->bases[c][i]);
    free(decomp->bases[c]);
  }
  free(decomp->bases);

  for (i = 0; i < decomp->num_primes; i++) {
    for (c = 0; c < decomp->num_conductors; c++)
      fmpq_mat_clear(decomp->hecke[i][c]);
    free(decomp->hecke[i]);
  }
    
  free(decomp->hecke);
  free(decomp->num);
  return;
}
