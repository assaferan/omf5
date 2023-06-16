/**********************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : decomposition.c
 *
 * Description: Decomposing M(Lambda) as a Hecke module
 *
 **********************************************************
 */

// System dependencies

#include <assert.h>

// Required packages dependencies

#include <flint/fmpq_mat.h>

// Self dependencies

#include "arith.h"
#include "decomposition.h"
#include "fmpq_poly.h"
#include "genus.h"
#include "hecke.h"

// initialization
void decomposition_init(decomposition_t decomp, slong num)
{
  decomp->num = (slong*)malloc(num*sizeof(slong));
  decomp->bases = (fmpq_mat_t**)malloc(num*sizeof(fmpq_mat_t*));
  // we should have num_primes here!!
  decomp->hecke = (fmpq_mat_t**)malloc(0*sizeof(fmpq_mat_t*));
  decomp->num_primes = 0;
  decomp->num_conductors = num;
  
  return;
}

// This function tries to decompose a subspace fixed by the Hecke algebra to irreducible subspaces,
// using a finite list of Hecke operators -
// Here basis_V is a matrix whose rows specify a basis for the subspace,
// ps is a list of primes, whose associated Hecke operators we wish to use,
// idx is the index of the prime we are currently using to decompose the space
// num_ps is the number of primes in ps,
// and c specifies the conductor
//
// The function proceeds recursively.
// It returns true if all the subspaces are irreducible as Hecke modules,
// and false if using all the primes in ps, there are still submodules which might be reducible.
// (due to failure of multiplicity one)
// TODO - get rid of the unnecessary recursion
bool decomposition_finite_subspace(decomposition_t decomp, const genus_t genus, const fmpq_mat_t basis_V,
				   const int* ps, slong idx, slong num_ps, slong c)
{
  slong i, dim_V, next_idx;
  fmpq_mat_t T, fT, W;
  fmpq_poly_t f;
  fmpq_poly_factor_t fac;
  bool is_complete, is_complete_W;
  
  dim_V = fmpq_mat_nrows(basis_V);
  if (dim_V == 0) {
    return true;
  }

  // If the index points beyond the length of our list of primes, we have failed to completely decompose the space
  // using these primes.
  if (idx >= num_ps) {
    // Wrapping up what we did compute - add the resulting basis to the appropriate place in decomp
    (decomp->num[c])++;
    decomp->bases[c] = (fmpq_mat_t*)realloc(decomp->bases[c], (decomp->num[c])*sizeof(fmpq_mat_t));
    fmpq_mat_init(decomp->bases[c][decomp->num[c]-1], dim_V, fmpq_mat_ncols(basis_V));
    fmpq_mat_set(decomp->bases[c][decomp->num[c]-1], basis_V);
    return false;
  }

#ifdef DEBUG
  fprintf(stderr, "decomposing subspace with basis ");
  fmpq_mat_print(basis_V);
  fprintf(stderr, "idx = %ld, ps[idx] = %d, num_ps = %ld\n", idx, ps[idx], num_ps);
#endif // DEBUG
  
  fmpq_mat_init(T, dim_V, dim_V);
  fmpq_mat_init(fT, dim_V, dim_V);
  fmpq_poly_init(f);
#ifdef DEBUG_LEVEL_FULL
  fprintf(stderr, "idx = %ld\n",idx);
  fprintf(stderr, "decomp->hecke[idx] = 0x%lx\n", (unsigned long)decomp->hecke[idx]);
#endif // DEBUG_LEVEL_FULL

  // If we have not yet computed the Hecke operators for this prime, we compute them now.
  // Note that with current implementation, computing the Hecke operator for all conductors at once saves time,
  // hence we compute already for all of them.
  if (idx >= decomp->num_primes) {
#ifdef DEBUG
    fprintf(stderr, "computing more hecke operators\n");
#endif // DEBUG
    (decomp->num_primes)++;
    decomp->hecke = (fmpq_mat_t**)realloc(decomp->hecke, (decomp->num_primes)*sizeof(fmpq_mat_t*));
    decomp->hecke[idx] = (fmpq_mat_t*)malloc((genus->num_conductors)*sizeof(fmpq_mat_t));
    get_hecke_fmpq_mat_all_conductors(decomp->hecke[idx], genus, ps[idx], 1);
  }
#ifdef DEBUG
  fprintf(stderr, "idx = %ld, hecke operator is", idx);
  fmpq_mat_print(decomp->hecke[idx][c]);
  fprintf(stderr, "restricting hecke operator to space\n");
#endif // DEBUG
  // T is the restriction of the Hecke operator to the relevant subspace
  restrict_mat(T, decomp->hecke[idx][c], basis_V);
#ifdef DEBUG
  fprintf(stderr, "restricted matrix is");
  fmpq_mat_print(T);
#endif // DEBUG

  // We factor the characteristic polynomials to obtain the subspaces
  // which are stable by its action
  fmpq_mat_charpoly(f, T);
  fmpq_poly_factor(fac, f);

  is_complete = true;
  for (i = 0; i < fac->num; i++) {
    fmpq_poly_evaluate_fmpq_mat(fT, &(fac->p[i]), T);
    kernel_on(W, fT, basis_V);
    // W is the stable subspace corresponding to this irreducible factor 
    assert(W->rows != 0);
    // We lift the basis to vectors of the original space
    // fmpq_mat_mul(W, W, basis_V);
    // optimally, we will already compute the eigenvectors at this point.
    if (fac->exp[i] == 1) { // test for irreduciblity
      (decomp->num[c])++;
      decomp->bases[c] = (fmpq_mat_t*)realloc(decomp->bases[c], (decomp->num[c])*sizeof(fmpq_mat_t));
      fmpq_mat_init_set(decomp->bases[c][decomp->num[c]-1], W);
      is_complete_W = true;
    }
    else { // if W = V is not irreducible, we move on to the next prime. If W != V, we can recursively decompose W.
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

// This function tries to decompose the space of orthogonal modular forms
// using a finite list of primes, specifiying the relevant Hecke operators -
// genus is the genus of the lattice,
// ps is the list of primes,
// num_ps their number,
// and c specifies the conductor (twist of the spinor norm)
// It returns true if all the subspaces in the decomposition are irreducible Hecke modules,
// and false if after using all primes there are still modules which mght be reducible.
// The decomposition is stored in decomp.
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

  // We use the previous recursive function, starting at the first prime, and with our
  // subspace being the entire space
  is_complete = decomposition_finite_subspace(decomp, genus, basis_M, ps, 0, num_ps, c);
  
  fmpq_mat_clear(basis_M);
  
  return is_complete;
}

// Finally, we may use the decomposition bases on a finite list of primes
// to write a function that fully decomposes the space.
// In theory, we would like to proceed until the decomposition is complete.
// However, due to oldforms, we may have a failure of multiplicity one,
// resulting in an infinite recursion.
// Until we implement degeneracy maps to overcome this issue,
// we fix a bound on the primes used, whcih should be sufficient for
// our purposes.

void decompose(decomposition_t decomp, const genus_t genus, slong c)
{
  slong /*bound, */ num_ps; /* , prime_idx, idx, p, j;*/
  int* ps;
  //  bool is_complete;
  // !! TODO - figure out degeneracy maps. Until we do, we use a fixed bound.
  slong MAX_BOUND = 20;

  // is_complete = false;
  // bound = 10;
  
  // while (!(is_complete) && (bound <= MAX_BOUND)) {
  num_ps = primes_up_to_prime_to(&ps, /* bound */ MAX_BOUND, fmpz_get_si(genus->disc));
#ifdef DEBUG
  fprintf(stderr, "trying to decompose with bound = %ld\n", /*bound*/ MAX_BOUND);
#endif // DEBUG
  /*is_complete = */ decomposition_finite(decomp, genus, ps, num_ps, c);
  free(ps);
  // bound *= 2;
  // }

  return;
}

// clear memory allocared for the decomposition structure
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
