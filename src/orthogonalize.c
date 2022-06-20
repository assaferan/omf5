#include "orthogonalize.h"

// This is somewhat of a duplicate for cholesky,
// but this one keeps everything integral.
// Maybe replace by cholesky from flint?
void orthogonalize_gram(fmpz_mat_t D, const fmpz_mat_t Q)
{
  fmpz_mat_t L;
  fmpz_t prod_diag, d, inner_sum, tmp_prod;
  int i, j, k, r;
  slong n;

  fmpz_init(d);
  fmpz_init(prod_diag);
  fmpz_init(inner_sum);
  fmpz_init(tmp_prod);
  n = fmpz_mat_nrows(Q);

  fmpz_mat_init(L, n, n);
 
  fmpz_one(prod_diag);
  fmpz_zero(d);
  fmpz_zero(inner_sum);
  
  // This works but inefficiently - for some reason we get O(n^4) operations.
  // !! TODO - check it out later
  // Oh I see - we should do the L update in two passes...
  for (i = 0; i < n; i++)
    {
      fmpz_set(fmpz_mat_entry(L, i, i), prod_diag);
      fmpz_set(d, prod_diag);
      for (j = 0; j < i; j++)
	{
	  fmpz_zero(fmpz_mat_entry(L, i, j)); 
	  for (k = j; k < i; k++)
	    {
	      fmpz_zero(inner_sum);
	      for (r = 0; r <= k; r++) {
		fmpz_mul(tmp_prod, fmpz_mat_entry(Q, i, r), fmpz_mat_entry(L, k, j));
		fmpz_addmul(inner_sum, fmpz_mat_entry(L, k, r), tmp_prod);
	      }
	      fmpz_mul(inner_sum, inner_sum, fmpz_mat_entry(L,i,i));
	      fmpz_divexact(inner_sum, inner_sum, fmpz_mat_entry(D,0,k));
	      fmpz_neg(inner_sum, inner_sum);
	      fmpz_add(fmpz_mat_entry(L, i, j), fmpz_mat_entry(L, i, j), inner_sum); 
	    }
	  fmpz_gcd(d, d, fmpz_mat_entry(L, i, j));
	}
      for (j = 0; j <= i; j++) {
	fmpz_divexact(fmpz_mat_entry(L, i, j), fmpz_mat_entry(L, i, j), d);
      }
      fmpz_zero(fmpz_mat_entry(D, 0, i));
      for (j = 0; j <= i; j++)
	for ( k = 0; k <= i; k++) {
	  fmpz_mul(tmp_prod, fmpz_mat_entry(Q, j, k), fmpz_mat_entry(L, i, k));
	  fmpz_addmul(fmpz_mat_entry(D, 0, i), fmpz_mat_entry(L, i, j), tmp_prod);
	}
      
      fmpz_lcm(prod_diag, prod_diag, fmpz_mat_entry(D, 0, i));
      for (j = i+1; j < n; j++)
	fmpz_zero(fmpz_mat_entry(L, i, j));
    }

  // Recall that this is an even lattice, so all entries in D
  // are even, and we are more interested in their half values,
  // which corresponds to the quadratic form.
  for ( i = 0; i < n; i++)
    fmpz_fdiv_q_2exp(fmpz_mat_entry(D, 0, i), fmpz_mat_entry(D, 0, i), 1);
  
#ifdef DEBUG_LEVEL_FULL
  printf("L=\n");
  fmpz_mat_print_pretty(L);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  fmpz_mat_clear(L);

  fmpz_clear(d);
  fmpz_clear(prod_diag);
  fmpz_clear(inner_sum);
  fmpz_clear(tmp_prod);
  return;
}
