#include <assert.h>

#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>

#include "arith.h"
#include "jordan.h"

void jordan_data_init(jordan_data_t jordan, ulong n)
{
  jordan->num_blocks = 0;
  jordan->matrices = (fmpq_mat_t*)malloc(n*sizeof(fmpq_mat_t));
  jordan->grams = (fmpq_mat_t*)malloc(n*sizeof(fmpq_mat_t));
  jordan->exponents = (ulong*)malloc(n*sizeof(ulong));

  return;
}

void jordan_data_clear(jordan_data_t jordan)
{
  ulong i;
  for (i = 0; i < jordan->num_blocks; i++) {
    fmpq_mat_clear(jordan->matrices[i]);
    fmpq_mat_clear(jordan->grams[i]);
  }
  free(jordan->matrices);
  free(jordan->grams);
  free(jordan->exponents);
  jordan->num_blocks = 0;
  return;
}

// compute the inner product of S[idx1] and S[idx2] with Gram matrix G, and store it in res
void inner_product(fmpq_t res, const fmpz_mat_t G,
		   const fmpq_mat_t S, ulong idx1, ulong idx2) 
{
  ulong i,j,n;
  fmpq_t mult;

  fmpq_init(mult);
  n = fmpz_mat_nrows(G);

  assert((idx1 < n) && (idx2 < n));

  fmpq_zero(res);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      fmpq_mul_fmpz(mult, fmpq_mat_entry(S, idx1, i), fmpz_mat_entry(G, i, j));
      fmpq_addmul(res, mult, fmpq_mat_entry(S, idx2, j));
    }
  
  fmpq_clear(mult);
  
  return;
}

void jordan_decomposition(jordan_data_t jordan, const matrix_TYP* q, const fmpz_t p)
{
  int even;
  fmpq_mat_t S, G, F, mF;
  fmpz_mat_t B;
  ulong i, j, k, l, n, old_val, ii, m, val, i1, i2, nrows, row, col;
  ulong* blocks;
  fmpq_mat_t T, X, m_t;
  fmpq_t d, tl, tl2, ul, ul2, nrm;

#ifdef DEBUG_LEVEL_FULL
  fmpq_t tmp_rat;
#endif // DEBUG_LEVEL_FULL

  n = q->rows;

  blocks = (ulong*)malloc(n*sizeof(ulong));
  
  fmpq_mat_init(S, n, n);
  fmpq_mat_init(G, n, n);
  fmpq_mat_init(F, n, n);
  fmpz_mat_init(B, n, n);
  fmpq_mat_init(X, 1, n);

  fmpq_mat_init(T,2,2);
  fmpq_init(d);
  fmpq_init(tl);
  fmpq_init(tl2);
  fmpq_init(ul);
  fmpq_init(ul2);
  fmpq_init(nrm);

  jordan->num_blocks = 0;
  even = fmpz_equal_si(p, 2);
  fmpq_mat_one(S);

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      fmpz_set_si(fmpz_mat_entry(B, i, j), q->array.SZ[i][j]);
      fmpq_set_fmpz(fmpq_mat_entry(F, i, j), fmpz_mat_entry(B, i, j));
    }

  k = 0;
  // virtually infinity
  old_val = ULONG_MAX;

  while (k < n) {
#ifdef DEBUG_LEVEL_FULL
    printf("k = %lu\n", k);
#endif // DEBUG_LEVEL_FULL
    // G = SFS^t
    // !! TODO - can we write simply G = S*(this->_B)*S.transpose() ?
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
	inner_product(fmpq_mat_entry(G,i,j), B, S, i, j);
	
#ifdef DEBUG_LEVEL_FULL
    printf("G = \n");
    fmpq_mat_print(G);
#endif // DEBUG_LEVEL_FULL

    ii = k;
    // infty
    m = ULONG_MAX;

    for (i = k; i < n; i++) {
      if (!fmpq_is_zero(fmpq_mat_entry(G, i, i))) {
	val = fmpq_valuation(fmpq_mat_entry(G, i, i), p);
	if (val < m) {
	  m = val;
	  ii = i;
	}
      }
    }

    i1 = ii;
    i2 = ii;
    
    for (i = k; i < n; i++)
      for (j = i+1; j < n; j++) {
	if (!fmpq_is_zero(fmpq_mat_entry(G, i, j))) {
	  val = fmpq_valuation(fmpq_mat_entry(G, i, j), p);
	  if (val < m) {
	    m = val;
	    i1 = i;
	    i2 = j;
	  }
	}
      }

#ifdef DEBUG_LEVEL_FULL
    printf("i_pair = (%lu,%lu)\nm = %lu\n", i1, i2, m);
#endif // DEBUG_LEVEL_FULL

     if (m != old_val) {
       blocks[jordan->num_blocks] = k;
       old_val = m;
       jordan->exponents[jordan->num_blocks++] = m;
     }

     
#ifdef DEBUG_LEVEL_FULL
     printf("blocks = \n");
     printf("jordan->exponents = \n");
#endif // DEBUG_LEVEL_FULL

     if ((even) && (i1 != i2)) {
       fmpq_mat_swap_rows(S, NULL, i1, k);
       fmpq_mat_swap_rows(S, NULL, i2, k+1);
	 
       // T12 = S[k]*F*S[k+1]^t
       inner_product(fmpq_mat_entry(T,0,1), B, S, k, k+1);

       // multiply S[k] by p^val(T12,p)/T12
       // Check whether we have to change to rational here
       for (i = 0; i < n; i++) {
	 fmpq_mul_2exp(fmpq_mat_entry(S, k, i), fmpq_mat_entry(S, k, i),
		       fmpq_valuation(fmpq_mat_entry(T,0,1), p));
	 fmpq_div(fmpq_mat_entry(S, k, i), fmpq_mat_entry(S, k, i), fmpq_mat_entry(T,0,1));
       }
       
       inner_product(fmpq_mat_entry(T,0,0), B, S, k, k);
       inner_product(fmpq_mat_entry(T,1,1), B, S, k+1, k+1);
       inner_product(fmpq_mat_entry(T,0,1), B, S, k, k+1);
       fmpq_set(fmpq_mat_entry(T,1,0), fmpq_mat_entry(T,0,1));
       fmpq_mat_det(d, T);

       for (l = k+2; l < n; l++) {
	 inner_product(tl, B, S, k+1, l);
	 fmpq_mul(tl, tl, fmpq_mat_entry(T,0,1));
	 inner_product(tl2, B, S, k, l);
	 fmpq_mul(tl2, tl2, fmpq_mat_entry(T,1,1));
	 fmpq_sub(tl,tl,tl2);
	 inner_product(ul, B, S, k, l);
	 fmpq_mul(ul, ul, fmpq_mat_entry(T,0,1));
	 inner_product(ul2, B, S, k+1, l);
	 fmpq_mul(ul2, ul2, fmpq_mat_entry(T,0,0));
	 fmpq_sub(ul,ul,ul2);
	 fmpq_div(tl, tl, d);
	 fmpq_div(ul, ul, d);
	 for (i = 0; i < n; i++) {
	   fmpq_addmul(fmpq_mat_entry(S,l,i), tl, fmpq_mat_entry(S,k,i));
	   fmpq_addmul(fmpq_mat_entry(S,l,i), ul, fmpq_mat_entry(S,k+1,i));
	 }
       }
       k += 2;
     } else {
       if (i1 == i2) {
	 
#ifdef DEBUG_LEVEL_FULL
	 printf("swapping rows\n");
#endif // DEBUG_LEVEL_FULL
	 fmpq_mat_swap_rows(S, NULL, i1, k);
	   
#ifdef DEBUG_LEVEL_FULL
	 printf("S = \n");
	 fmpq_mat_print(S);
#endif // DEBUG_LEVEL_FULL
       } else {
	 // std::cerr << "adding rows" << std::endl;
	 for (j = 0; j < n; j++)
	   fmpq_add(fmpq_mat_entry(S,i1,j), fmpq_mat_entry(S,i1,j), fmpq_mat_entry(S,i2,j));
	     
#ifdef DEBUG_LEVEL_FULL
	 printf("S = \n");
	 fmpq_mat_print(S);
	 printf("swapping rows\n");
#endif // DEBUG_LEVEL_FULL
	 fmpq_mat_swap_rows(S, NULL, i1, k);
#ifdef DEBUG_LEVEL_FULL
	 printf("S = \n");
	 fmpq_mat_print(S);
#endif // DEBUG_LEVEL_FULL
       }
       inner_product(nrm, B, S, k, k);

#ifdef DEBUG_LEVEL_FULL
       printf("nrm = ");
       fmpq_print(nrm);
       printf("\n");
#endif // DEBUG_LEVEL_FULL
       for (i = 0; i < n; i++)
	 inner_product(fmpq_mat_entry(X, 0, i), B, S, k, i);
	 
#ifdef DEBUG_LEVEL_FULL
       printf("X = ");
       fmpq_mat_print(X);
#endif // DEBUG_LEVEL_FULL
       
       for (l = k+1; l < n; l++)
	 for (i = 0; i < n; i++) {
	   fmpq_div(tl,fmpq_mat_entry(X,0,l), nrm);
	   fmpq_mul(tl, tl, fmpq_mat_entry(S,k,i));
	   fmpq_sub(fmpq_mat_entry(S,l,i), fmpq_mat_entry(S,l,i), tl);
	 }
	 
#ifdef DEBUG_LEVEL_FULL
       printf("S = ");
       fmpq_mat_print(S);
#endif // DEBUG_LEVEL_FULL
       k += 1;
     }
  }

  blocks[jordan->num_blocks] = n;

#ifdef DEBUG_LEVEL_FULL
  printf("blocks = ");
  for (i = 0; i <= jordan->num_blocks; i++)
    printf("%lu ", blocks[i]);
  printf("\n");
#endif // DEBUG_LEVEL_FULL

  for (i = 0; i < jordan->num_blocks; i++) {
    nrows = blocks[i+1]-blocks[i];
    fmpq_mat_init(jordan->matrices[i], nrows, n);
    for (row = 0; row < nrows; row++)
      for (col = 0; col < n; col++)
	fmpq_set(fmpq_mat_entry(jordan->matrices[i], row, col),
		 fmpq_mat_entry(S, blocks[i]+row, col));
  }

  for (i = 0; i < jordan->num_blocks; i++) {
    fmpq_mat_init(m_t, fmpq_mat_ncols(jordan->matrices[i]), fmpq_mat_nrows(jordan->matrices[i]));
    fmpq_mat_transpose(m_t, jordan->matrices[i]);
    fmpq_mat_init(mF, fmpq_mat_nrows(jordan->matrices[i]), fmpq_mat_ncols(F));
    fmpq_mat_mul(mF, jordan->matrices[i], F);
    fmpq_mat_init(jordan->grams[i], fmpq_mat_nrows(mF), fmpq_mat_nrows(mF));
    fmpq_mat_mul(jordan->grams[i], mF, m_t);
    
#ifdef DEBUG_LEVEL_FULL
    printf("m = ");
    fmpq_mat_print(jordan->matrices[i]);
    printf("F = ");
    fmpq_mat_print(F);
    printf("m^t = ");
    fmpq_mat_print(m_t);
    fmpq_init(tmp_rat);
    fmpq_mul(tmp_rat, fmpq_mat_entry(jordan->matrices[i], 0, 0), fmpq_mat_entry(F, 0, 0));
    printf("tmp_rat = ");
    fmpq_print(tmp_rat);
    printf("\n");
    fmpq_clear(tmp_rat);
    printf("m*F = ");
    fmpq_mat_print(mF);
    printf("m*F*m^t = ");
    fmpq_mat_print(jordan->grams[i]);
#endif // DEBUG_LEVEL_FULL
    
    fmpq_mat_clear(m_t);
    fmpq_mat_clear(mF);
  }
#ifdef DEBUG_LEVEL_FULL
  printf("jordan->matrices = \n");
 
  for (i = 0; i < jordan->num_blocks; i++)
    fmpq_mat_print(jordan->matrices[i]);

  printf("jordan->grams = \n");
  
  for ( i = 0; i < jordan->num_blocks; i++)
    fmpq_mat_print(jordan->grams[i]);
  
  printf("jordan->exponents = \n");
  for (i = 0; i < jordan->num_blocks; i++)
    printf("%lu ", jordan->exponents[i]);
  printf("\n");

#endif // DEBUG_LEVEL_FULL
  
  fmpq_clear(nrm);
  fmpq_clear(tl);
  fmpq_clear(tl2);
  fmpq_clear(ul);
  fmpq_clear(ul2);
  fmpq_clear(d);
  fmpq_mat_clear(T);
  fmpq_mat_clear(X);
  fmpz_mat_clear(B);
  fmpq_mat_clear(F);
  fmpq_mat_clear(G);
  fmpq_mat_clear(S);
  free(blocks);
 
  return;
}
