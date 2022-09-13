#include "arith.h"
#include "orthogonalize.h"
#include "qf_inv.h"

void qf_inv_init(qf_inv_t hasse, ulong n)
{
  slong i = 0;
  slong NUM_INIT = 16;
  
  hasse->num_stored = 0;
  hasse->I = 0;
  hasse->num_alloc = NUM_INIT;

  hasse->primes = (fmpz_t*)malloc(NUM_INIT * sizeof(fmpz_t));
  hasse->symbols = (slong*)malloc(NUM_INIT * sizeof(slong));

  hasse->n = n;

  for (i = 0; i < NUM_INIT; i++)
    fmpz_init(hasse->primes[i]);
  
  return;
}

void qf_inv_clear(qf_inv_t hasse)
{
  slong i;
  
  for (i = 0; i < hasse->num_alloc; i++)
    fmpz_clear(hasse->primes[i]);
  free(hasse->primes);
  free(hasse->symbols);
}

void qf_inv_add(qf_inv_t hasse, const fmpz_t p, slong s)
{
  ulong i;
  
  // right now implementing a silly search
  for (i = 0; i < hasse->num_stored; i++)
    if (fmpz_equal(hasse->primes[i], p)) return;
  
  if (hasse->num_stored == hasse->num_alloc) {
    hasse->num_alloc <<= 1;
    hasse->primes = (fmpz_t*)realloc(hasse->primes, hasse->num_alloc * sizeof(fmpz_t));
    hasse->symbols = (slong*)realloc(hasse->symbols, hasse->num_alloc * sizeof(slong));
  }
  for (i = hasse->num_stored; i < hasse->num_alloc; i++)
    fmpz_init(hasse->primes[i]);
  
  fmpz_set(hasse->primes[hasse->num_stored], p);
  hasse->symbols[hasse->num_stored] = s;
  hasse->num_stored++;

  return;
}

slong hasse(fmpz_mat_t D, fmpz_t p)
{
  slong s, i, n;
  fmpz_t prod;

  fmpz_init_set_ui(prod, 1);
  n = fmpz_mat_ncols(D);
  s = 1;

  for (i = 0; i < n; i++)
    fmpz_mul(prod, prod, fmpz_mat_entry(D,0,i));
  
  for (i = 0; i < n-1; i++) {
    fmpz_fdiv_q(prod, prod, fmpz_mat_entry(D,0,i));
    s *= hilbert_symbol(fmpz_mat_entry(D,0,i), prod, p);
  }
  
  return s;
}

void invariants(fmpz_t prod, qf_inv_t F, const fmpz_mat_t Q)
{
  fmpz_factor_t facs;
  fmpz_mat_t D;
  slong n, i, j;
  fmpz_t d, two;

  n = fmpz_mat_nrows(Q);
  fmpz_init(d);
  fmpz_mat_init(D, 1, n);
  fmpz_factor_init(facs);
  fmpz_init_set_ui(two, 2);

  orthogonalize_gram(D, Q);

  qf_inv_init(F, n);
  qf_inv_add(F, two, hasse(D,two));
  
  for (i = 0; i < n; i++) {
    fmpz_set(d, fmpz_mat_entry(D,0,i));
    if (d < 0) (F->I)++;
    fmpz_factor(facs, d);
    for (j = 0; j < facs->num; j++)
      if (facs->exp[j] % 2 == 1)
	qf_inv_add(F, &(facs->p[j]), hasse(D, &(facs->p[j])));
  }

  fmpz_set_ui(prod, 1);
  
  for (i = 0; i < n; i++)
    fmpz_mul(prod, prod, fmpz_mat_entry(D,0,i));

  fmpz_clear(prod);
  fmpz_clear(two);
  fmpz_factor_clear(facs);
  fmpz_mat_clear(D);
  fmpz_clear(d);
  
  return;
}

int witt_to_hasse(fmpz_mat_t hasse,
		  const fmpz_t det,
		  const qf_inv_t finite)
{
  slong c_table[8] = {2, 1, 1, -2, -2, -1, -1, 2};
  slong c_mask = c_table[finite->n % 8];
  fmpz_t c, minus_one;
  slong i;
  int hilbert, num_hasse;

  fmpz_init_set_si(c, c_mask % 2);
  fmpz_addmul_si(c, det, c_mask / 2);

  fmpz_init_set_si(minus_one, -1);

  num_hasse = 0;
  for (i = 0; i < finite->num_stored; i++) {
    hilbert = hilbert_symbol(minus_one, c, finite->primes[i]);
    if (hilbert != finite->symbols[i]) {
      fmpz_set(fmpz_mat_entry(hasse, 0, num_hasse), finite->primes[i]);
      num_hasse++;
    }
  }

  fmpz_clear(minus_one);
  fmpz_clear(c);
  return num_hasse;
}
