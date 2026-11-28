/*****************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : hash.c
 *
 * Description: Hash table for quickly testing
 *              isometry between lattices.
 *
 *****************************************************************
 */

// System dependencies

#include <assert.h>

// Required packages dependencies

#include <carat/symm.h>
#include <flint/ulong_extras.h>

// Self dependencies

#include <time.h>

#include "aut_grp.h"
#include "hash.h"
#include "isometry.h"
#include "matrix_tools.h"
#include "nbr_data.h"
#include "neighbor.h"

// an initial number of terms to use from the theta series
// TODO : We should be able to determine thes in a smarter way

#define INITIAL_THETA_PREC  25

/* hash_vec is a hash function on a vector of integers */
// hash_vec is defined in the end, as different choices of hash functions might lead to extremely
// different performance times
hash_t hash_vec(const int* vec, unsigned int n);

/* a helper function to add an entry to the table */
int _add(hash_table_t table, const square_matrix_t key, hash_t val, int do_push_back);

/* hash a form using its theta series */
// This should be the declaration, but for some reason short_vectors doesn't restrict the pointer to be constant
// Since we don't want to copy the matrix, we leave it at that.
hash_t hash_form(const square_matrix_t q, W32 theta_prec)
{
  matrix_TYP *Q, *SV;
  square_matrix_t red;
  hash_t x;
  int norm;
  int* num_short;

  assert(square_matrix_is_positive_definite(q));
  // The numbers of vectors of each norm do not depend on the basis, and the search for them is very
  // slow in a basis which is far from reduced, so they are counted in a reduced basis.
  square_matrix_set(red, q);
  if (square_matrix_is_far_from_reduced(q) && !greedy_reduce_try(red))
    square_matrix_reduce(red);
  Q = matrix_TYP_init_set_square_matrix(red);
  num_short = (int*)malloc(theta_prec * sizeof(int));
  for (norm = 2; norm <= 2*theta_prec; norm += 2) {
    SV = short_vectors(Q, norm, norm, 0, 1, &(num_short[norm/2 - 1]));
    free_mat(SV);
  }

  x = hash_vec(num_short, theta_prec);

  free_mat(Q);
  free(num_short);
  
  return x;
}

/* initialize a hash table of given size */
void hash_table_init(hash_table_t table, hash_t hash_size)
{
  hash_t i, log_hash_size, n;

  log_hash_size = 0;
  n = hash_size;
  while (n >>= 1) {
    log_hash_size++;
  }

  // It is easier to work with powers of 2 for memory allocation
  table->capacity = 1LL << log_hash_size;
  table->mask = (1LL << (log_hash_size + 1)) - 1;
  table->keys = (square_matrix_t*) malloc(table->capacity * sizeof(square_matrix_t));
  table->vals = (hash_t*) malloc(table->capacity * sizeof(hash_t));
  table->key_ptr = (hash_t*) malloc(2*table->capacity * sizeof(int));
  table->counts = (W32*) malloc(2*table->capacity * sizeof(int));
  table->probs = (fmpq_t*) malloc(table->capacity * sizeof(fmpq_t));
  for (i = 0; i < table->capacity; i++)
    fmpq_init(table->probs[i]);
  table->num_stored = 0;
  table->theta_prec = INITIAL_THETA_PREC;
  table->red_on_isom = TRUE;
  table->isom_cache = NULL;

  for (i = 0; i < 2*table->capacity; i++) {
    table->key_ptr[i] = -1;
    table->counts[i] = 0;
  }
  
  return;
}

/* Compute the time spent for an isometry testing, red_cost returns the amount of type spent on reduction */
// At the moment, done by sampling neighbors, and testing isometries between them.
// The neighbors are those of at most ISOM_COST_SAMPLE lattices of the table, spread over it: the tests
// cost about a millisecond each, and NUM_ISOMS tests at every lattice took most of the time of the
// genus at large discriminants (20 s of 32 s at the discriminant 997).
double get_isom_cost(const hash_table_t table, double* red_cost)
{
  double isom_cost;
  clock_t cputime;
  hash_t idx, offset, i, sample, num_sample;
  
#ifdef NBR_DATA
  nbr_data_t nbr_man;
  fmpz_mat_t nbr_fmpz, nbr_isom;
#else
  neighbor_manager_t nbr_man;
  bool has_ended;
  vector_t x;
  Z64 p;
#endif // NBR_DATA
  
  isometry_t s;
  square_matrix_t nbr;

#define NUM_ISOMS 10
#define ISOM_COST_SAMPLE 20
  
#ifdef NBR_DATA
  fmpz_mat_init(nbr_isom, QF_RANK, QF_RANK);
  fmpz_mat_init(nbr_fmpz, QF_RANK, QF_RANK);
#endif // NBR_DATA

  isom_cost = 0;
  if (red_cost != NULL)
    *red_cost = 0;
  idx = 0;
  num_sample = (table->num_stored < ISOM_COST_SAMPLE) ? table->num_stored : ISOM_COST_SAMPLE;
  for (sample = 0; sample < num_sample; sample++) {
    offset = (sample * table->num_stored) / num_sample;
#ifdef NBR_DATA
    nbr_data_init(nbr_man, table->keys[offset], 3, 1);
#else
    // The neighbors are at 3, or at the first odd prime at which the lattice
    // has an isotropic vector that nbr_process_init can start from.
    p = 3;
    while (!get_isotropic_vector(x, table->keys[offset], p))
      p = n_nextprime(p, true);
    isometry_init(s);
    do {
      nbr_process_init(nbr_man, table->keys[offset], p, idx, s);
      idx++;
      has_ended = nbr_process_has_ended(nbr_man);
      if (has_ended)
	nbr_process_clear(nbr_man);
    } while (has_ended);
    isometry_clear(s);
#endif // NBR_DATA
    for (i = 0; i <  NUM_ISOMS; i++) {
#ifdef NBR_DATA
      nbr_data_build_neighbor(nbr_fmpz, nbr_isom, nbr_man);
      square_matrix_set_fmpz_mat(nbr, nbr_fmpz);
#else
      nbr_process_build_nb(nbr, nbr_man);
#endif // NBR_DATA
      isometry_init(s);
      if (red_cost != NULL) {
	cputime = clock();
	greedy_checked(nbr, s);
	(*red_cost) += (clock() - cputime);
      }
      cputime = clock();
      is_isometric(s, table->keys[i % table->num_stored], nbr);
      isom_cost += (clock() - cputime);
      isometry_clear(s);
#ifdef NBR_DATA
      nbr_data_get_next_neighbor(nbr_man);
      if (nbr_data_has_ended(nbr_man)) {
	idx++;
	nbr_data_clear(nbr_man);
	nbr_data_init(nbr_man, table->keys[offset], 3, 1);
      }
#else
      nbr_process_advance(nbr_man);
      isometry_init(s);
      while (nbr_process_has_ended(nbr_man)) {
	nbr_process_clear(nbr_man);
	nbr_process_init(nbr_man, table->keys[offset], p, idx, s);
	idx++;
      }
      isometry_clear(s);
#endif // NBR_DATA
    }
#ifdef NBR_DATA
    nbr_data_clear(nbr_man);
#else
    nbr_process_clear(nbr_man);
#endif // NBR_DATA
  }
  isom_cost /= (NUM_ISOMS * num_sample);
  if (red_cost != NULL)
    (*red_cost) /= (NUM_ISOMS * num_sample);

  return isom_cost;
}

/* get_total_cost returns the average time spent when looking up an isometry class.
   This is done by computing the time spent for computing the theta series,
   using the collision counts and the sizes of automorphism groups (probs) to estimate the average number
   of isometry tests on every call, and computing the weighted average based on the given timing
   of a single isometry test. We also input the time of checking an isometry with reduction.
   Given the timings, we decide whether to reduce when performing isometry testing (returned in red_on_isom),
   and output the total cost accordingly. */
double get_total_cost(const hash_table_t table, W32 theta_prec, double isom_cost, double red_isom_cost, fmpq_t* wt_cnts,
		      W32* counts, hash_t* vals, bool* red_on_isom, const fmpq_t* probs)
{
  hash_t i, offset, idx;
  double theta_cost, total_cost;
  clock_t cputime;
  
  for (i = 0; i < 2*table->capacity; i++) {
    fmpq_zero(wt_cnts[i]);
    counts[i] = 0;
  }
  
  theta_cost = 0;
  for (offset = 0; offset < table->num_stored; offset++) {
    cputime = clock();
    vals[offset] = hash_form(table->keys[offset], theta_prec);
    theta_cost += (clock() - cputime);
    idx = vals[offset] & table->mask;
    fmpq_add(wt_cnts[idx], wt_cnts[idx], probs[offset]);
    counts[vals[offset] & table->mask]++;
  }
  theta_cost /= table->num_stored;
    
  total_cost = 0;
  for (i = 0; i < 2 * table->capacity; i++) {
    total_cost += counts[i] * fmpq_get_d(wt_cnts[i]);
  }
  total_cost -= 1;
#ifdef DEBUG_LEVEL_FULL
  printf("Expecting average number of %f calls to is_isometric.\n", total_cost);
#endif // DEBUG_LEVEL_FULL

  *red_on_isom = true;

  // !! TODO - for large p this seems necessary to prevent overflow
  // check also if changing to 64-bit helps
  /*
  for (offset = 0; offset < table->num_stored; offset++)
    *red_on_isom &= (counts[vals[offset] & table->mask] == 1);

  *red_on_isom = !(*red_on_isom);
  */

  if (*red_on_isom)
    total_cost *= red_isom_cost;
  else
    total_cost *= isom_cost;
  
  total_cost += theta_cost;
    
  return total_cost;
}

/*  recalibrate the hash table to use a different number of terms from the theta series, to optimize performance. */
// This can be done better, by first computing the theta series and
// only then going over them. Since this is not a critical path,
// we postpone doing that.

// table is not const, since we're using it and clearing it in the end
void hash_table_recalibrate(hash_table_t new_table, const hash_table_t table)
{
  hash_t offset, i;
  W32 theta_prec;
  double greedy_cost, red_isom_cost, isom_cost;
  double total_cost;
  double old_cost;
  fmpq_t* wt_cnts;
  fmpq_t mass;

  fmpq_init(mass);
  fmpq_zero(mass);

  hash_table_init(new_table, table->capacity);
  
  wt_cnts = (fmpq_t*) malloc(2*table->capacity*sizeof(fmpq_t));
  for (i = 0; i < 2 * table->capacity; i++) {
    fmpq_init(wt_cnts[i]);
  }

  for (offset = 0; offset < table->num_stored; offset++) {
    fmpq_add(mass, mass, table->probs[offset]);
  }

  for (offset = 0; offset < table->num_stored; offset++) {
    fmpq_div(new_table->probs[offset], table->probs[offset], mass);
  }
  
  // the cost of checking isometries with random neighbors, reducing them first. get_total_cost always
  // reduces them (red_on_isom), so the cost without reduction is not used and is not measured.
  red_isom_cost = get_isom_cost(table, &greedy_cost);
  isom_cost = red_isom_cost + greedy_cost;

#ifdef DEBUG_LEVEL_FULL
  printf("Expected cost of reduced isometries is %f\n", red_isom_cost + greedy_cost);
#endif // DEBUG_LEVEL_FULL
  
  theta_prec = 0;

  total_cost = (table->num_stored - 1) * (red_isom_cost + greedy_cost);
  old_cost = total_cost + 1;

  // finding the minimum by brute force
  while (total_cost < old_cost) {
    theta_prec++;
    old_cost = total_cost;
    total_cost = get_total_cost(table, theta_prec, isom_cost, red_isom_cost + greedy_cost, wt_cnts,
				new_table->counts, new_table->vals, &(new_table->red_on_isom), new_table->probs);
  }
  // we are now passed the peak, have to go back one
  theta_prec--;

#ifdef DEBUG_LEVEL_FULL
  printf("Recalibrated with theta_prec = %u and red_on_isom = %d \n", theta_prec, new_table->red_on_isom);
#endif // DEBUG_LEVEL_FULL
  
  new_table->theta_prec = theta_prec;

  // revert the used counts
  for (i = 0; i < 2*new_table->capacity; i++) {
    new_table->counts[i] = 0;
  }
  
  // adding all keys
  for (offset = 0; offset < table->num_stored; offset++)
    hash_table_add(new_table, table->keys[offset]);

  free(wt_cnts);
  return;
}

/* insert a key to the hash table at a specific index */
// key is non-const to prevent copying the matrices
int hash_table_insert(hash_table_t table, const square_matrix_t key, hash_t val,
		      int index, int do_push_back)
{
  int offset;
  aut_grp_t aut_grp;
  
  if (table->num_stored == table->capacity) {
    hash_table_expand(table);
    return _add(table, key, val, do_push_back);
  }

  offset = table->num_stored;
  ++table->num_stored;

  if (do_push_back) {
    square_matrix_set(table->keys[offset],key);
    table->vals[offset] = val;
    aut_grp_init_square_matrix(aut_grp, key);
    fmpq_set_si(table->probs[offset], 1, aut_grp->order);
    aut_grp_clear(aut_grp);
  }

  table->key_ptr[index] = offset;
  table->counts[(val & table->mask)]++;
  
  return 1; 
}

/* a helper function to add an entry to the hash table */
int _add(hash_table_t table, const square_matrix_t key, hash_t val, int do_push_back)
{
  int offset, i;
  hash_t index;
  isometry_t isom;

  index = val & table->mask;
  i = 0;
  while (i <= table->counts[val & table->mask]) {
    offset = table->key_ptr[index];
    if (offset == -1) 
      return hash_table_insert(table, key, val, index, do_push_back);
    // we first check equal hashes before testing an isometry
    if ((table->vals[offset] & table->mask) == (val & table->mask)) {
      i++;
      if (table->vals[offset] == val) {
	if (is_isometric(isom, table->keys[offset], key))
	  return 0;
      }
    }
    index = (index + 1) & table->mask;
  }

  offset = table->key_ptr[index];
  while (offset != -1) {
    index = (index + 1) & table->mask;
    offset = table->key_ptr[index];
  }
  return hash_table_insert(table, key, val, index, do_push_back);
  
}

/* add a lattice to the table */
int hash_table_add(hash_table_t table, const square_matrix_t key)
{
  return _add(table, key, hash_form(key, table->theta_prec), 1);
}

/* a modified version of exists to be able to run over all with the same
hash value */
bool hash_table_get_key(square_matrix_t new_key, const hash_table_t table, const square_matrix_t key, int* index)
{
  int offset;
  hash_t value;

  value = hash_form(key, table->theta_prec);
  
  if (*index == -1) 
    *index = value & table->mask;
  else
    *index = (*index + 1) & table->mask;
  
  offset = table->key_ptr[*index];
  if (offset == -1)
    return false;

  while (hash_form(table->keys[offset], table->theta_prec) != value) {
    *index = (*index + 1) & table->mask;
    offset = table->key_ptr[*index];
    if (offset == -1)
      return false;
  }

  square_matrix_set(new_key, table->keys[offset]);
  
  return true;
}

/* Checks whether there exists a lattice with the same hash in the table.
   If check_isom is set, checks whether there is an isometric one. */
int hash_table_exists(const hash_table_t table, const square_matrix_t key, int check_isom)
{
  int offset, i;
  hash_t index, value;
  isometry_t isom;

  value = hash_form(key, table->theta_prec);
  
  index = value & table->mask;
  i = 0;
  while (i < table->counts[value & table->mask]) {
    offset = table->key_ptr[index];
    if (offset == -1)
      return 0;
    if ((table->vals[offset] & table->mask) == (value & table->mask)) {
      if ((table->counts[value & table->mask] == i + 1) && (!check_isom))
	return 1;
      i++;
      if (table->vals[offset] == value) {
	if (is_isometric(isom, table->keys[offset], key))
	   return 1;
      }
    }

    index = (index + 1) & table->mask;
  }

  return 0;
}

/* Returns the index in the hash table of a lattice with the same hash value,
   or an isometric lattice if check_isom is set. */

int hash_table_indexof(const hash_table_t table, const square_matrix_t key, int check_isom, double* theta_time, double* isom_time, slong* num_isom)
{
  int offset, i;
  hash_t index, value;
  clock_t cputime;
  isometry_t s;
  square_matrix_t key_copy;
  bool is_reduced;

  cputime = clock();
  // a key which is far from reduced is reduced once, for its theta series (which does not depend on
  // the basis) and for the isometry tests; another only when it is tested for isometry
  square_matrix_set(key_copy, key);
  is_reduced = false;
  if (table->red_on_isom && square_matrix_is_far_from_reduced(key)) {
    isometry_init(s);
    greedy_checked(key_copy, s);
    isometry_clear(s);
    is_reduced = true;
  }
  value = hash_form(key_copy, table->theta_prec);
  (*theta_time) += clock() - cputime;
  index = value & table->mask;
  i = 0;
  while (i < table->counts[value & table->mask]) {
    offset = table->key_ptr[index];
    assert(offset != -1);
    if ((table->vals[offset] & table->mask) == (value & table->mask)) {
      if ((table->counts[value & table->mask] == i + 1) && (!check_isom))
	return offset;
      i++;
      if (table->vals[offset] == value) {
	(*num_isom)++;
	cputime = clock();
	if (table->red_on_isom && !is_reduced) {
	  isometry_init(s);
	  greedy_checked(key_copy, s);
	  isometry_clear(s);
	  is_reduced = true;
	}
	if (is_isometric(s, table->keys[offset], key_copy)) {
	  (*isom_time) += clock() - cputime;
	  return offset;
	}
	(*isom_time) += clock() - cputime;
      }
    }

    index = (index + 1) & table->mask;
  }

  assert(false);

  return -1;
}

// ---- The cache of hash_table_index_and_isom ----
//
// The neighbors of the lattices of a genus reduce to few forms: at D = 148 and the weight (2, 1), 4.4
// million neighbors with -prec=100 reduce to 397 forms. The class of a form is unique, and the isometry
// which is_isometric finds from it to a key of the table depends only on the two forms, so the cache
// keeps, for each reduced form found, the offset of its class and that isometry: the lookup then gives
// what the search gives. It holds at most ISOM_CACHE_MAX forms (about 35 MB), and is cleared with the table.
// is_isometric gives CARAT no automorphisms, so its search does not use rand() and depends on nothing but
// the two forms; if it is ever given generators of a group, the cache has to be checked again.

#define ISOM_CACHE_MAX (1 << 16)
#define ISOM_CACHE_LEN (QF_RANK * (QF_RANK + 1) / 2)

typedef struct {
  slong cap, n;
  Z64* forms;              // ISOM_CACHE_LEN entries for each slot, the upper triangle of the form
  int* offsets;            // the offset of the class, or -1 for an empty slot
  isometry_struct* isoms;  // the isometry from the form to the key of the class
} isom_cache_struct;

static void isom_cache_form(Z64* v, const square_matrix_t q)
{
  int r, c, k = 0;
  for (r = 0; r < QF_RANK; r++)
    for (c = r; c < QF_RANK; c++)
      v[k++] = q[r][c];
}

static slong isom_cache_slot(const isom_cache_struct* C, const Z64* v)
{
  ulong h = UWORD(1469598103934665603);
  slong slot;
  int k;
  for (k = 0; k < ISOM_CACHE_LEN; k++)
    h = (h ^ (ulong)v[k]) * UWORD(1099511628211);
  for (slot = (slong)(h & (ulong)(C->cap - 1)); C->offsets[slot] >= 0; slot = (slot + 1) & (C->cap - 1)) {
    for (k = 0; (k < ISOM_CACHE_LEN) && (C->forms[ISOM_CACHE_LEN * slot + k] == v[k]); k++);
    if (k == ISOM_CACHE_LEN)
      break;
  }
  return slot;
}

static void isom_cache_alloc(isom_cache_struct* C, slong cap)
{
  slong i;
  C->cap = cap;
  C->n = 0;
  C->forms = (Z64*)malloc(cap * ISOM_CACHE_LEN * sizeof(Z64));
  C->offsets = (int*)malloc(cap * sizeof(int));
  C->isoms = (isometry_struct*)malloc(cap * sizeof(isometry_struct));
  for (i = 0; i < cap; i++)
    C->offsets[i] = -1;
}

static void isom_cache_free(isom_cache_struct* C)
{
  free(C->forms);
  free(C->offsets);
  free(C->isoms);
}

static void isom_cache_clear(hash_table_t table)
{
  if (table->isom_cache != NULL) {
    isom_cache_free((isom_cache_struct*)table->isom_cache);
    free(table->isom_cache);
    table->isom_cache = NULL;
  }
}

// the offset of the class of the reduced form q and the isometry from q to its key, or -1
static int isom_cache_find(isometry_t isom, const hash_table_t table, const square_matrix_t q)
{
  const isom_cache_struct* C = (const isom_cache_struct*)table->isom_cache;
  Z64 v[ISOM_CACHE_LEN];
  slong slot;
  if (C == NULL)
    return -1;
  isom_cache_form(v, q);
  slot = isom_cache_slot(C, v);
  if (C->offsets[slot] < 0)
    return -1;
  isom[0] = C->isoms[slot];
  return C->offsets[slot];
}

static void isom_cache_add(const hash_table_t table, const square_matrix_t q, int offset, const isometry_t isom)
{
  // the table is const for its keys; the cache it points to is not
  hash_table* T = (hash_table*)table;
  isom_cache_struct* C = (isom_cache_struct*)T->isom_cache;
  Z64 v[ISOM_CACHE_LEN];
  slong slot, i;
  int k;
  if (C == NULL) {
    C = (isom_cache_struct*)malloc(sizeof(isom_cache_struct));
    isom_cache_alloc(C, 1024);
    T->isom_cache = C;
  }
  if (2 * (C->n + 1) > C->cap) {
    // at most half full: twice as many slots, up to twice ISOM_CACHE_MAX
    isom_cache_struct D;
    if (C->cap >= 2 * ISOM_CACHE_MAX)
      return;
    isom_cache_alloc(&D, 2 * C->cap);
    for (i = 0; i < C->cap; i++)
      if (C->offsets[i] >= 0) {
	slot = isom_cache_slot(&D, C->forms + ISOM_CACHE_LEN * i);
	for (k = 0; k < ISOM_CACHE_LEN; k++)
	  D.forms[ISOM_CACHE_LEN * slot + k] = C->forms[ISOM_CACHE_LEN * i + k];
	D.offsets[slot] = C->offsets[i];
	D.isoms[slot] = C->isoms[i];
	D.n++;
      }
    isom_cache_free(C);
    *C = D;
  }
  isom_cache_form(v, q);
  slot = isom_cache_slot(C, v);
  if (C->offsets[slot] >= 0)
    return;
  for (k = 0; k < ISOM_CACHE_LEN; k++)
    C->forms[ISOM_CACHE_LEN * slot + k] = v[k];
  C->offsets[slot] = offset;
  C->isoms[slot] = isom[0];
  C->n++;
}

/* Returns also an isometry between the lattices */

int hash_table_index_and_isom(const hash_table_t table, const square_matrix_t key, isometry_t isom,
			      double* theta_time, double* isom_time, slong* num_isom)
{
  int offset, i;
  hash_t index, value;
  clock_t cputime;
  isometry_t s;
  square_matrix_t key_copy;

  cputime = clock();
  // the key is reduced once, for its theta series (which does not depend on the basis) and for the
  // isometry tests; s takes key to key_copy
  square_matrix_set(key_copy, key);
  if (table->red_on_isom) {
    isometry_init(s);
    greedy_checked(key_copy, s);
#ifdef DEBUG_LEVEL_FULL
    assert(isometry_is_isom(s, key, key_copy));
#endif // DEBUG_LEVEL_FULL
  }
  offset = isom_cache_find(isom, table, key_copy);
  if (offset >= 0) {
    (*num_isom)++;
    if (table->red_on_isom) {
      isometry_muleq_left(isom, s);
      isometry_clear(s);
    }
    (*theta_time) += clock() - cputime;
    return offset;
  }
  value = hash_form(key_copy, table->theta_prec);
  (*theta_time) += clock() - cputime;
  index = value & table->mask;
  i = 0;
  while (i < table->counts[value & table->mask]) {
    offset = table->key_ptr[index];
    assert(offset != -1);
    if ((table->vals[offset] & table->mask) == (value & table->mask)) {
      i++;
      if (table->vals[offset] == value) {
	(*num_isom)++;
	cputime = clock();
	if (is_isometric(isom, key_copy, table->keys[offset])) {
#ifdef DEBUG_LEVEL_FULL
	  assert(isometry_is_isom(isom, key_copy, table->keys[offset]));
#endif // DEBUG_LEVEL_FULL
	  isom_cache_add(table, key_copy, offset, isom);
	  if (table->red_on_isom) {
	    isometry_muleq_left(isom, s);
#ifdef DEBUG_LEVEL_FULL
	    assert(isometry_is_isom(isom, key, table->keys[offset]));
#endif // DEBUg_LEVEL_FULL
	    isometry_clear(s);
	  }
	  (*isom_time) += clock() - cputime;
	  return offset;
	}
	(*isom_time) += clock() - cputime;
      }
    }

    index = (index + 1) & table->mask;
  }

  if (table->red_on_isom)
    isometry_clear(s);

  // the key is not in the table
  return -1;
}

/* Expand the size of the table in order to support more entries. */

void hash_table_expand(hash_table_t table)
{
  int i, stored;
  
  table->capacity <<= 1;
  table->mask = (table->capacity << 1)-1;
  // problem - this might invalidate existing references to the HashMap (e.g. mother)
  table->keys = (square_matrix_t*)realloc(table->keys, (table->capacity)*sizeof(square_matrix_t));
  table->vals = (hash_t*)realloc(table->vals, (table->capacity)*sizeof(hash_t));
  table->probs = (fmpq_t*)realloc(table->probs, (table->capacity)*sizeof(fmpq_t));

  // initializing the new ones
  for (i = (table->capacity)>>1; i < table->capacity; i++)
    fmpq_init(table->probs[i]);
  
  free(table->key_ptr);
  table->key_ptr = (hash_t*)malloc((table->capacity << 1)*sizeof(hash_t));
  free(table->counts);
  table->counts = (W32*)malloc((table->capacity << 1)*sizeof(W32));

  for (i = 0; i < 2*table->capacity; i++)
    table->key_ptr[i] = -1;

  for (i = 0; i < 2*table->capacity; i++)
    table->counts[i] = 0;

  stored = table->num_stored;
  table->num_stored = 0;

  for (i=0; i<stored; i++)
    _add(table, table->keys[i], table->vals[i], 0);

  return;
}

/* clear the memory allocated for the table */

void hash_table_clear(hash_table_t table)
{
  hash_t i;
  
  free(table->key_ptr);
  free(table->vals);
  free(table->keys);
  for (i = 0; i < table->capacity; i++)
    fmpq_clear(table->probs[i]);
  free(table->probs);
  free(table->counts);
  isom_cache_clear(table);

  return;
}

/*************************************************************************************
 * Here we list several different hash functions that can be used for the hash_table
 *
 *************************************************************************************/

hash_t RSHash(const int* vec, unsigned int n)
{
  hash_t b = 378551;
  hash_t a = 63689;
  hash_t hash = 0;
  unsigned int i = 0;
  for (i = 0; i < n; i++) {
    hash = hash * a + vec[i];
    a *= b;
  }
  return hash;
}

hash_t hash_vec(const int* vec, unsigned int n)
{
  return RSHash(vec, n);
}

/* hash_t hash_vec(int* vec, int n) */
/* { */
/*   hash_t fnv = FNV_OFFSET; */
/*   for (int i = 0; i < n; i++) */
/*     // too slow */
/*     // fnv = (fnv ^ vec[i]) * FNV_PRIME; */
/*     // This one follows a hashing scheme by D J Bernstein */
/*     fnv = ((fnv << 5) + fnv) + vec[i]; */
/*   return fnv; */
/* } */
