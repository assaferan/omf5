#include "carat/symm.h"

#include "hash.h"
#include "matrix_tools.h"

// TODO : We should be able to determine thes in a smarter way

#define INITIAL_THETA_PREC  25
#define MAX_THETA_PREC     100

// hash_vec is defined in the end, as different choices of hash functions might lead to extremely
// different performance times
hash_t hash_vec(const int* vec, unsigned int n);

hash_t hash_form(matrix_TYP* Q, W32 theta_prec)
{
  hash_t x;
  int norm;
  int* num_short;

  num_short = (int*)malloc(theta_prec * sizeof(int));
  for (norm = 2; norm <= 2*theta_prec; norm += 2) {
    short_vectors(Q, norm, norm, 0, 1, &(num_short[norm/2 - 1]));
  }

  // x = num_short[0] - num_short[1] + num_short[2] ;
  x = hash_vec(num_short, theta_prec);

  return x;
}

hash_table* create_hash(hash_t hash_size)
{
  hash_table* table;
  hash_t i, log_hash_size, n;

  table = (hash_table*) malloc(sizeof(hash_table));

  log_hash_size = 0;
  n = hash_size;
  while (n >>= 1) {
    log_hash_size++;
  }
  
  table->capacity = 1LL << log_hash_size;
  table->mask = (1LL << (log_hash_size + 1)) - 1;
  table->keys = (matrix_TYP**) malloc(table->capacity * sizeof(matrix_TYP*));
  table->vals = (hash_t*) malloc(table->capacity * sizeof(hash_t));
  table->key_ptr = (hash_t*) malloc(2*table->capacity * sizeof(int));
  table->counts = (W32*) malloc(2*table->capacity * sizeof(int));
  table->probs = (fmpq_t*) malloc(table->capacity * sizeof(fmpq_t));
  for (i = 0; i < table->capacity; i++)
    fmpq_init(table->probs[i]);
  table->num_stored = 0;
  table->theta_prec = INITIAL_THETA_PREC;

  for (i = 0; i < 2*table->capacity; i++) {
    table->key_ptr[i] = -1;
    table->counts[i] = 0;
  }
  
  return table;
}

// This can be done better, by first computing the theta series and
// only then going over them. Since this is not a critical path,
// we postpone doing that.

hash_table* recalibrate_hash(hash_table* table)
{
  hash_table* new_table;
  hash_t offset, i;
  W32 theta_prec;
  double total_cost, isom_cost, theta_cost, old_cost;
  clock_t cputime;
  fmpq_t* wt_cnts;
  hash_t idx;

  wt_cnts = (fmpq_t*) malloc(2*table->capacity*sizeof(fmpq_t));
  for (i = 0; i < 2 * table->capacity; i++)
    fmpq_init(wt_cnts[i]);
  
  #define NUM_ISOMS 100
  
  cputime = clock();
  for (i = 0; i <  NUM_ISOMS; i++)
    is_isometric(table->keys[0], table->keys[i % table->num_stored]);
  isom_cost = (clock() - cputime) / NUM_ISOMS;
  
  theta_prec = 0;

  total_cost = (table->num_stored - 1) * isom_cost;
  old_cost = total_cost + 1;

  // finding the minimum by brute force
  while (total_cost < old_cost) {
    theta_prec++;

    for (i = 0; i < 2*table->capacity; i++) {
      fmpq_zero(wt_cnts[i]);
      table->counts[i] = 0;
    }
    theta_cost = 0;
    for (offset = 0; offset < table->num_stored; offset++) {
      cputime = clock();
      table->vals[offset] = hash_form(table->keys[offset], theta_prec);
      theta_cost += (clock() - cputime);
      idx = table->vals[offset] & table->mask;
      fmpq_add(wt_cnts[idx], wt_cnts[idx], table->probs[offset]);
      table->counts[table->vals[offset] & table->mask]++;
    }
    theta_cost /= table->num_stored;
    
    old_cost = total_cost;
    total_cost = 0;
    for (i = 0; i < 2 * table->capacity; i++) {
      total_cost += table->counts[i] * fmpq_get_d(wt_cnts[i]);
    }
    total_cost -= 1;
    total_cost *= isom_cost;
    total_cost += theta_cost;
  }

#ifdef DEBUG
  printf("Recalibrated with theta_prec = %u\n", theta_prec);
#endif // DEBUG
  
  new_table = create_hash(table->num_stored);
  new_table->theta_prec = theta_prec;
  for (offset = 0; offset < table->num_stored; offset++)
    add(new_table, table->keys[offset]);

  free_hash(table);
  
  return new_table;
}

int insert(hash_table* table, matrix_TYP* key, hash_t val,
	   int index, int do_push_back)
{
  int offset;
  bravais_TYP* aut_grp;
  
  if (table->num_stored == table->capacity) {
    expand(table);
    return _add(table, key, val, do_push_back);
  }

  offset = table->num_stored;
  ++table->num_stored;

  if (do_push_back) {
    table->keys[offset] = key;
    table->vals[offset] = val;
    aut_grp = automorphism_group(key);
    fmpq_set_si(table->probs[offset], 1, aut_grp->order);
  }

  table->key_ptr[index] = offset;
  table->counts[(val & table->mask)]++;

  return 1; 
}

int _add(hash_table* table, matrix_TYP* key, hash_t val, int do_push_back)
{
  int offset, i;
  hash_t index;

  index = val & table->mask;
  i = 0;
  while (i <= table->counts[val & table->mask]) {
    offset = table->key_ptr[index];
    if (offset == -1) 
      return insert(table, key, val, index, do_push_back);
    // we first check equal hashes before testing an isometry
    if ((table->vals[offset] & table->mask) == (val & table->mask)) {
      i++;
      if (table->vals[offset] == val) {
	if (is_isometric(key,table->keys[offset]))
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
  return insert(table, key, val, index, do_push_back);
  
}

int add(hash_table* table, matrix_TYP* key)
{
  return _add(table, key, hash_form(key, table->theta_prec), 1);
}

/* a modified version of exists to be able to run over all with the same
hash value */
matrix_TYP* get_key(hash_table* table, matrix_TYP* key, int* index)
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
    return NULL;

  while (hash_form(table->keys[offset], table->theta_prec) != value) {
    *index = (*index + 1) & table->mask;
    offset = table->key_ptr[*index];
    if (offset == -1)
      return NULL;
  }
  
  return table->keys[offset];
}

int exists(hash_table* table, matrix_TYP* key, int check_isom)
{
  int offset, i;
  hash_t index, value;

  value = hash_form(key, table->theta_prec);
  
  index = value & table->mask;
  i = 0;
  while (i < table->counts[value & table->mask]) {
    offset = table->key_ptr[index];
    if (offset == -1)
      return 0;
    if ((table->vals[offset] & table->mask) == (value & table->mask)) {
      if ((table->counts[value & table->mask] == 1) && (!check_isom))
	return 1;
      i++;
      if (table->vals[offset] == value) {
	if (is_isometric(key,table->keys[offset]))
	   return 1;
      }
    }

    index = (index + 1) & table->mask;
  }

  return 0;
}

int indexof(hash_table* table, matrix_TYP* key, int check_isom)
{
  int offset, i;
  hash_t index, value;

  value = hash_form(key, table->theta_prec);
  index = value & table->mask;
  i = 0;
  while (i < table->counts[value & table->mask]) {
    offset = table->key_ptr[index];
    if (offset == -1) {
      printf("Error! Key not found!\n");
      return -1;
    }
    if ((table->vals[offset] & table->mask) == (value & table->mask)) {
      if ((table->counts[value & table->mask] == 1) && (!check_isom))
	return offset;
      i++;
      if (table->vals[offset] == value) {
	if (is_isometric(key,table->keys[offset]))
	   return offset;
      }
    }

    index = (index + 1) & table->mask;
  }

  printf("Error! Key not found!\n");
  exit(-1);

  return -1;
}

void expand(hash_table* table)
{
  int i, stored;
  
  table->capacity <<= 1;
  table->mask = (table->capacity << 1)-1;
  // problem - this might invalidate existing references to the HashMap (e.g. mother)
  table->keys = realloc(table->keys, (table->capacity)*sizeof(matrix_TYP*));
  table->vals = realloc(table->vals, (table->capacity)*sizeof(matrix_TYP*));

  free(table->key_ptr);
  table->key_ptr = (hash_t*)malloc((table->capacity << 1)*sizeof(int));
  free(table->counts);
  table->counts = (W32*)malloc((table->capacity << 1)*sizeof(int));

  for (i = 0; i < 2*table->capacity; i++)
    table->key_ptr[i] = -1;

  for (i = 0; i < 2*table->capacity; i++)
    table->counts[i] = 0;

  stored = table->num_stored;
  table->num_stored = 0;

  for (i=0; i<stored; i++)
    _add(table, table->keys[i], table->vals[i], 0);

}

void free_hash(hash_table* table)
{
  hash_t i;
  
  free(table->key_ptr);
  free(table->vals);
  free(table->keys);
  for (i = 0; i < table->capacity; i++)
    fmpq_clear(table->probs[i]);
  free(table->probs);
  free(table->counts);
  free(table);
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
