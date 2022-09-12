#include "carat/symm.h"

#include "hash.h"
#include "matrix_tools.h"

#define FNV_OFFSET 0x811c9dc5
#define FNV_PRIME 0x01000193

W64 hash_vec(int* vec, int n)
{
  W64 fnv = FNV_OFFSET;
  for (int i = 0; i < n; i++)
      fnv = (fnv ^ vec[i]) * FNV_PRIME;

  return fnv;
}

W64 hash_form(matrix_TYP* Q)
{
  W64 x;
  int norm;
  int num_short[3];
  
  for (norm = 2; norm <= 6; norm += 2) {
    short_vectors(Q, norm, norm, 0, 1, &(num_short[norm/2 - 1]));
  }

  x = num_short[0] - num_short[1] + num_short[2] ;
  // x = hash_vec(num_short, 3);

  return x;
}

hash_table* create_hash(int hash_size)
{
  hash_table* table;
  int i, log_hash_size, n;

  table = (hash_table*) malloc(sizeof(hash_table));

  log_hash_size = 0;
  n = hash_size;
  while (n >>= 1) {
    log_hash_size++;
  }
  
  table->capacity = 1LL << log_hash_size;
  table->mask = (1LL << (log_hash_size + 1)) - 1;
  table->keys = (matrix_TYP**) malloc(table->capacity * sizeof(matrix_TYP*));
  table->vals = (int*) malloc(table->capacity * sizeof(int));
  table->key_ptr = (int*) malloc(2*table->capacity * sizeof(int));
  table->counts = (int*) malloc(2*table->capacity * sizeof(int));
  table->num_stored = 0;

  for (i = 0; i < 2*table->capacity; i++) {
    table->key_ptr[i] = -1;
    table->counts[i] = 0;
  }
  
  return table;
}

int insert(hash_table* table, matrix_TYP* key, W64 val,
	   int index, int do_push_back)
{
  int offset;
  
  if (table->num_stored == table->capacity) {
    expand(table);
    return _add(table, key, val, do_push_back);
  }

  offset = table->num_stored;
  ++table->num_stored;

  if (do_push_back) {
    table->keys[offset] = key;
    table->vals[offset] = val;
  }

  table->key_ptr[index] = offset;
  table->counts[(val & table->mask)]++;

  return 1; 
}

int _add(hash_table* table, matrix_TYP* key, W64 val, int do_push_back)
{
  int offset, i;
  W64 index;

  index = val & table->mask;
  i = 0;
  while (i <= table->counts[val & table->mask]) {
    offset = table->key_ptr[index];
    if (offset == -1) 
      return insert(table, key, val, index, do_push_back);
    if (is_isometric(key,table->keys[offset]))
      return 0;

    index = (index + 1) & table->mask;
    i++;
  }

  while (offset != -1) {
    index = (index + 1) & table->mask;
    offset = table->key_ptr[index];
  }
  return insert(table, key, val, index, do_push_back);
  
}

int add(hash_table* table, matrix_TYP* key)
{
  return _add(table, key, hash_form(key), 1);
}

/* a modified version of exists to be able to run over all with the same
hash value */
matrix_TYP* get_key(hash_table* table, matrix_TYP* key, int* index)
{
  int offset;
  W64 value;

  value = hash_form(key);
  
  if (*index == -1) 
    *index = value & table->mask;
  else
    *index = (*index + 1) & table->mask;
  
  offset = table->key_ptr[*index];
  if (offset == -1)
    return NULL;

  while (hash_form(table->keys[offset]) != value) {
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
  W64 index, value;

  value = hash_form(key);
  
  index = value & table->mask;
  i = 0;
  while (i < table->counts[value & table->mask]) {
    offset = table->key_ptr[index];
    if (offset == -1)
      return 0;
    if ((table->counts[value & table->mask] == 1) && (!check_isom))
      return 1;
    if (is_isometric(key,table->keys[offset]))
      return 1;

    index = (index + 1) & table->mask;
    i++;
  }

  return 0;
}

int indexof(hash_table* table, matrix_TYP* key, int check_isom)
{
  int offset, i;
  W64 index, value;

  value = hash_form(key);
  index = value & table->mask;
  i = 0;
  while (i < table->counts[value & table->mask]) {
    offset = table->key_ptr[index];
    if (offset == -1) {
      printf("Error! Key not found!\n");
      return -1;
    }
    
    if ((table->counts[value & table->mask] == 1) && (!check_isom))
      return offset;
    
    if (is_isometric(key,table->keys[offset]))
      return offset;

    index = (index + 1) & table->mask;
    i++;
  }

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
  table->key_ptr = (int*)malloc((table->capacity << 1)*sizeof(int));
  free(table->counts);
  table->counts = (int*)malloc((table->capacity << 1)*sizeof(int));

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
  free(table->key_ptr);
  free(table->vals);
  free(table->keys);
  free(table->counts);
  free(table);
}
