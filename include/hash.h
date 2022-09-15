#ifndef __HASH_H__
#define __HASH_H__

#include "carat/matrix.h"
#include "flint/fmpq.h"

typedef uint64_t W64;
typedef uint32_t W32;
typedef uint16_t W16;

// We do not anticipate more than 2^32 forms in the genus...
// Maybe changing to W16 could improve that, have to check

typedef W32 hash_t;

struct hash_table_t {
  matrix_TYP** keys;
  hash_t* vals;
  hash_t* key_ptr;
  W32* counts;
  fmpq_t* probs;

  hash_t num_stored;
  hash_t mask; 
  hash_t capacity;
  W32 theta_prec;

  int red_on_isom;
};

typedef struct hash_table_t hash_table;

/* hash the form Q into an index between 0 and hash_size */
hash_t hash_form(matrix_TYP* Q, W32 theta_prec);

hash_table* create_hash(hash_t hash_size);

hash_table* recalibrate_hash(hash_table* table);

void free_hash(hash_table* table);

int add(hash_table* table, matrix_TYP* key);

matrix_TYP* get_key(hash_table* table, matrix_TYP* key, int* index);

int exists(hash_table* table, matrix_TYP* key, int check_isom);

int indexof(const hash_table* table, matrix_TYP* key, int check_isom, double* theta_time, double* isom_time, int* num_isom);

void expand(hash_table* table);

int _add(hash_table* table, matrix_TYP* key, hash_t val, int do_push_back);

int insert(hash_table* table, matrix_TYP* key, hash_t val,
	   int index, int do_push_back);

#endif // __HASH_H__
