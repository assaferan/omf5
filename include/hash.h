#ifndef __HASH_H__
#define __HASH_H__

#include <stdint.h>

#include <carat/matrix.h>
#include <flint/fmpq.h>

#include "typedefs.h"

// We do not anticipate more than 2^32 forms in the genus...
// Maybe changing to W16 could improve that, have to check

typedef W32 hash_t;

typedef struct {
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
} hash_table;

typedef hash_table hash_table_t[1];

/* hash the form Q into an index between 0 and hash_size */
hash_t hash_form(matrix_TYP* Q, W32 theta_prec);

void hash_table_init(hash_table_t table, hash_t hash_size);

void hash_table_recalibrate(hash_table_t new_table, hash_table_t old_table);

void hash_table_clear(hash_table_t table);

int hash_table_add(hash_table_t table, matrix_TYP* key);

matrix_TYP* hash_table_get_key(const hash_table_t table, matrix_TYP* key, int* index);

int hash_table_exists(const hash_table_t table, matrix_TYP* key, int check_isom);

int hash_table_indexof(const hash_table_t table, matrix_TYP* key, int check_isom, double* theta_time, double* isom_time, int* num_isom);

void hash_table_expand(hash_table_t table);

int hash_table_insert(hash_table_t table, matrix_TYP* key, hash_t val,
		      int index, int do_push_back);

#endif // __HASH_H__
