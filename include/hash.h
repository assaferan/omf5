#ifndef __HASH_H__
#define __HASH_H__

#include <stdint.h>

#include <carat/matrix.h>
#include <flint/fmpq.h>

#include "isometry.h"
#include "square_matrix.h"
#include "typedefs.h"

// We do not anticipate more than 2^32 forms in the genus...
// Maybe changing to W16 could improve that, have to check

typedef W32 hash_t;

typedef struct {
  square_matrix_t* keys;
  hash_t* vals;
  hash_t* key_ptr;
  W32* counts;
  fmpq_t* probs;

  hash_t num_stored;
  hash_t mask; 
  hash_t capacity;
  W32 theta_prec;

  bool red_on_isom;
} hash_table;

typedef hash_table hash_table_t[1];

/* hash the form Q into an index between 0 and hash_size */
hash_t hash_form(const square_matrix_t Q, W32 theta_prec);

void hash_table_init(hash_table_t table, hash_t hash_size);

void hash_table_recalibrate(hash_table_t new_table, const hash_table_t old_table);

void hash_table_clear(hash_table_t table);

int hash_table_add(hash_table_t table, const square_matrix_t key);

bool hash_table_get_key(square_matrix_t val, const hash_table_t table, const square_matrix_t key, int* index);

int hash_table_exists(const hash_table_t table, const square_matrix_t key, int check_isom);

int hash_table_indexof(const hash_table_t table, const square_matrix_t key, int check_isom, double* theta_time, double* isom_time, int* num_isom);

int hash_table_index_and_isom(const hash_table_t table, const square_matrix_t key, isometry_t isom,
			      double* theta_time, double* isom_time, int* num_isom);

void hash_table_expand(hash_table_t table);

int hash_table_insert(hash_table_t table, const square_matrix_t key, hash_t val,
		      int index, int do_push_back);

#endif // __HASH_H__
