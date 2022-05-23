#ifndef __HASH_H__
#define __HASH_H__

#include "carat/matrix.h"

struct hash_table_t {
  matrix_TYP** keys;
  int* vals;
  int* key_ptr;

  int num_stored;
  int mask; // uint64_t ?
  int capacity;
};

typedef struct hash_table_t hash_table;

/* hash the form Q into an index between 0 and hash_size */
int hash_form(matrix_TYP* Q);

hash_table* create_hash(int hash_size);

void free_hash(hash_table* table);

int add(hash_table* table, matrix_TYP* key);

matrix_TYP* get_key(hash_table* table, matrix_TYP* key, int* index);

int exists(hash_table* table, matrix_TYP* key);

int indexof(hash_table* table, matrix_TYP* key);

void expand(hash_table* table);

int _add(hash_table* table, matrix_TYP* key, int val, int do_push_back);

int insert(hash_table* table, matrix_TYP* key, int val,
	   int index, int do_push_back);

#endif // __HASH_H__
