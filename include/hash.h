#ifndef __HASH_H__
#define __HASH_H__

#include "carat/matrix.h"

typedef uint64_t W64;

struct hash_table_t {
  matrix_TYP** keys;
  int* vals;
  int* key_ptr;
  int* counts; 

  W64 num_stored;
  W64 mask; 
  W64 capacity;
};

typedef struct hash_table_t hash_table;

/* hash the form Q into an index between 0 and hash_size */
W64 hash_form(matrix_TYP* Q);

hash_table* create_hash(int hash_size);

void free_hash(hash_table* table);

int add(hash_table* table, matrix_TYP* key);

matrix_TYP* get_key(hash_table* table, matrix_TYP* key, int* index);

int exists(hash_table* table, matrix_TYP* key, int check_isom);

int indexof(hash_table* table, matrix_TYP* key, int check_isom);

void expand(hash_table* table);

int _add(hash_table* table, matrix_TYP* key, W64 val, int do_push_back);

int insert(hash_table* table, matrix_TYP* key, W64 val,
	   int index, int do_push_back);

#endif // __HASH_H__
