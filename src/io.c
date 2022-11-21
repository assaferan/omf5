#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "arith.h"
#include "io.h"
#include "isometry.h"
#include "square_matrix.h"

#define BUF_SIZE 4096 // usually optimal for x86 architecture

// can change this one to be dynamically allocated, but it's small enough
#define MAX_GEN_SIZE 2000 

// parse a comma-separated array of integers
bool parse_array_int(Z64* arr, const char* arr_str, int num_entries)
{
  int idx, len, offset;
  char *original;
  char *cpy;
  char *token;

  idx = 0;
  offset = 0;
  len = strlen(arr_str);
  if ((len > 0) && (arr_str[0] == '[')) {
    assert(arr_str[len-1] == ']');
    len -= 2;
    offset = 1;
  }
  original = (char*)malloc((len+1)*sizeof(char));
  memcpy(original, arr_str+offset, len+1);
  cpy = original;
  token = strsep(&cpy, ",");
  while(token != NULL) {
    arr[idx++] = atoi(token);
    token = strsep(&cpy, ",");
  }
  free(original);
  return (idx == num_entries);
}

// parse a comma-separated array of integers
bool parse_array_rat(Z64* num, int* denom,
		     const char* arr_str, int num_entries)
{
  int idx, len, offset;
  char *original;
  char *cpy;
  char *token;

  idx = 0;
  offset = 0;
  len = strlen(arr_str);
  if ((len > 0) && (arr_str[0] == '[')) {
    assert(arr_str[len-1] == ']');
    len -= 2;
    offset = 1;
  }
  original = (char*)malloc((len+1)*sizeof(char));
  memcpy(original, arr_str+offset, len+1);
  cpy = original;
  token = strsep(&cpy, "/");
  while(token != NULL) {
    num[idx] = atoi(token);
    token = strsep(&cpy, ",");
    denom[idx++] = atoi(token);
    token = strsep(&cpy, "/");
  }
  free(original);
  return (idx == num_entries);
}

bool parse_matrix(Z64* Q_coeffs, const char* mat_str)
{
  return parse_array_int(Q_coeffs, mat_str, 15);
}

bool parse_int_matrix_full(square_matrix_t Q, const char* mat_str)
{
  int idx, len;
  char *original;
  char *cpy;
  char *token;

  idx = 0;
  len = strlen(mat_str);
  original = (char*)malloc((len+1)*sizeof(char));
  memcpy(original, mat_str, len+1);
  cpy = original;
  token = strsep(&cpy, ",");
  while(token != NULL) {
    parse_array_int(Q[idx++], token, QF_RANK);
    token = strsep(&cpy, ",");
  }
  free(original);
  return (idx == QF_RANK);
}

bool parse_isom(isometry_t isom, const char* mat_str)
{
  int idx, len, i, j;
  char *original;
  char *cpy;
  char *token;
  int denom[QF_RANK];
  square_matrix_t isom_num;
  int old_isom_den;
  int isom_den;

  idx = 0;
  len = strlen(mat_str);
  original = (char*)malloc((len+1)*sizeof(char));
  memcpy(original, mat_str, len+1);
  cpy = original;
  token = strsep(&cpy, ",");
  isom_den = 1;
  while(token != NULL) {
    parse_array_rat(isom_num[idx], denom, token, QF_RANK);
    old_isom_den = isom_den;
    for (j = 0; j < QF_RANK; j++)
      isom_den = lcm(isom_den, denom[j]);
    for (j = 0; j < QF_RANK; j++)
      isom_num[idx][j] *= (isom_den / denom[j]);
    for (i = 0; i < idx; i++)
      for (j = 0; j < QF_RANK; j++)
	isom_num[i][j] *= (isom_den / old_isom_den);
    idx++;
    token = strsep(&cpy, ",");
  }
  isometry_init_set_square_matrix(isom, isom_num, isom_den);
  free(original);
  return (idx == QF_RANK);
}

bool parse_mat_and_isom(square_matrix_t Q, isometry_t isom,
			const char* mat_str)
{
  bool success;
  int len;
  char *original;
  char *isom_str;
  char *intmat_str;

  len = strlen(mat_str);
  original = (char*)malloc((len+1)*sizeof(char));
  memcpy(original, mat_str, len+1);
  isom_str = original;
  intmat_str = strsep(&isom_str, "]");
  success = parse_int_matrix_full(Q, intmat_str);
  if (!success) return false;
  success = parse_isom(isom, isom_str);
  return success;
}

size_t read_genus(square_matrix_t** p_genus, const char* fname, size_t disc)
{
  FILE* genus_file;
  char buffer[BUF_SIZE];
  // this one can be smaller
  char matrix_buffer[BUF_SIZE];
  size_t num_read = BUF_SIZE;
  size_t depth = 0;
  size_t i;
  size_t cur_disc = 0;
  Z64 coeffs[MAX_GEN_SIZE][15];
  size_t genus_size = 0;
  size_t mat_buf_idx = 0;

  genus_file = fopen(fname, "r");

  while ((num_read == BUF_SIZE) && (cur_disc <= disc) ){
    num_read = fread(buffer, sizeof(char), BUF_SIZE, genus_file);
    for (i = 0; (i < num_read) && (cur_disc <= disc); i++) {
      switch(buffer[i]) {
      case '[':
	depth++;
	break;
      case ']':
	depth--;
	if (depth == 1) {
	  if (cur_disc == disc) {
	    if (mat_buf_idx != 0) {
	      matrix_buffer[mat_buf_idx] = '\0';
	      parse_matrix(coeffs[genus_size++], matrix_buffer);
	      mat_buf_idx = 0;
	    }
	  }
	  cur_disc++;
	}
	break;
      default:
	if (cur_disc == disc) {
	  switch(depth) {
	  case 2:
	    if (mat_buf_idx != 0) {
	      matrix_buffer[mat_buf_idx] = '\0';
	      parse_matrix(coeffs[genus_size++], matrix_buffer);
	      mat_buf_idx = 0;
	    }
	    break;
	  case 3:
	    matrix_buffer[mat_buf_idx++] = buffer[i];
	    break;
	  default:
	    break;
	  }
	}
	break;
      }
    }
  }

  fclose(genus_file);

  *p_genus = (square_matrix_t*)malloc(genus_size * sizeof(square_matrix_t));

  for (i = 0; i < genus_size; i++) {
    // at the moment genus files are in sage (=GG) format
    square_matrix_init_set_symm((*p_genus)[i], coeffs[i], "GG");
  }
  
  return genus_size;
}

size_t read_genus_and_isom(square_matrix_t** p_genus,
			   isometry_t** p_isom,
			   const char* fname, size_t disc)
{
  FILE* genus_file;
  char buffer[BUF_SIZE];
  // this one can be smaller
  char matrix_buffer[BUF_SIZE];
  size_t num_read = BUF_SIZE;
  size_t depth = 0;
  size_t i;
  size_t cur_disc = 0;
  square_matrix_t mats[MAX_GEN_SIZE];
  isometry_t isoms[MAX_GEN_SIZE];
  size_t genus_size = 0;
  size_t mat_buf_idx = 0;

  genus_file = fopen(fname, "r");

  while ((num_read == BUF_SIZE) && (cur_disc <= disc) ){
    num_read = fread(buffer, sizeof(char), BUF_SIZE, genus_file);
    for (i = 0; (i < num_read) && (cur_disc <= disc); i++) {
      switch(buffer[i]) {
      case '[':
	depth++;
	break;
      case ']':
	depth--;
	if (depth == 1) {
	  if (cur_disc == disc) {
	    if (mat_buf_idx != 0) {
	      matrix_buffer[mat_buf_idx] = '\0';
	      parse_mat_and_isom(mats[genus_size], isoms[genus_size],
				 matrix_buffer);
	      genus_size++;
	      mat_buf_idx = 0;
	    }
	  }
	  cur_disc++;
	}
	break;
      default:
	if (cur_disc == disc) {
	  switch(depth) {
	  case 2:
	    if (mat_buf_idx != 0) {
	      matrix_buffer[mat_buf_idx] = '\0';
	      parse_mat_and_isom(mats[genus_size], isoms[genus_size],
				 matrix_buffer);
	      genus_size++;
	      mat_buf_idx = 0;
	    }
	    break;
	  case 3:
	    matrix_buffer[mat_buf_idx++] = buffer[i];
	    break;
	  default:
	    break;
	  }
	}
	break;
      }
    }
  }

  fclose(genus_file);

  *p_genus = (square_matrix_t*)malloc(genus_size * sizeof(square_matrix_t));
  *p_isom = (isometry_t*)malloc(genus_size * sizeof(isometry_t));

  for (i = 0; i < genus_size; i++) {
    square_matrix_init((*p_genus)[i]);
    square_matrix_set((*p_genus)[i], mats[i]);
    isometry_init_set((*p_isom)[i], isoms[i]);
  }
  
  return genus_size;
}
