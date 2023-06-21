/*****************************************************************
 *
 * Package : omf5 - orthogonal modular forms of rank 5
 * Filename : io.c
 *
 * Description: Functions for handling input and output
 *
 *****************************************************************
 */

// System dependencies

#include <assert.h>
#include <stdio.h>
#include <string.h>

// Self dependencies

#include "arith.h"
#include "io.h"
#include "isometry.h"
#include "square_matrix.h"

// size of buffer to use when reading from files
#define BUF_SIZE 4096 // usually optimal for x86 architecture

// maximal number of genus representatives (class number)
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

// parse a comma-separated array of ratinoal numbers
bool parse_array_rat(Z64* num, int* denom,
		     const char* arr_str, int num_entries)
{
  int idx, len, offset;
  char *original;
  char *cpy;
  char *token;
  char *slash_pos, *comma_pos;

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
  original[len] = '\0';
  cpy = original;
  token = cpy;
  while(token != NULL) {
    slash_pos = strstr(cpy, "/");
    comma_pos = strstr(cpy, ",");
    if ((slash_pos != NULL) && (comma_pos != NULL) && (slash_pos < comma_pos)) {
      token[slash_pos - cpy] = '\0';
      num[idx] = atoi(token);
      slash_pos++;
      slash_pos[comma_pos - slash_pos] = '\0';
      denom[idx++] = atoi(slash_pos);
      cpy = comma_pos + 2;
    }
    else if (comma_pos!= NULL) { 
      token[comma_pos-cpy] = '\0';
      num[idx] = atoi(token);
      denom[idx++] = 1;
      cpy = comma_pos + 2;
    } else if (slash_pos != NULL) {
      token[slash_pos - cpy] = '\0';
      num[idx] = atoi(token);
      slash_pos++;
      denom[idx++] = atoi(slash_pos);
      cpy = NULL;
    } else {
      num[idx] = atoi(token);
      denom[idx++] = 1;
      cpy = NULL;
    }
    token = cpy;
  }
  free(original);
  return (idx == num_entries);
}

// parse a symmetric 5x5 matrix, given as an array of 15 integers
bool parse_matrix(Z64* Q_coeffs, const char* mat_str)
{
  return parse_array_int(Q_coeffs, mat_str, 15);
}

// parse a 5x5 matrix of integers, given as an array of 5 arrays of integers 
bool parse_int_matrix_full(square_matrix_t Q, const char* mat_str)
{
  int idx, len;
  char *original;
  char *cpy;
  char *token;

  idx = 0;
  len = strlen(mat_str)-2;
  original = (char*)malloc((len+1)*sizeof(char));
  memcpy(original, mat_str+1, len);
  original[len] = '\0';
  cpy = original;
  token = cpy;
  while(token != NULL) {
    cpy = strstr(cpy, "],");
    if (cpy != NULL) {
      token[cpy-token+1] = '\0';
      cpy += 3;
    }
    parse_array_int(Q[idx++], token, QF_RANK);
    token = cpy;
  }
  free(original);
  return (idx == QF_RANK);
}

// parse an isometry (a 5x5 matrix of rationals), given as an array of 5 arrays of rationals) 
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
  bool success;

  idx = 0;
  len = strlen(mat_str)-2;
  original = (char*)malloc((len+1)*sizeof(char));
  memcpy(original, mat_str+1, len);
  original[len] = '\0';
  cpy = original;
  token = cpy;
  isom_den = 1;
  while(token != NULL) {
    cpy = strstr(cpy, "],");
    if (cpy != NULL) {
      token[cpy-token+1] = '\0';
      cpy += 3;
    }
    success = parse_array_rat(isom_num[idx], denom, token, QF_RANK);
    if (!success) return false;
    old_isom_den = isom_den;
    for (j = 0; j < QF_RANK; j++)
      isom_den = lcm(isom_den, denom[j]);
    for (j = 0; j < QF_RANK; j++)
      isom_num[idx][j] *= (isom_den / denom[j]);
    for (i = 0; i < idx; i++)
      for (j = 0; j < QF_RANK; j++)
	isom_num[i][j] *= (isom_den / old_isom_den);
    idx++;
    token = cpy;
  }
  isometry_init_set_square_matrix(isom, isom_num, isom_den);
  free(original);
  return (idx == QF_RANK);
}

// parse a pair of a matrix and an isometry from input
bool parse_mat_and_isom(square_matrix_t Q, isometry_t isom,
			const char* mat_str)
{
  bool success;
  int len;
  char *isom_str;
  char *intmat_str;

  len = strlen(mat_str)-2;
  intmat_str = (char*)malloc((len+1)*sizeof(char));
  memcpy(intmat_str, mat_str+1, len);
  intmat_str[len] = '\0';
  isom_str = strstr(intmat_str, "],[") + 2;
  intmat_str[isom_str-intmat_str-1] = '\0';
  success = parse_int_matrix_full(Q, intmat_str);
  if (!success) return false;
  success = parse_isom(isom, isom_str);
  free(intmat_str);
  return success;
}

// read the genus from a file with all the genus representatives for all discriminants
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

// read the genus and the isometries between the genus representatives from a file
// containing the genus representatives and the isometries for all discriminants.
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
	if (depth >= 3) {
	  if (cur_disc == disc) {
	    matrix_buffer[mat_buf_idx++] = buffer[i];
	  }
	}
	break;
      case ']':
	depth--;
	if (depth >= 2) {
	  if (cur_disc == disc) {
	    matrix_buffer[mat_buf_idx++] = buffer[i];
	  }
	}
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
	  case 4:
	  case 5:
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
