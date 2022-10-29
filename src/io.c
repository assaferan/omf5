#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "io.h"
#include "square_matrix.h"

#define BUF_SIZE 4096 // usually optimal for x86 architecture

// can change this one to be dynamically allocated, but it's small enough
#define MAX_GEN_SIZE 2000 

bool parse_matrix(const char* mat_str, int* Q_coeffs)
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
    Q_coeffs[idx++] = atoi(token);
    token = strsep(&cpy, ",");
  }
  free(original);
  return (idx == 15);
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
  int coeffs[MAX_GEN_SIZE][15];
  size_t genus_size = 0;
  size_t mat_buf_idx = 0;
  bool mat_read;

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
	      mat_read = parse_matrix(matrix_buffer, coeffs[genus_size++]);
	      assert(mat_read);
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
	      mat_read = parse_matrix(matrix_buffer, coeffs[genus_size++]);
	      assert(mat_read);
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
