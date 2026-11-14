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
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Self dependencies

#include "arith.h"
#include "io.h"
#include "isometry.h"
#include "square_matrix.h"

// size of buffer to use when reading from files
#define BUF_SIZE 4096 // usually optimal for x86 architecture

// skip the white space at the beginning of a string
static const char* skip_spaces(const char* str)
{
  while (isspace((unsigned char)(*str)))
    str++;
  return str;
}

// parse a comma-separated array of integers, possibly enclosed in brackets.
// The array has to consist of exactly num_entries integers, and no more than
// num_entries entries are written to arr.
bool parse_array_int(Z64* arr, const char* arr_str, int num_entries)
{
  int idx;
  const char* str;
  char* end;
  long long val;
  bool has_bracket;

  idx = 0;
  str = skip_spaces(arr_str);
  has_bracket = (*str == '[');
  if (has_bracket)
    str++;
  while (true) {
    errno = 0;
    val = strtoll(str, &end, 10);
    // every entry has to be an integer which fits in 64 bits
    if ((end == str) || (errno == ERANGE))
      return false;
    // there is no room for this entry
    if (idx == num_entries)
      return false;
    arr[idx++] = val;
    str = skip_spaces(end);
    if (*str != ',')
      break;
    str++;
  }
  if (has_bracket) {
    if (*str != ']')
      return false;
    str = skip_spaces(str+1);
  }
  // nothing may follow the array
  return (*str == '\0') && (idx == num_entries);
}

// parse a symmetric 5x5 matrix, given as an array of 15 integers
bool parse_matrix(Z64* Q_coeffs, const char* mat_str)
{
  return parse_array_int(Q_coeffs, mat_str, 15);
}

// The parsers below read from the string *p_str, and move it past what they
// have read. When they fail, they may have moved it.

// read the character c, possibly after white space
static bool parse_char(const char** p_str, char c)
{
  const char* str;

  str = skip_spaces(*p_str);
  if (*str != c)
    return false;
  *p_str = str + 1;
  return true;
}

// read an integer which fits in 64 bits
static bool parse_int(const char** p_str, Z64* val)
{
  char* end;
  long long num;

  errno = 0;
  num = strtoll(*p_str, &end, 10);
  if ((end == *p_str) || (errno == ERANGE))
    return false;
  *val = num;
  *p_str = end;
  return true;
}

// read a rational number, whose denominator is positive and fits in an int
static bool parse_rat(const char** p_str, Z64* num, Z64* denom)
{
  if (!parse_int(p_str, num))
    return false;
  *denom = 1;
  if (!parse_char(p_str, '/'))
    return true;
  return parse_int(p_str, denom) && (*denom > 0) && (*denom <= INT_MAX);
}

// read a 5x5 matrix, given as an array of 5 arrays of 5 numbers. The numbers
// are integers if denom is NULL. Otherwise they are rational numbers, with
// the numerators read to num and the denominators to denom.
static bool parse_square_matrix(const char** p_str, square_matrix_t num,
				square_matrix_t denom)
{
  int i, j;

  if (!parse_char(p_str, '['))
    return false;
  for (i = 0; i < QF_RANK; i++) {
    if ((i > 0) && (!parse_char(p_str, ',')))
      return false;
    if (!parse_char(p_str, '['))
      return false;
    for (j = 0; j < QF_RANK; j++) {
      if ((j > 0) && (!parse_char(p_str, ',')))
	return false;
      if (denom == NULL) {
	if (!parse_int(p_str, &(num[i][j])))
	  return false;
      }
      else if (!parse_rat(p_str, &(num[i][j]), &(denom[i][j])))
	return false;
    }
    if (!parse_char(p_str, ']'))
      return false;
  }
  return parse_char(p_str, ']');
}

// read an isometry, given as a 5x5 matrix of rational numbers
static bool parse_isom(const char** p_str, isometry_t isom)
{
  square_matrix_t num, denom;
  Z64 isom_denom, mult, x, y;
  int i, j;

  if (!parse_square_matrix(p_str, num, denom))
    return false;

  // the common denominator has to fit in an int
  isom_denom = 1;
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++) {
      isom_denom *= denom[i][j] / gcdext(isom_denom, denom[i][j], &x, &y);
      if (isom_denom > INT_MAX)
	return false;
    }

  // and the numerators over it have to fit in 64 bits
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++) {
      mult = isom_denom / denom[i][j];
      if ((num[i][j] > INT64_MAX / mult) || (num[i][j] < INT64_MIN / mult))
	return false;
      num[i][j] *= mult;
    }

  // The denominators of isometries are multiplied when the isometries are
  // composed, so the common denominator is reduced.
  mult = isom_denom;
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      mult = gcdext(mult, num[i][j], &x, &y);
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < QF_RANK; j++)
      num[i][j] /= mult;
  isom_denom /= mult;

  isometry_init_set_square_matrix(isom, num, (int)isom_denom);
  return true;
}

// parse a lattice with an isometry to it, given as an array of the Gram matrix
// of the lattice, which is a symmetric 5x5 matrix of integers, and of the
// isometry, which is a 5x5 matrix of rational numbers
static bool parse_mat_and_isom(square_matrix_t Q, isometry_t isom,
			       const char* mat_str)
{
  const char* str;
  int i, j;

  str = mat_str;
  if ((!parse_char(&str, '[')) || (!parse_square_matrix(&str, Q, NULL)))
    return false;
  for (i = 0; i < QF_RANK; i++)
    for (j = 0; j < i; j++)
      if (Q[i][j] != Q[j][i])
	return false;
  if ((!parse_char(&str, ',')) || (!parse_isom(&str, isom)) ||
      (!parse_char(&str, ']')))
    return false;
  // nothing may follow the array
  return (*skip_spaces(str) == '\0');
}

// read the genus of the given discriminant from a file with the lattices of
// the genera of all discriminants, and with an isometry to each of the
// lattices if p_isom is not NULL. Exits with an error if the genus is not in
// the file, or is not given as expected.
static size_t read_genus_file(square_matrix_t** p_genus, isometry_t** p_isom,
			      const char* fname, size_t disc)
{
  FILE* genus_file;
  char buffer[BUF_SIZE];
  // the text of the lattice which is being read
  char lattice[BUF_SIZE];
  size_t num_read;
  size_t depth = 0;
  size_t i;
  size_t cur_disc = 0;
  Z64 coeffs[15];
  size_t genus_size = 0;
  size_t capacity = 0;
  size_t len = 0;
  bool is_done = false;
  bool is_found = false;
  bool is_valid;
  char c;

  *p_genus = NULL;
  if (p_isom != NULL)
    *p_isom = NULL;

  genus_file = fopen(fname, "r");
  if (genus_file == NULL) {
    fprintf(stderr, "Error: could not open the genus file %s\n", fname);
    exit(EXIT_FAILURE);
  }

  // The file is an array of genera, each of which is an array of lattices.
  // The brackets of the file are at depth 1, those of a genus at depth 2 and
  // those of a lattice at depth 3. The genus read is the one at index disc.
  while ((!is_done) &&
	 ((num_read = fread(buffer, sizeof(char), BUF_SIZE, genus_file)) > 0)) {
    for (i = 0; (i < num_read) && (!is_done); i++) {
      c = buffer[i];
      if (c == '[')
	depth++;
      if (cur_disc == disc) {
	// the text of a lattice is everything between its brackets, with them
	if (depth >= 3) {
	  if (len == BUF_SIZE - 1) {
	    fprintf(stderr, "Error: in the genus file %s, a lattice of discriminant %zu is longer than %d characters\n",
		    fname, disc, BUF_SIZE - 1);
	    exit(EXIT_FAILURE);
	  }
	  lattice[len++] = c;
	}
	// between the lattices of a genus there are only commas and white space
	else if ((depth == 2) && (c != '[') && (c != ']') && (c != ',') &&
		 (!isspace((unsigned char)c))) {
	  fprintf(stderr, "Error: in the genus file %s, the character '%c' is between the lattices of discriminant %zu\n",
		  fname, c, disc);
	  exit(EXIT_FAILURE);
	}
      }
      if (c == ']') {
	// the file has ended, or has no array of genera
	if (depth <= 1)
	  is_done = true;
	else {
	  depth--;
	  // the lattice has ended
	  if ((depth == 2) && (cur_disc == disc)) {
	    lattice[len] = '\0';
	    len = 0;
	    if (genus_size == capacity) {
	      capacity = (capacity == 0) ? 16 : 2 * capacity;
	      *p_genus = (square_matrix_t*)realloc(*p_genus, capacity * sizeof(square_matrix_t));
	      if (p_isom != NULL)
		*p_isom = (isometry_t*)realloc(*p_isom, capacity * sizeof(isometry_t));
	    }
	    if (p_isom != NULL)
	      is_valid = parse_mat_and_isom((*p_genus)[genus_size], (*p_isom)[genus_size], lattice);
	    else
	      is_valid = parse_matrix(coeffs, lattice);
	    if (!is_valid) {
	      fprintf(stderr, "Error: in the genus file %s, the lattice %s of discriminant %zu is not %s\n",
		      fname, lattice, disc, (p_isom != NULL) ?
		      "a symmetric 5x5 matrix of integers with a 5x5 matrix of rational numbers" :
		      "a list of 15 integers");
	      exit(EXIT_FAILURE);
	    }
	    // at the moment genus files are in sage (=GG) format
	    if (p_isom == NULL)
	      square_matrix_init_set_symm((*p_genus)[genus_size], coeffs, "GG");
	    genus_size++;
	  }
	  // the genus has ended
	  if (depth == 1) {
	    is_found = (cur_disc == disc);
	    is_done = is_found;
	    cur_disc++;
	  }
	}
      }
    }
  }

  fclose(genus_file);

  if ((!is_found) && (cur_disc == disc) && (depth >= 2)) {
    fprintf(stderr, "Error: the genus file %s ends inside the genus of discriminant %zu\n",
	    fname, disc);
    exit(EXIT_FAILURE);
  }
  if (!is_found) {
    fprintf(stderr, "Error: the genus file %s has no genus of discriminant %zu\n",
	    fname, disc);
    exit(EXIT_FAILURE);
  }

  return genus_size;
}

// read the genus from a file with all the genus representatives for all discriminants
size_t read_genus(square_matrix_t** p_genus, const char* fname, size_t disc)
{
  return read_genus_file(p_genus, NULL, fname, disc);
}

// read the genus and the isometries between the genus representatives from a file
// containing the genus representatives and the isometries for all discriminants.
size_t read_genus_and_isom(square_matrix_t** p_genus,
			   isometry_t** p_isom,
			   const char* fname, size_t disc)
{
  return read_genus_file(p_genus, p_isom, fname, disc);
}
