#include <assert.h>

#include <carat/typedef.h>

#include "io.h"
#include "tests.h"

#define MAX_STR_LEN 256

bool handle_flag_int(const char* flag_name, const char* param_str, int* flag_val);
bool handle_flag_int_seq(const char* flag_name, const char* param_str, int** flag_val, int* num_vals);
int print_param_desc(char* argv[]);

int main(int argc, char* argv[])
{
  int prec, p, c, disc;
  Z64 Q_coeffs[15];
  STATUS test_res;
  char* input_type;
  char* genus_fname;
  int max_args;
  bool do_tests, is_valid, is_prec, is_format, is_p, is_c, is_genus, is_disc, is_idxs;
  bool has_quad, has_format, has_prec, has_p, has_hecke, has_c;
  bool has_genus, has_disc, has_large, has_row, has_isom, has_nonlifts, has_idxs;
  int i, j;
  genus_t genus;
  int* idxs;
  int num_idxs;

  idxs = NULL;
  input_type = NULL;
  genus_fname = NULL;
  max_args = 8;
  num_idxs = 0;
  
  if ((argc > max_args) || (argc == 1)) {
    // print correct usage
    return print_param_desc(argv);
  }

  do_tests = false;
  has_quad = has_format = has_prec  = has_p = has_hecke = has_c = false;
  has_row = has_genus = has_disc = has_large = has_isom = has_nonlifts = has_idxs = false;
  
  for (i = 1; i < argc; i++) {
    is_valid = false;
    
    // checkng to see if the flag is -tests
    if (strcmp(argv[i], "-tests") == 0) {
      do_tests = true;
      is_valid = true;
    }

    if (strcmp(argv[i], "-hecke") == 0) {
      has_hecke = true;
      is_valid = true;
    }

    if (strcmp(argv[i], "-isom") == 0) {
      has_isom = true;
      is_valid = true;
    }

    if (strcmp(argv[i], "-nonlifts") == 0) {
      has_nonlifts = true;
      is_valid = true;
    }
    
    if (strcmp(argv[i], "-large") == 0) {
      has_large = true;
      is_valid = true;
    }

    if (strcmp(argv[i], "-row") == 0) {
      has_row = true;
      is_valid = true;
    }

    // checking whether this is a matrix input
    if (strncmp(argv[i], "-quad=", 6) == 0) {
      is_valid = parse_matrix(Q_coeffs, argv[i]+6);
      has_quad = is_valid;
    }

    // checking for format of the matrix input
    if (strncmp(argv[i], "-format=", 8) == 0) {
      input_type = argv[i]+8;
      is_format = strncmp(input_type,"GG",2) || strncmp(input_type, "A", 1);
      if (is_format)
	has_format = true;
      is_valid = (is_valid) || (is_format);
    }

    // checking for a file with the genus representatives
    if (strncmp(argv[i], "-genus=", 7) == 0) {
      genus_fname = argv[i]+7;
      is_genus = true;
      has_genus = true;
      is_valid = (is_valid) || (is_genus);
    }

    is_prec = handle_flag_int("prec", argv[i], &prec);
    is_valid = (is_valid) || (is_prec);
    if (is_prec)
      has_prec = true;

    is_p = handle_flag_int("p", argv[i], &p);
    is_valid = (is_valid) || (is_p);
    if (is_p)
      has_p = true;

    is_c = handle_flag_int("cond", argv[i], &c);
    is_valid = (is_valid) || (is_c);
    if (is_c)
      has_c = true;

    is_disc = handle_flag_int("disc", argv[i], &disc);
    is_valid = (is_valid) || (is_disc);
    if (is_disc)
      has_disc = true;

    is_idxs = handle_flag_int_seq("idxs", argv[i], &idxs, &num_idxs);
    is_valid = (is_valid) || (is_idxs);
    if (is_idxs) {
      has_idxs = true;
      printf("Read %d indices. These are their values:", num_idxs);
      for (j = 0; j < num_idxs; j++)
	printf("%d,", idxs[j]);
      printf("\b\n");
    }

    if (!is_valid)
      return print_param_desc(argv);
  }

  test_res = SUCCESS;
  
  if (do_tests) {
    test_res <<= 1;
    test_res |= test_greedy_overflow();
    test_res <<= 1;
    test_res |= test_61();
    test_res <<= 1;
    test_res |= test_69();
  }

  if (has_quad && has_format)
    compute_genus(genus, Q_coeffs, input_type);
  else if (has_genus && has_disc)
    genus_init_file(genus, genus_fname, disc, has_isom);

  has_genus = (has_quad && has_format) || (has_genus && has_disc);

  if (has_genus) {
    test_res <<= 1;
    if (has_prec)
      test_res |= compute_eigenvalues_up_to(genus, prec, has_nonlifts, num_idxs, idxs);
    else
      if (has_p)
	if (has_hecke) {
	  if (has_row) {
	    // !! TODO - right now it doesn't matter which column we take,
	    // so we take 0 index,
	    // might change in the future.
	    if (has_c)
	      test_res |= compute_hecke_col(genus, p, c);
	    else
	      test_res |= compute_hecke_col_all_conds(genus, p, 0);
	  } else {
	    if (has_c)
	      test_res |= compute_hecke_matrix(genus, p, c);
	    else
	      test_res |= compute_hecke_matrix_all_conds(genus, p);
	  }
	}
	else
	  test_res |= compute_eigenvalues(genus, p, has_nonlifts, num_idxs, idxs);
      else {
	if (has_hecke)
	  if (has_large)
	    test_res |= compute_first_large_hecke(genus);
	  else
	    test_res |= compute_first_hecke_matrix_all_conds(genus);
	else
	  test_res |= compute_eigenvectors(genus, has_nonlifts, num_idxs, idxs);
      }
  }
  else {
    if (!do_tests)
      return print_param_desc(argv);
  }

  if (has_genus)
    genus_clear(genus);

  if (has_idxs)
    free(idxs);
  
  return test_res;
}


bool handle_flag_int(const char* flag_name, const char* param_str, int* flag_val)
{
  size_t flag_len;
  char* full_flag_name;
  
  flag_len = strlen(flag_name) + 2;
  full_flag_name = (char*)malloc((flag_len+1)*sizeof(char));

  if (full_flag_name == NULL) {
    printf("Error, flag name is too long!\n");
    return false;
  }

  strncpy(full_flag_name, "-", 2);
  strncat(full_flag_name, flag_name, MAX_STR_LEN);
  strncat(full_flag_name, "=", 2);

  assert(flag_len == strlen(full_flag_name));
  
  if (strncmp(param_str, full_flag_name, flag_len) == 0) {
    *flag_val = atoi(param_str + flag_len);
    free(full_flag_name);
    return true;
  }

  free(full_flag_name);
  return false;
}

bool handle_flag_int_seq(const char* flag_name, const char* param_str, int** flag_val, int* num_vals)
{
  size_t flag_len;
  char* full_flag_name;
  char number_str[8];
  size_t str_idx, num_idx;
  char next_ch;
  size_t INITIAL_NUM_VALS_ALLOC = 8;
  size_t num_allocated = 0;
  
  flag_len = strlen(flag_name) + 2;
  full_flag_name = (char*)malloc((flag_len+1)*sizeof(char));

  if (full_flag_name == NULL) {
    printf("Error, flag name is too long!\n");
    return false;
  }

  strncpy(full_flag_name, "-", 2);
  strncat(full_flag_name, flag_name, MAX_STR_LEN);
  strncat(full_flag_name, "=", 2);

  assert(flag_len == strlen(full_flag_name));

  str_idx = flag_len;
  if (strncmp(param_str, full_flag_name, flag_len) == 0) {
    *num_vals = 0;
    *flag_val = (int *)malloc(INITIAL_NUM_VALS_ALLOC*sizeof(int));
    num_allocated = INITIAL_NUM_VALS_ALLOC;
    next_ch = param_str[str_idx];
    // while the string didn't end
    while (next_ch != '\0') {
      // we read off the next number
      num_idx = 0;
      while ((next_ch != ',') && (next_ch != '\0')) {
	number_str[num_idx] = next_ch;
	str_idx++;
	num_idx++;
	next_ch = param_str[str_idx];
      }
      number_str[num_idx] = '\0';
      // skipping the comma
      if (next_ch == ',') {
	str_idx++;
	next_ch = param_str[str_idx];
      }
      assert((*num_vals) < num_allocated);
      (*flag_val)[(*num_vals)] = atoi(number_str);
      (*num_vals)++;
      if ((*num_vals) >= num_allocated) {
	num_allocated *= 2;
	*flag_val = (int *)realloc(*flag_val, num_allocated*sizeof(int));
      }
    }
    free(full_flag_name);
    return true;
  }

  free(full_flag_name);
  return false;
}

int print_param_desc(char* argv[])
{
  fprintf(stderr, "Usage: %s [-tests] [-quad=Q] [-format=f] [-prec=L] [-hecke] [-row] [-p=p] [-genus=g] [-isom] [-disc=d] [-cond=c] \n", argv[0]);
  fprintf(stderr, "The genus can be specified in one of two ways. Either via Q and f - \n");
  fprintf(stderr, "[Q] is the quinary quadratic form (lattice) given as 15 comma-separated integers in a format specified by f,\n");
  fprintf(stderr, "[f] is either 'GG' or 'A', the former for the format in Rama-Toranria webpage, the latter for the Magma format,\n");
  fprintf(stderr, "in which case, the genus will be computed using p-neighbors, or via g and d - \n");
  fprintf(stderr, "[g] is the name of a file containing the list of genera,\n");
  fprintf(stderr, "[d] is the discriminant of the lattice, so that g[d] is the relevant genus,\n");
  fprintf(stderr, "[p] is a prime at which to compute the Hecke matrix/eigenvalue, \n"); 
  fprintf(stderr, "[L] is the preicision up to which to compute the hecke matrices/eigenvalues (a_p for p <= L and a_{p^2} for p^2 <= L),\n");
  fprintf(stderr, "[c] is the conductor of the spinor norm character. If not specified, the program will compute all of them. At the moment, only relevant for computing a column of the Hecke matrix.");
  fprintf(stderr, "If the flag -hecke is supplied, computes a column of the Hecke matrix, otherwise computes the Hecke eigenvalues of forms that are non-lifts. If p is not supplied, computes the Hecke matrix of the first prime not dividing the discriminant.\n");
  fprintf(stderr, "If either L or p is not supplied, only decomposes the space, and finds Hecke eigenvectors.\n");
  fprintf(stderr, "If the flag -row is supplied in addition to -hecke, computes a single row.\n");
  fprintf(stderr, "If the flag -isom is supplied in addition to -genus=g, the genus includes the isometries.\n");
  fprintf(stderr, "If the flag -tests is supplied, additionally runs standard tests.\n");
  return -1;
}
