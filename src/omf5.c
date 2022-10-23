#include <assert.h>

#include <carat/typedef.h>

#include "tests.h"

#define MAX_STR_LEN 256

bool handle_flag_int(const char* flag_name, const char* param_str, int* flag_val);
int print_param_desc(char* argv[]);
int parse_matrix(const char* mat_str, int* Q_coeffs);

int main(int argc, char* argv[])
{
  int form_idx, prec, p;
  int Q_coeffs[15];
  STATUS test_res;
  char* input_type;
  int max_args;
  bool do_tests, is_valid, is_prec, is_format, is_form_idx, is_p;
  bool has_quad, has_format, has_prec, has_form_idx, has_p;
  int i;

  input_type = NULL;
  max_args = 6;
  
  if ((argc > max_args) || (argc == 1)) {
    // print correct usage
    return print_param_desc(argv);
  }

  do_tests = false;
  has_quad = has_format = has_prec = has_form_idx = false;
  
  for (i = 1; i < argc; i++) {
    is_valid = false;
    
    // checkng to see if the flag is -tests
    if (strcmp(argv[i], "-tests") == 0) {
      do_tests = true;
      is_valid = true;
    }

    // checking whether this is a matrix input
    if (strncmp(argv[i], "-quad=", 6) == 0) {
      is_valid = parse_matrix(argv[i]+6, Q_coeffs);
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

    is_prec = handle_flag_int("prec", argv[i], &prec);
    is_valid = (is_valid) || (is_prec);
    if (is_prec)
      has_prec = true;

    is_p = handle_flag_int("p", argv[i], &p);
    is_valid = (is_valid) || (is_p);
    if (is_p)
      has_p = true;
    
    is_form_idx = handle_flag_int("form_idx", argv[i], &form_idx);
    is_valid = (is_valid) || (is_form_idx);
    if (is_form_idx)
      has_form_idx = true;

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
    // return test_res;
  }

  if (has_quad && has_format) {
    test_res <<= 1;
    if (has_prec && has_form_idx)
      test_res |= compute_eigenvalues_up_to(Q_coeffs,form_idx,
					    prec, input_type);
    else
      if (has_p)
	test_res |= compute_eigenvalues(Q_coeffs, p, input_type);
      else
	test_res |= compute_eigenvectors(Q_coeffs, input_type);
  }
  else
    if (!do_tests)
      return print_param_desc(argv);
  
  
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

int print_param_desc(char* argv[])
{
  printf("Usage: %s [-tests] [-quad=Q] [-format=f] [-prec=L] [-form_idx=idx]\n", argv[0]);
  printf("[Q] is the quinary quadratic form (lattice) given as 15 comma-separated integers in a format specified by f,\n");
  printf("[f] is either 'GG' or 'A', the former for the format in Rama-Toranria webpage, the latter for the Magma format,\n");
  printf("[L] is the preicision up to which to compute the hecke eigenvalues (a_p for p <= L and a_{p^2} for p^2 <= L),\n");
  printf("[idx] is the index of the form in the decomposition to eigenvectors.\n");
  printf("If either L or i is not supplied, only decomposes the space, and finds Hecke eigenvectors.\n");
  printf("If the flag -tests is supplied, additionally runs standard tests.\n");
  return -1;
}

int parse_matrix(const char* mat_str, int* Q_coeffs)
{
  int idx, len;
  char *original;
  char *token;

  idx = 0;
  len = strlen(mat_str);
  original = (char*)malloc((len+1)*sizeof(char));
  memcpy(original, mat_str, len+1);
  token = strsep(&original, ",");
  while(token != NULL) {
    Q_coeffs[idx++] = atoi(token);
    token = strsep(&original, ",");
  }
  free(original);
  return (idx == 15);
}
