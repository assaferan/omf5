#include <carat/typedef.h>

#include "tests.h"

int parse_matrix(const char* mat_str, int* Q_coeffs)
{
  int idx, len;
  char *original;
  char *token;

  idx = 0;
  len = strlen(mat_str);
  original = (char*)malloc((len+1)*sizeof(char));
  strncpy(original, mat_str, len);
  token = strsep(&original, ",");
  while(token != NULL) {
    Q_coeffs[idx++] = atoi(token);
    token = strsep(&original, ",");
  }
  free(original);
  return (idx == 15);
}

int main(int argc, char* argv[])
{
  int form_idx, prec;
  int Q_coeffs[15];
  STATUS test_res;
  char input_type[3];
  int max_args;
  int forbidden_args;

  forbidden_args = 3;
  max_args = 4;
  
  if (argc == 1) {
    test_res = test_greedy() << 1;
    test_res |= test_61();
    test_res <<= 1;
    test_res |= test_69();
    return test_res;
  }

  if (argc > 1) {
    if (strcmp(argv[argc-1], "GG") == 0) {
      strcpy(input_type, "GG");
      forbidden_args += 1;
      max_args += 1;
    }
    else if (strcmp(argv[argc-1], "A") == 0) {
      strcpy(input_type, "A");
      forbidden_args += 1;
      max_args += 1;
    }
    else
      strcpy(input_type, "A");
  }
  
  if (argc == max_args - 2) {
    parse_matrix(argv[1], Q_coeffs);
    return compute_eigenvectors(Q_coeffs, input_type);
  }

  if (argc == max_args) {
    parse_matrix(argv[1], Q_coeffs);
    form_idx = atoi(argv[2]);
    prec = atoi(argv[3]);
    return compute_eigenvalues_up_to(Q_coeffs, form_idx, prec, input_type);
  }
  
  if ((argc > max_args) || (argc == forbidden_args))
    return -1;

}
