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
  int test_res;
  int Q_coeffs[15];

  if (argc == 1) {
    test_res = test_61();
    // test_res |= test_69() << 1;
    return test_res;
  }
  
  if (argc == 2) {
    parse_matrix(argv[1], Q_coeffs);
    return compute_eigenvectors(Q_coeffs);
  }

  if (argc == 4) {
    parse_matrix(argv[1], Q_coeffs);
    form_idx = atoi(argv[2]);
    prec = atoi(argv[3]);
    return compute_eigenvalues_up_to(Q_coeffs, form_idx, prec);
  }
  
  if ((argc > 4) || (argc == 3))
    return -1;

}
