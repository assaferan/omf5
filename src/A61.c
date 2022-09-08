#include <carat/typedef.h>

#include "tests.h"

int main(int argc, char* argv[])
{
  int test_res;

  if (argc > 1)
    return -1;

  test_res = test_61();

  if (test_res == TRUE)
    return 0;
  else
    return -1;
  
}
