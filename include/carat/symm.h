#ifndef CARAT_SYMM_H
#define CARAT_SYMM_H

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: dsylv.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *dsylv(matrix_TYP *M);
extern int definite_test(matrix_TYP *M);

/*-------------------------------------------------------------*\
| FILE: rest_short.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *rest_short(matrix_TYP *mat, int *restvec, int rkgv,
     int zaehler, int nenner, int find_opt, int count_opt, int *anz);

/*-------------------------------------------------------------*\
| FILE: short.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *short_vectors(matrix_TYP *mat, int length,
     int lengthmin, int find_opt, int count_opt, int *anz);

/*-------------------------------------------------------------*\
| FILE: shortest.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *shortest(matrix_TYP *mat, int *min_norm);

#endif
