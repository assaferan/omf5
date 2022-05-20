#ifndef CARAT_REDUCTION_H
#define CARAT_REDUCTION_H

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: mink_red.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *mink_red(matrix_TYP *G, matrix_TYP *T);

/*-------------------------------------------------------------*\
| FILE: pair_red.c 
\*-------------------------------------------------------------*/
extern void pr_red(int **G, int **T, int n);
extern matrix_TYP *pair_red(matrix_TYP *Gram, matrix_TYP *Tr);

/*-------------------------------------------------------------*\
| FILE: reduction_sort.c 
\*-------------------------------------------------------------*/
extern void reduction_sort(int **G,int **T,int n);

/*-------------------------------------------------------------*\
| FILE: short_reduce.c 
\*-------------------------------------------------------------*/
extern matrix_TYP *short_reduce(matrix_TYP *A, matrix_TYP *SV, matrix_TYP *Trf);
extern matrix_TYP *pr_short_red(matrix_TYP *A, matrix_TYP *Trf);

#endif
