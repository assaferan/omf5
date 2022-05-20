#ifndef CARAT_HYPERBOLIC_H
#define CARAT_HYPERBOLIC_H

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: hyp_isom.c
\*-------------------------------------------------------------*/
extern matrix_TYP *hyperbolic_isometry(matrix_TYP *x1, matrix_TYP *x2,
     matrix_TYP *S);

/*-------------------------------------------------------------*\
| FILE: hyp_stabilizer.c
\*-------------------------------------------------------------*/
extern bravais_TYP *hyperbolic_stabilizer(matrix_TYP *x, matrix_TYP *S);

#endif
