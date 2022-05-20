#ifndef CARAT_AUTGRP_H
#define CARAT_AUTGRP_H

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE:  autgrp.c
\*-------------------------------------------------------------*/
extern bravais_TYP *autgrp( matrix_TYP **Fo, 
    int Foanz,
    matrix_TYP *SV,
    matrix_TYP **Erz,
    int Erzanz,
    int *options);

extern bravais_TYP *perfect_normal_autgrp( matrix_TYP *Fo,
    matrix_TYP *SV,
    matrix_TYP **Erz,
    int Erzanz,
    int * options,
    matrix_TYP **P, 
    int Panz,
    matrix_TYP **Pbase,
    int Pdim);

/*-------------------------------------------------------------*\
| FILE:  isometry.c
\*-------------------------------------------------------------*/
extern matrix_TYP *isometry(matrix_TYP **F1,
     matrix_TYP **F2,
     int Fanz, 
     matrix_TYP *SV1,
     matrix_TYP *SV2,
     matrix_TYP **Erz,
     int Erzanz,
     int *options);

extern matrix_TYP *perfect_normal_isometry(matrix_TYP *F1,
      matrix_TYP *F2,
      matrix_TYP *SV1,
      matrix_TYP *SV2,
      matrix_TYP **Erz,
      int Erzanz,
      int *options,
      matrix_TYP **P,
      int Panz, 
      matrix_TYP **base,
      int Pdim);

/*-------------------------------------------------------------*\
| FILE:  pr_aut.c
\*-------------------------------------------------------------*/
extern bravais_TYP *pr_aut(matrix_TYP **Fo,
      int Foanz,
      matrix_TYP **Erz,
      int Erzanz,
      int *options);

/*-------------------------------------------------------------*\
| FILE:  pr_isom.c
\*-------------------------------------------------------------*/
extern matrix_TYP *pr_isom(matrix_TYP **F1,
      matrix_TYP **F2,
      int Fanz,
      matrix_TYP **Erz,
      int Erzanz,
      int *options);

#endif
