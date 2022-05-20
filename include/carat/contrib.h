#ifndef CARAT_CONTRIB_H
#define CARAT_CONTRIB_H

#include "typedef.h"

/*-------------------------------------------------------------*\
| FILE: torsionfree.c
\*-------------------------------------------------------------*/

int *torsionfree(bravais_TYP *R,
                 int *order_out,
                 int *number_of_conjugacy_classes);

/***************************************************************
| FILE: suche_kand.c
****************************************************************/

matrix_TYP *suche_kand (bravais_TYP *Gen_A, bravais_TYP *Gen_B);

#endif
