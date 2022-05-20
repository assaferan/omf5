#ifndef CARAT_ZZ_H
#define CARAT_ZZ_H

#include <stdlib.h>		/* for atoi */
#include "typedef.h"
#include "tools.h"
#include "matrix.h"
#include "getput.h"
#include "datei.h"

extern bravais_TYP **q2z(bravais_TYP *G,
                         int *number,
                         int ADFLAG,
                         QtoZ_TYP *INZ,
                         int quiet,
                         int graph);

extern void free_QtoZ(QtoZ_TYP *inz,
                      int flag);

extern void ZZ(bravais_TYP * group,
                     matrix_TYP * gram,
                     int *divisors,
	             QtoZ_TYP *inzidenz,
                     const char *options,
                     FILE* outputfile,
                     int super_nr,
                     int konst_flag);

extern bravais_TYP **get_groups(bravais_TYP **ADGROUPS,
                         int ad_no,
                         int *number);

extern int NUMBER;		/* Abbruch nach NUMBER Zentrierungen */
extern int LEVEL;		/* Abbruch nach Iterationszahl LEVEL */

#endif /* CARAT_ZZ_H */
